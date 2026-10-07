from pathlib import Path
from dataclasses import dataclass
from typing import Dict, Generator, Iterable, List, Optional
import csv
import gzip
import zipfile
import xml.etree.ElementTree as ET

from src.constants import DataSourceName
from src.interfaces.input_adapter import InputAdapter
from src.models.datasource_version_info import DatasourceVersionInfo
from src.models.metabolite_harmonization import (
    MetaboliteChemProps,
    MetaboliteIdentifier,
    MetaboliteName,
    MetaboliteStructureComponent,
)
from src.shared.chebi_mass import validated_chebi_mass_values
from src.shared.metabolite_structure_chemistry import (
    DERIVED_INCHI_KEY_METHOD,
    calculate_smiles_chemistry,
    derive_inchi_key_from_smiles,
)


HMDB_STRUCTURES_SDF_MEMBER = "structures.sdf"
HMDB_METABOLITES_XML_MEMBER = "hmdb_metabolites.xml"
LIPIDMAPS_SDF_MEMBER = "structures.sdf"
SDF_MOL_BLOCK_INPUT_FIELD = "sdf_mol_block"


@dataclass(frozen=True)
class _SdfRecord:
    tags: Dict[str, str]
    mol_block: str

    def get(self, key: str, default=None):
        return self.tags.get(key, default)


class _ChemPropsAdapter(InputAdapter):
    def get_all(self) -> Generator[List[MetaboliteIdentifier], None, None]:
        batch: List[MetaboliteIdentifier] = []
        for node in self._iter_nodes():
            batch.append(node)
            if len(batch) >= self.batch_size:
                yield batch
                batch = []
        if batch:
            yield batch

    def _iter_nodes(self) -> Iterable[MetaboliteIdentifier]:
        raise NotImplementedError


class HmdbMetaboliteChemPropsAdapter(_ChemPropsAdapter):
    def __init__(
        self,
        data_source=None,
        metabolites_data_source=None,
        structures_zip_file: Optional[str] = None,
        metabolites_zip_file: Optional[str] = None,
        max_records: Optional[int] = None,
    ):
        if data_source is not None:
            structures_zip_file = str(data_source.file("structures.zip"))
            self.version_info = data_source.version_info()
        else:
            self.version_info = DatasourceVersionInfo(version=None, version_date=None, download_date=None)
        if structures_zip_file is None:
            raise ValueError("HmdbMetaboliteChemPropsAdapter requires data_source or structures_zip_file")
        if metabolites_data_source is not None:
            metabolites_zip_file = str(
                metabolites_data_source.file("hmdb_metabolites.zip")
            )
        self.structures_zip_file = Path(structures_zip_file)
        self.metabolites_zip_file = (
            Path(metabolites_zip_file) if metabolites_zip_file is not None else None
        )
        self.max_records = max_records

    def get_datasource_name(self) -> DataSourceName:
        return DataSourceName.HMDB

    def get_version(self) -> DatasourceVersionInfo:
        return self.version_info

    def _iter_nodes(self) -> Iterable[MetaboliteIdentifier]:
        count = 0
        sdf_source_ids: set[str] = set()
        with zipfile.ZipFile(self.structures_zip_file) as archive:
            with archive.open(HMDB_STRUCTURES_SDF_MEMBER) as handle:
                for record in _iter_sdf_tag_records(handle):
                    source_id = _prefixed_id("HMDB", record.get("DATABASE_ID"))
                    if source_id is None:
                        continue
                    if source_id in sdf_source_ids:
                        continue
                    sdf_source_ids.add(source_id)
                    inchi_key = _clean_text(record.get("INCHI_KEY"))
                    yield MetaboliteIdentifier(
                        id=source_id,
                        chem_props=[
                            _metabolite_chem_props(
                                source="HMDB",
                                source_id=source_id,
                                iso_smiles=_clean_text(record.get("SMILES")),
                                inchi_key_prefix=_inchi_key_prefix(inchi_key),
                                inchi_key=inchi_key,
                                inchi=_clean_text(record.get("INCHI_IDENTIFIER")),
                                mw=_clean_text(record.get("MOLECULAR_WEIGHT")),
                                monoisotopic_mass=_clean_text(record.get("EXACT_MASS")),
                                common_name=_clean_text(record.get("GENERIC_NAME")),
                                molecular_formula=_clean_text(record.get("FORMULA")),
                            )
                        ],
                    )
                    count += 1
                    if self.max_records is not None and count >= self.max_records:
                        return
        if (
            self.metabolites_zip_file is None
            or (self.max_records is not None and count >= self.max_records)
        ):
            return

        xml_source_ids: set[str] = set()
        with zipfile.ZipFile(self.metabolites_zip_file) as archive:
            with archive.open(HMDB_METABOLITES_XML_MEMBER) as handle:
                for _event, elem in ET.iterparse(handle, events=("end",)):
                    if _xml_local_name(elem.tag) != "metabolite":
                        continue
                    values = {
                        _xml_local_name(child.tag): _clean_text(child.text)
                        for child in list(elem)
                    }
                    elem.clear()
                    source_id = _prefixed_id("HMDB", values.get("accession"))
                    if (
                        source_id is None
                        or source_id in sdf_source_ids
                        or source_id in xml_source_ids
                    ):
                        continue
                    xml_source_ids.add(source_id)
                    props = _hmdb_xml_chem_props(source_id, values)
                    if props is None:
                        continue
                    yield MetaboliteIdentifier(id=source_id, chem_props=[props])
                    count += 1
                    if self.max_records is not None and count >= self.max_records:
                        return


class ChebiMetaboliteChemPropsAdapter(_ChemPropsAdapter):
    def __init__(
        self,
        data_source=None,
        chebi_sdf_file: Optional[str] = None,
        max_records: Optional[int] = None,
    ):
        if data_source is not None:
            chebi_sdf_file = str(data_source.file())
            self.version_info = data_source.version_info()
        else:
            self.version_info = DatasourceVersionInfo(version=None, version_date=None, download_date=None)
        if chebi_sdf_file is None:
            raise ValueError("ChebiMetaboliteChemPropsAdapter requires data_source or chebi_sdf_file")
        self.chebi_sdf_file = Path(chebi_sdf_file)
        self.max_records = max_records

    def get_datasource_name(self) -> DataSourceName:
        return DataSourceName.ChEBI

    def get_version(self) -> DatasourceVersionInfo:
        return self.version_info

    def _iter_nodes(self) -> Iterable[MetaboliteIdentifier]:
        count = 0
        with _open_text_or_gzip(self.chebi_sdf_file) as handle:
            for record in _iter_sdf_tag_records(handle):
                source_id = _prefixed_id("CHEBI", record.get("ChEBI ID"))
                if source_id is None:
                    continue
                name_field = next((field for field in ("ChEBI NAME", "ChEBI Name")
                                   if _clean_text(record.get(field))), None)
                chebi_name = _clean_text(record.get(name_field)) if name_field else None
                inchi_key = _record_text(record, "INCHIKEY", "InChIKey")
                smiles = _clean_text(record.get("SMILES"))
                molecular_formula = _record_text(record, "FORMULA", "Formulae")
                mass, monoisotopic_mass = validated_chebi_mass_values(
                    molecular_formula,
                    _clean_text(record.get("MASS")),
                    _record_text(record, "MONOISOTOPIC_MASS", "Monoisotopic Mass"),
                    smiles,
                )
                yield MetaboliteIdentifier(
                    id=source_id,
                    names=([MetaboliteName(value=chebi_name, source="ChEBI",
                                           source_field=name_field)] if chebi_name else []),
                    chem_props=[
                        _metabolite_chem_props(
                            source="ChEBI",
                            source_id=source_id,
                            iso_smiles=smiles,
                            inchi_key_prefix=_inchi_key_prefix(inchi_key),
                            inchi_key=inchi_key,
                            inchi=_record_text(record, "INCHI", "InChI"),
                            mw=mass,
                            monoisotopic_mass=monoisotopic_mass,
                            common_name=chebi_name,
                            molecular_formula=molecular_formula,
                        )
                    ],
                )
                count += 1
                if self.max_records is not None and count >= self.max_records:
                    return


class LipidMapsMetaboliteChemPropsAdapter(_ChemPropsAdapter):
    def __init__(
        self,
        data_source=None,
        sdf_zip_file: Optional[str] = None,
        max_records: Optional[int] = None,
    ):
        if data_source is not None:
            sdf_zip_file = str(data_source.file("LMSD.sdf.zip"))
            self.version_info = data_source.version_info()
        else:
            self.version_info = DatasourceVersionInfo(version=None, version_date=None, download_date=None)
        if sdf_zip_file is None:
            raise ValueError("LipidMapsMetaboliteChemPropsAdapter requires data_source or sdf_zip_file")
        self.sdf_zip_file = Path(sdf_zip_file)
        self.max_records = max_records

    def get_datasource_name(self) -> DataSourceName:
        return DataSourceName.LipidMaps

    def get_version(self) -> DatasourceVersionInfo:
        return self.version_info

    def _iter_nodes(self) -> Iterable[MetaboliteIdentifier]:
        count = 0
        with zipfile.ZipFile(self.sdf_zip_file) as archive:
            with archive.open(LIPIDMAPS_SDF_MEMBER) as handle:
                for record in _iter_sdf_tag_records(handle):
                    source_id = _prefixed_id("LIPIDMAPS", record.get("LM_ID"))
                    if source_id is None:
                        continue
                    inchi_key = _clean_text(record.get("INCHI_KEY"))
                    reported_smiles = _clean_text(record.get("SMILES"))
                    mol_smiles = None
                    mol_inchi_key_derivation = None
                    if reported_smiles is None:
                        mol_smiles, mol_inchi_key_derivation = _derive_lipidmaps_mol_chemistry(
                            record.mol_block
                        )
                    yield MetaboliteIdentifier(
                        id=source_id,
                        chem_props=[
                            _metabolite_chem_props(
                                precomputed_inchi_key_derivation=mol_inchi_key_derivation,
                                source="LipidMaps",
                                source_id=source_id,
                                iso_smiles=reported_smiles,
                                canonical_smiles=mol_smiles,
                                inchi_key_prefix=_inchi_key_prefix(inchi_key),
                                inchi_key=inchi_key,
                                inchi=_clean_text(record.get("INCHI")),
                                mw=_clean_text(record.get("MASS")),
                                monoisotopic_mass=_clean_text(record.get("EXACT_MASS")),
                                common_name=_first_clean(
                                    record.get("COMMON_NAME"),
                                    record.get("NAME"),
                                    record.get("SYSTEMATIC_NAME"),
                                    record.get("ABBREVIATION"),
                                ),
                                molecular_formula=_clean_text(record.get("FORMULA")),
                            )
                        ],
                    )
                    count += 1
                    if self.max_records is not None and count >= self.max_records:
                        return


class PubchemMetaboliteChemPropsAdapter(_ChemPropsAdapter):
    def __init__(
        self,
        data_source=None,
        molecular_info_file: Optional[str] = None,
        max_records: Optional[int] = None,
    ):
        if data_source is not None:
            molecular_info_file = str(data_source.file("cid_molecular_info.tsv"))
            self.version_info = data_source.version_info()
        else:
            self.version_info = DatasourceVersionInfo(version=None, version_date=None, download_date=None)
        if molecular_info_file is None:
            raise ValueError("PubchemMetaboliteChemPropsAdapter requires data_source or molecular_info_file")
        self.molecular_info_file = Path(molecular_info_file)
        self.max_records = max_records

    def get_datasource_name(self) -> DataSourceName:
        return DataSourceName.PubChem

    def get_version(self) -> DatasourceVersionInfo:
        return self.version_info

    def _iter_nodes(self) -> Iterable[MetaboliteIdentifier]:
        count = 0
        with self.molecular_info_file.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                source_id = _clean_text(row.get("pubchem_id"))
                if source_id is None:
                    continue
                isomeric_smiles = _clean_text(row.get("isomeric_smiles"))
                yield MetaboliteIdentifier(
                    id=source_id,
                    names=([MetaboliteName(value=title, source="PubChem", source_field="title")]
                           if (title := _clean_text(row.get("title"))) else []),
                    chem_props=[
                        _metabolite_chem_props(
                            source="PubChem",
                            source_id=source_id,
                            iso_smiles=isomeric_smiles,
                            canonical_smiles=_clean_text(row.get("canonical_smiles")),
                            isomeric_smiles=isomeric_smiles,
                            inchi_key_prefix=_clean_text(row.get("inchi_key_prefix")),
                            inchi_key=_clean_text(row.get("inchikey")),
                            inchi=_clean_text(row.get("inchi")),
                            mw=_clean_text(row.get("molecular_weight")),
                            monoisotopic_mass=_clean_text(row.get("monoisotopic_mass")),
                            iupac_name=_clean_text(row.get("iupac_name")),
                            molecular_formula=_clean_text(row.get("molecular_formula")),
                        )
                    ],
                )
                count += 1
                if self.max_records is not None and count >= self.max_records:
                    return


def _iter_sdf_tag_records(handle) -> Iterable[_SdfRecord]:
    lines: List[str] = []
    for raw_line in handle:
        if isinstance(raw_line, bytes):
            line = raw_line.decode("utf-8", "replace").rstrip("\r\n")
        else:
            line = raw_line.rstrip("\r\n")
        if line == "$$$$":
            yield _parse_sdf_record(lines)
            lines = []
        else:
            lines.append(line)
    if lines:
        yield _parse_sdf_record(lines)


def _parse_sdf_record(lines: List[str]) -> _SdfRecord:
    first_tag_index = next(
        (index for index, line in enumerate(lines) if line.startswith("> <") and line.endswith(">")),
        len(lines),
    )
    mol_lines = lines[:first_tag_index]
    while mol_lines and not mol_lines[-1]:
        mol_lines.pop()
    return _SdfRecord(
        tags=_parse_sdf_tags(lines[first_tag_index:]),
        mol_block=("\n".join(mol_lines) + "\n") if mol_lines else "",
    )


def _parse_sdf_tags(lines: List[str]) -> Dict[str, str]:
    tags: Dict[str, List[str]] = {}
    current_tag = None
    for line in lines:
        if line.startswith("> <") and line.endswith(">"):
            current_tag = line[3:-1]
            tags.setdefault(current_tag, [])
            continue
        if current_tag is not None:
            tags[current_tag].append(line)
    return {
        tag: "\n".join(value_lines).strip()
        for tag, value_lines in tags.items()
        if "\n".join(value_lines).strip()
    }


def _open_text_or_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _clean_text(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    value = str(value).strip()
    return value or None


def _xml_local_name(tag: str) -> str:
    return tag.split("}", 1)[-1]


def _hmdb_xml_chem_props(
    source_id: str,
    values: Dict[str, Optional[str]],
) -> Optional[MetaboliteChemProps]:
    smiles = _clean_text(values.get("smiles"))
    inchi_key = _clean_text(values.get("inchikey"))
    inchi = _clean_text(values.get("inchi"))
    mw = _clean_text(values.get("average_molecular_weight"))
    monoisotopic_mass = _clean_text(values.get("monisotopic_molecular_weight"))
    molecular_formula = _clean_text(values.get("chemical_formula"))
    if not any((smiles, inchi_key, inchi, mw, monoisotopic_mass, molecular_formula)):
        return None
    return _metabolite_chem_props(
        source="HMDB",
        source_id=source_id,
        iso_smiles=smiles,
        inchi_key_prefix=_inchi_key_prefix(inchi_key),
        inchi_key=inchi_key,
        inchi=inchi,
        mw=mw,
        monoisotopic_mass=monoisotopic_mass,
        common_name=_clean_text(values.get("name")),
        molecular_formula=molecular_formula,
    )


def _record_text(record: Dict[str, str], *field_names: str) -> Optional[str]:
    return _first_clean(*(record.get(field_name) for field_name in field_names))


def _first_clean(*values: Optional[str]) -> Optional[str]:
    for value in values:
        value = _clean_text(value)
        if value:
            return value
    return None


def _prefixed_id(prefix: str, value: Optional[str]) -> Optional[str]:
    value = _clean_text(value)
    if not value:
        return None
    if prefix == "CHEBI" and value.upper().startswith("CHEBI:"):
        value = value.split(":", 1)[1].strip()
    if prefix == "LIPIDMAPS" and value.upper().startswith("LIPIDMAPS:"):
        value = value.split(":", 1)[1].strip()
    if prefix == "HMDB" and value.upper().startswith("HMDB:"):
        value = value.split(":", 1)[1].strip()
    return f"{prefix}:{value}"


def _inchi_key_prefix(inchi_key: Optional[str]) -> Optional[str]:
    inchi_key = _clean_text(inchi_key)
    if not inchi_key:
        return None
    return inchi_key.split("-", 1)[0]


def _metabolite_chem_props(
    *,
    precomputed_inchi_key_derivation: Optional[dict] = None,
    **values,
) -> MetaboliteChemProps:
    structure_candidates = (
        ("iso_smiles", _clean_text(values.get("iso_smiles"))),
        ("isomeric_smiles", _clean_text(values.get("isomeric_smiles"))),
        ("canonical_smiles", _clean_text(values.get("canonical_smiles"))),
    )
    structure_input_field, structure_smiles = next(
        ((field_name, value) for field_name, value in structure_candidates if value),
        (None, None),
    )
    structure_values = calculate_smiles_chemistry(structure_smiles, structure_input_field)
    component_values = structure_values.pop("structure_components", [])
    derived_values = (
        precomputed_inchi_key_derivation
        if precomputed_inchi_key_derivation is not None
        else _derive_inchi_key_from_smiles(
            iso_smiles=values.get("iso_smiles"),
            isomeric_smiles=values.get("isomeric_smiles"),
            canonical_smiles=values.get("canonical_smiles"),
        )
    )
    return MetaboliteChemProps(
        **values,
        **derived_values,
        **structure_values,
        structure_components=[
            MetaboliteStructureComponent(**component)
            for component in component_values
        ],
    )


def _derive_lipidmaps_mol_chemistry(mol_block: str) -> tuple[Optional[str], dict]:
    if not mol_block.strip():
        return None, {}

    try:
        import rdkit
        from rdkit import Chem, rdBase
        from rdkit.Chem import inchi
    except ImportError as exc:
        raise RuntimeError(
            "RDKit is required to derive chemistry from LipidMaps mol blocks"
        ) from exc

    derivation_metadata = {
        "derived_inchi_key_input_field": SDF_MOL_BLOCK_INPUT_FIELD,
        "derived_inchi_key_method": DERIVED_INCHI_KEY_METHOD,
        "derived_inchi_key_method_version": rdkit.__version__,
    }
    try:
        with rdBase.BlockLogs():
            molecule = Chem.MolFromMolBlock(
                mol_block,
                sanitize=True,
                removeHs=False,
                strictParsing=False,
            )
    except Exception as exc:
        return None, {
            **derivation_metadata,
            "derived_inchi_key_error": f"mol_parse_error:{type(exc).__name__}",
        }
    if molecule is None:
        return None, {
            **derivation_metadata,
            "derived_inchi_key_error": "mol_parse_failed",
        }
    if molecule.GetNumAtoms() == 0:
        return None, {}

    try:
        canonical_smiles = Chem.MolToSmiles(
            molecule,
            canonical=True,
            isomericSmiles=True,
        )
    except Exception as exc:
        return None, {
            **derivation_metadata,
            "derived_inchi_key_error": f"smiles_generation_error:{type(exc).__name__}",
        }

    has_query_features = any(
        atom.GetAtomicNum() == 0 or atom.HasQuery()
        for atom in molecule.GetAtoms()
    ) or any(bond.HasQuery() for bond in molecule.GetBonds())
    if has_query_features:
        return canonical_smiles, {
            **derivation_metadata,
            "derived_inchi_key_error": "query_structure_not_supported",
        }

    try:
        with rdBase.BlockLogs():
            derived_inchi_key = _clean_text(inchi.MolToInchiKey(molecule))
    except Exception as exc:
        return canonical_smiles, {
            **derivation_metadata,
            "derived_inchi_key_error": f"inchi_key_generation_error:{type(exc).__name__}",
        }
    if derived_inchi_key is None:
        return canonical_smiles, {
            **derivation_metadata,
            "derived_inchi_key_error": "inchi_key_generation_failed",
        }
    return canonical_smiles, {
        **derivation_metadata,
        "derived_inchi_key_prefix": _inchi_key_prefix(derived_inchi_key),
        "derived_inchi_key": derived_inchi_key,
    }


def _derive_inchi_key_from_smiles(
    iso_smiles: Optional[str] = None,
    isomeric_smiles: Optional[str] = None,
    canonical_smiles: Optional[str] = None,
) -> dict:
    return derive_inchi_key_from_smiles(
        iso_smiles=iso_smiles,
        isomeric_smiles=isomeric_smiles,
        canonical_smiles=canonical_smiles,
    )
