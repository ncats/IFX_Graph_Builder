from pathlib import Path
from typing import Dict, Generator, Iterable, List, Optional, Set, Tuple, Union
import csv
import re

from src.constants import DataSourceName
from src.interfaces.input_adapter import InputAdapter
from src.models.datasource_version_info import DatasourceVersionInfo
from src.models.metabolite_harmonization import (
    MetaboliteChemProps,
    MetaboliteIdentifier,
    MetaboliteIdentifierMappingDetail,
    MetaboliteIdentifierMappingEdge,
    MetaboliteName,
)


ID_FIELD_PREFIXES = {
    "pubchem_cid": "PUBCHEM.COMPOUND",
    "chebi_id": "CHEBI",
    "hmdb_id": "HMDB",
    "lipidmaps_id": "LIPIDMAPS",
    "kegg_id": "KEGG.COMPOUND",
}


class RefMetMetaboliteEquivalenceAdapter(InputAdapter):
    def __init__(
        self,
        data_source=None,
        refmet_csv_file: Optional[str] = None,
        max_records: Optional[int] = None,
    ):
        if data_source is not None:
            refmet_csv_file = str(data_source.file("refmet.csv"))
            self.version_info = data_source.version_info()
        else:
            self.version_info = DatasourceVersionInfo(
                version=None,
                version_date=None,
                download_date=None,
            )
        if refmet_csv_file is None:
            raise ValueError("RefMetMetaboliteEquivalenceAdapter requires data_source or refmet_csv_file")
        self.refmet_csv_file = Path(refmet_csv_file)
        self.max_records = max_records

    def get_datasource_name(self) -> DataSourceName:
        return DataSourceName.RefMet

    def get_version(self) -> DatasourceVersionInfo:
        return self.version_info

    def get_all(self) -> Generator[
        List[Union[MetaboliteIdentifier, MetaboliteIdentifierMappingEdge]],
        None,
        None,
    ]:
        yield from self._iter_node_batches()
        yield from self._iter_edge_batches()

    def _iter_node_batches(self) -> Generator[List[MetaboliteIdentifier], None, None]:
        batch: List[MetaboliteIdentifier] = []
        emitted_ids: Set[str] = set()
        emitted_labeled_records: Set[Tuple[str, Tuple[str, ...]]] = set()
        emitted_xref_labels: Set[Tuple[str, str]] = set()

        for record in self._iter_metabolite_records():
            primary_id = self._refmet_id(record["refmet_id"])
            if primary_id is None:
                continue

            name_key = tuple(name.value for name in self._names(record))
            label_key = (primary_id, name_key)
            if label_key not in emitted_labeled_records:
                batch.append(self._primary_node(primary_id, record))
                emitted_labeled_records.add(label_key)
                emitted_ids.add(primary_id)

            for source_field, values in record["xrefs"].items():
                for value in values:
                    node_id = self._external_id(source_field, value)
                    if node_id is None or node_id == primary_id:
                        continue
                    label = record.get("refmet_name")
                    if label:
                        label_key = (node_id, label)
                        if label_key in emitted_xref_labels:
                            continue
                        emitted_xref_labels.add(label_key)
                        names = [MetaboliteName(value=label, source="RefMet", source_field="refmet_name")]
                    else:
                        if node_id in emitted_ids:
                            continue
                        names = []
                    batch.append(MetaboliteIdentifier(id=node_id, names=names))
                    emitted_ids.add(node_id)
                    if len(batch) >= self.batch_size:
                        yield batch
                        batch = []

            if len(batch) >= self.batch_size:
                yield batch
                batch = []

        if batch:
            yield batch

    def _iter_edge_batches(self) -> Generator[List[MetaboliteIdentifierMappingEdge], None, None]:
        batch: List[MetaboliteIdentifierMappingEdge] = []
        emitted_details: Set[Tuple[str, str, str, str, str]] = set()

        for record in self._iter_metabolite_records():
            primary_id = self._refmet_id(record["refmet_id"])
            if primary_id is None:
                continue

            for source_field, values in record["xrefs"].items():
                for value in values:
                    end_id = self._external_id(source_field, value)
                    if end_id is None or end_id == primary_id:
                        continue
                    detail_key = (primary_id, end_id, "RefMet", source_field, primary_id)
                    if detail_key in emitted_details:
                        continue
                    emitted_details.add(detail_key)
                    batch.append(
                        MetaboliteIdentifierMappingEdge(
                            start_node=MetaboliteIdentifier(id=primary_id),
                            end_node=MetaboliteIdentifier(id=end_id),
                            details=[
                                MetaboliteIdentifierMappingDetail(
                                    source="RefMet",
                                    source_field=source_field,
                                    source_id=primary_id,
                                )
                            ],
                        )
                    )
                    if len(batch) >= self.batch_size:
                        yield batch
                        batch = []

        if batch:
            yield batch

    def _iter_metabolite_records(self) -> Iterable[Dict]:
        with self.refmet_csv_file.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            count = 0
            for row in reader:
                normalized = {
                    (key or "").strip(): self._clean_text(value)
                    for key, value in row.items()
                }
                if not normalized.get("refmet_id"):
                    continue
                yield {
                    "refmet_id": normalized.get("refmet_id"),
                    "refmet_name": normalized.get("refmet_name"),
                    "formula": normalized.get("formula"),
                    "exactmass": normalized.get("exactmass"),
                    "inchi_key": normalized.get("inchi_key"),
                    "xrefs": {
                        field: self._unique_clean_values([normalized.get(field)])
                        for field in ID_FIELD_PREFIXES
                    },
                }
                count += 1
                if self.max_records is not None and count >= self.max_records:
                    return

    @classmethod
    def _primary_node(cls, primary_id: str, record: Dict) -> MetaboliteIdentifier:
        return MetaboliteIdentifier(
            id=primary_id,
            names=cls._names(record),
            chem_props=cls._chem_props(primary_id, record),
        )

    @classmethod
    def _chem_props(cls, primary_id: str, record: Dict) -> List[MetaboliteChemProps]:
        formula = cls._clean_text(record.get("formula"))
        exactmass = cls._clean_text(record.get("exactmass"))
        inchi_key = cls._clean_text(record.get("inchi_key"))
        if not any((formula, exactmass, inchi_key)):
            return []
        return [
            MetaboliteChemProps(
                source="RefMet",
                source_id=primary_id,
                inchi_key_prefix=(inchi_key.split("-", 1)[0] if inchi_key else None),
                inchi_key=inchi_key,
                monoisotopic_mass=exactmass,
                molecular_formula=formula,
            )
        ]

    @staticmethod
    def _names(record: Dict) -> List[MetaboliteName]:
        value = RefMetMetaboliteEquivalenceAdapter._clean_text(record.get("refmet_name"))
        if not value:
            return []
        return [
            MetaboliteName(
                value=value,
                source="RefMet",
                source_field="refmet_name",
            )
        ]

    @classmethod
    def _external_id(cls, source_field: str, value: Optional[str]) -> Optional[str]:
        prefix = ID_FIELD_PREFIXES.get(source_field)
        value = cls._clean_text(value)
        if not prefix or not value:
            return None
        if prefix == "CHEBI" and value.upper().startswith("CHEBI:"):
            value = value.split(":", 1)[1].strip()
        if prefix == "PUBCHEM.COMPOUND":
            value = re.sub(r"^CID", "", value, flags=re.IGNORECASE)
        return f"{prefix}:{value}"

    @classmethod
    def _refmet_id(cls, value: Optional[str]) -> Optional[str]:
        value = cls._clean_text(value)
        if not value:
            return None
        if value.upper().startswith("REFMET:"):
            return f"REFMET:{value.split(':', 1)[1].strip()}"
        return f"REFMET:{value}"

    @staticmethod
    def _clean_text(value: Optional[str]) -> Optional[str]:
        if value is None:
            return None
        value = value.strip()
        return value or None

    @classmethod
    def _unique_clean_values(cls, values: Iterable[Optional[str]]) -> List[str]:
        result = []
        seen = set()
        for value in values:
            value = cls._clean_text(value)
            if not value or value in seen:
                continue
            seen.add(value)
            result.append(value)
        return result
