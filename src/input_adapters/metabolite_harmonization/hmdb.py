from datetime import date
from pathlib import Path
from typing import Dict, Generator, Iterable, List, Optional, Set, Tuple, Union
import xml.etree.ElementTree as ET
import zipfile

from src.constants import DataSourceName
from src.interfaces.input_adapter import InputAdapter
from src.models.datasource_version_info import DatasourceVersionInfo
from src.models.metabolite_harmonization import (
    MetaboliteIdentifier,
    MetaboliteIdentifierMappingDetail,
    MetaboliteIdentifierMappingEdge,
    MetaboliteName,
)


HMDB_XML_MEMBER = "hmdb_metabolites.xml"
HMDB_VERSION = "5.0"
HMDB_VERSION_DATE = date(2021, 11, 17)

ID_FIELD_PREFIXES = {
    "secondary_accessions": "HMDB",
    "chebi_id": "CHEBI",
    "kegg_id": "KEGG.COMPOUND",
    "chemspider_id": "ChemSpider",
    "drugbank_id": "DRUGBANK",
    "drugbank_metabolite_id": "DRUGBANK.METABOLITE",
    "phenol_explorer_compound_id": "PhenolExplorer.COMPOUND",
    "phenol_explorer_metabolite_id": "PhenolExplorer.METABOLITE",
    "foodb_id": "FoodDB",
    "knapsack_id": "KNApSAcK",
    "pubchem_compound_id": "PUBCHEM.COMPOUND",
    "cas_registry_number": "CAS",
    "biocyc_id": "BioCyc",
    "bigg_id": "BiGG",
    "wikipedia": "Wikipedia",
    "nugowiki": "NugoWiki",
    "metagene": "MetaGene",
    "metlin_id": "METLIN",
    "het_id": "PDB.HET",
    "lipidmaps_id": "LIPIDMAPS",
}


class HmdbMetaboliteEquivalenceAdapter(InputAdapter):
    def __init__(
        self,
        data_source=None,
        hmdb_zip_file: Optional[str] = None,
        max_records: Optional[int] = None,
    ):
        if data_source is not None:
            hmdb_zip_file = str(data_source.file("hmdb_metabolites.zip"))
            self.version_info = data_source.version_info()
        else:
            self.version_info = DatasourceVersionInfo(
                version=HMDB_VERSION,
                version_date=HMDB_VERSION_DATE,
                download_date=None,
            )
        if hmdb_zip_file is None:
            raise ValueError("HmdbMetaboliteEquivalenceAdapter requires data_source or hmdb_zip_file")
        self.hmdb_zip_file = Path(hmdb_zip_file)
        self.max_records = max_records

    def get_datasource_name(self) -> DataSourceName:
        return DataSourceName.HMDB

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
        emitted_stubs: Set[str] = set()
        emitted_primary_labels: Set[str] = set()
        emitted_xref_labels: Set[Tuple[str, str]] = set()

        for record in self._iter_metabolite_records():
            primary_id = self._hmdb_id(record["accession"])
            if primary_id is None:
                continue

            if primary_id not in emitted_primary_labels:
                batch.append(self._primary_node(primary_id, record))
                emitted_primary_labels.add(primary_id)
                emitted_stubs.add(primary_id)

            for source_field, values in record["xrefs"].items():
                for value in values:
                    node_id = self._external_id(source_field, value)
                    if node_id is None or node_id == primary_id:
                        continue
                    label = record.get("name")
                    if label:
                        label_key = (node_id, label)
                        if label_key in emitted_xref_labels:
                            continue
                        emitted_xref_labels.add(label_key)
                        names = [MetaboliteName(value=label, source="HMDB", source_field="name")]
                    else:
                        if node_id in emitted_stubs:
                            continue
                        names = []
                    batch.append(MetaboliteIdentifier(id=node_id, names=names))
                    emitted_stubs.add(node_id)
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
            primary_id = self._hmdb_id(record["accession"])
            if primary_id is None:
                continue

            for source_field, values in record["xrefs"].items():
                for value in values:
                    end_id = self._external_id(source_field, value)
                    if end_id is None or end_id == primary_id:
                        continue
                    detail_key = (primary_id, end_id, "HMDB", source_field, primary_id)
                    if detail_key in emitted_details:
                        continue
                    emitted_details.add(detail_key)
                    batch.append(
                        MetaboliteIdentifierMappingEdge(
                            start_node=MetaboliteIdentifier(id=primary_id),
                            end_node=MetaboliteIdentifier(id=end_id),
                            details=[
                                MetaboliteIdentifierMappingDetail(
                                    source="HMDB",
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
        with zipfile.ZipFile(self.hmdb_zip_file) as archive:
            with archive.open(HMDB_XML_MEMBER) as handle:
                count = 0
                for _event, elem in ET.iterparse(handle, events=("end",)):
                    if self._local_name(elem.tag) != "metabolite":
                        continue
                    yield self._parse_metabolite(elem)
                    count += 1
                    elem.clear()
                    if self.max_records is not None and count >= self.max_records:
                        break

    @classmethod
    def _parse_metabolite(cls, elem) -> Dict:
        record = {
            "accession": None,
            "name": None,
            "status": None,
            "synonyms": [],
            "xrefs": {field: [] for field in ID_FIELD_PREFIXES},
        }
        for child in list(elem):
            tag = cls._local_name(child.tag)
            if tag == "accession":
                record["accession"] = cls._clean_text(child.text)
            elif tag == "name":
                record["name"] = cls._clean_text(child.text)
            elif tag == "status":
                record["status"] = cls._clean_text(child.text)
            elif tag == "synonyms":
                record["synonyms"] = cls._unique_clean_values(
                    cls._clean_text(synonym.text)
                    for synonym in list(child)
                    if cls._local_name(synonym.tag) == "synonym"
                )
            elif tag == "secondary_accessions":
                record["xrefs"]["secondary_accessions"] = cls._unique_clean_values(
                    cls._clean_text(accession.text)
                    for accession in list(child)
                    if cls._local_name(accession.tag) == "accession"
                )
            elif tag in ID_FIELD_PREFIXES:
                value = cls._clean_text(child.text)
                if value:
                    record["xrefs"][tag] = [value]
        return record

    @staticmethod
    def _primary_node(primary_id: str, record: Dict) -> MetaboliteIdentifier:
        names = []
        if record.get("name"):
            names.append(
                MetaboliteName(
                    value=record["name"],
                    source="HMDB",
                    source_field="name",
                )
            )
        synonyms = [
            MetaboliteName(
                value=synonym,
                source="HMDB",
                source_field="synonym",
            )
            for synonym in record.get("synonyms", [])
        ]
        return MetaboliteIdentifier(id=primary_id, names=names, synonyms=synonyms,
                                    hmdb_status=record.get("status"))

    @classmethod
    def _external_id(cls, source_field: str, value: Optional[str]) -> Optional[str]:
        if source_field == "secondary_accessions":
            return cls._hmdb_id(value)
        prefix = ID_FIELD_PREFIXES.get(source_field)
        value = cls._clean_text(value)
        if not prefix or not value:
            return None
        if prefix == "CHEBI" and value.upper().startswith("CHEBI:"):
            value = value.split(":", 1)[1].strip()
        elif prefix == "PUBCHEM.COMPOUND":
            value = value.removeprefix("CID").strip()
        return f"{prefix}:{value}"

    @classmethod
    def _hmdb_id(cls, value: Optional[str]) -> Optional[str]:
        value = cls._clean_text(value)
        if not value:
            return None
        if value.upper().startswith("HMDB:"):
            return f"HMDB:{value.split(':', 1)[1].strip()}"
        return f"HMDB:{value}"

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

    @staticmethod
    def _local_name(tag: str) -> str:
        return tag.split("}", 1)[-1]
