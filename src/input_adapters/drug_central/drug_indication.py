"""DrugCentral drug indications joined to human protein interactions."""

from collections import OrderedDict, defaultdict
from collections.abc import Generator
import hashlib

from src.constants import DataSourceName, Prefix
from src.input_adapters.drug_central.drug_node import DrugCentralAdapter
from src.interfaces.input_adapter import InputAdapter
from src.models.datasource_version_info import DatasourceVersionInfo
from src.models.disease import Disease, DiseaseAssociationDetail, ProteinDiseaseEdge
from src.models.node import EquivalentId, Node, Relationship
from src.models.protein import Protein


class DrugCentralIndicationAdapter(InputAdapter, DrugCentralAdapter):
    batch_size = 10000

    def __init__(self, data_source):
        DrugCentralAdapter.__init__(self, data_source)

    def get_datasource_name(self) -> DataSourceName:
        return DataSourceName.DrugCentral

    def get_version(self) -> DatasourceVersionInfo:
        return self.version_info

    def get_all(self) -> Generator[list[Node | Relationship], None, None]:
        drug_names: dict[str, str] = {}
        for row in self.rows("drugcentral_structures.smiles.tsv", ("ID", "DRUG_NAME")):
            struct_id = row["ID"].strip()
            if not struct_id:
                raise ValueError("DrugCentral structures export contains a blank ID")
            drug_names[struct_id] = row["DRUG_NAME"]

        accessions_by_drug: dict[str, set[str]] = defaultdict(set)
        for row in self.rows(
            "drugcentral_drug_target_interaction.tsv.gz",
            ("STRUCT_ID", "ACCESSION", "ORGANISM"),
        ):
            if row["ORGANISM"] != "Homo sapiens":
                continue
            struct_id = row["STRUCT_ID"].strip()
            accession = row["ACCESSION"].strip()
            if not struct_id or not accession:
                raise ValueError(
                    "DrugCentral interaction export has blank STRUCT_ID or ACCESSION"
                )
            if struct_id not in drug_names:
                raise ValueError(
                    f"DrugCentral interaction STRUCT_ID {struct_id} has no structure"
                )
            accessions_by_drug[struct_id].add(accession)

        diseases_by_id: OrderedDict[str, Disease] = OrderedDict()
        edge_by_key: OrderedDict[tuple[str, str], ProteinDiseaseEdge] = OrderedDict()
        seen_detail_keys: set[tuple[str, str, str, str | None, str | None]] = set()
        seen_join_rows: set[tuple[str, ...]] = set()

        for row in self.rows(
            "drugcentral_indications.tsv",
            (
                "STRUCT_ID", "RELATIONSHIP_NAME", "CONCEPT_NAME", "UMLS_CUI",
                "SNOMED_CONCEPTID", "DOID",
            ),
        ):
            if row["RELATIONSHIP_NAME"] != "indication":
                continue
            struct_id = row["STRUCT_ID"].strip()
            if not struct_id:
                raise ValueError("DrugCentral indication export contains a blank STRUCT_ID")
            accessions = accessions_by_drug.get(struct_id)
            if not accessions:
                continue
            drug_name = drug_names[struct_id].strip() or None
            disease_name = row["CONCEPT_NAME"].strip()
            if not disease_name:
                continue
            umls_cui = row["UMLS_CUI"].strip() or None
            disease_id = self._disease_id(disease_name, umls_cui)
            diseases_by_id.setdefault(
                disease_id, Disease(id=disease_id, name=disease_name)
            )
            snomed_id = (
                EquivalentId(
                    id=row["SNOMED_CONCEPTID"].strip(), type=Prefix.SNOMEDCT
                ).id_str()
                if row["SNOMED_CONCEPTID"].strip() else None
            )
            doid = row["DOID"].strip() or None

            for accession_text in sorted(accessions):
                join_key = (
                    struct_id, drug_name or "", row["CONCEPT_NAME"],
                    row["UMLS_CUI"], row["SNOMED_CONCEPTID"], row["DOID"],
                    accession_text,
                )
                if join_key in seen_join_rows:
                    continue
                seen_join_rows.add(join_key)
                for accession in self._split_accessions(accession_text):
                    protein_id = EquivalentId(
                        id=accession, type=Prefix.UniProtKB
                    ).id_str()
                    edge_key = (protein_id, disease_id)
                    detail_key = (
                        protein_id, disease_id, drug_name or "", snomed_id, doid
                    )
                    if detail_key in seen_detail_keys:
                        continue
                    seen_detail_keys.add(detail_key)
                    detail = DiseaseAssociationDetail(
                        source="DrugCentral Indication",
                        source_id=(
                            EquivalentId(id=umls_cui, type=Prefix.UMLS).id_str()
                            if umls_cui else None
                        ),
                        drug_name=drug_name,
                        snomed_id=snomed_id,
                        doid=doid,
                    )
                    if edge_key not in edge_by_key:
                        edge_by_key[edge_key] = ProteinDiseaseEdge(
                            start_node=Protein(id=protein_id),
                            end_node=diseases_by_id[disease_id],
                            details=[detail],
                        )
                    else:
                        edge_by_key[edge_key].details.append(detail)

        disease_values = list(diseases_by_id.values())
        for index in range(0, len(disease_values), self.batch_size):
            yield disease_values[index:index + self.batch_size]
        edge_values = list(edge_by_key.values())
        for index in range(0, len(edge_values), self.batch_size):
            yield edge_values[index:index + self.batch_size]

    @staticmethod
    def _split_accessions(raw_accessions: str) -> list[str]:
        return [
            token.strip()
            for token in raw_accessions.split("|")
            if token and token.strip()
        ]

    @staticmethod
    def _disease_id(disease_name: str, umls_cui: str | None) -> str:
        if umls_cui:
            return EquivalentId(id=umls_cui, type=Prefix.UMLS).id_str()
        digest = hashlib.sha1(disease_name.strip().lower().encode("utf-8")).hexdigest()[:16]
        return f"DrugCentral:INDICATION:{digest}"
