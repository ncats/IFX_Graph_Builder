from typing import Generator, List, Union

from src.input_adapters.pharos_arango.tcrd.protein import PharosArangoAdapter
from src.models.datasource_version_info import DataSourceDetails
from src.models.disease import Disease, DiseaseParentEdge, ProteinDiseaseEdge, DiseaseAssociationDetail, TINXImportanceEdge
from src.models.protein import Protein


def disease_query() -> str:
    return """FOR d IN `Disease`
    RETURN d
    """


def associated_disease_ids_query(collection_name: str) -> str:
    return f"""FOR rel IN `{collection_name}`
    COLLECT disease_id = rel.end_id
    RETURN disease_id
    """


def disease_parent_query() -> str:
    # Current graphs use this collection for the DISEASES projection. Legacy
    # graphs used the same name for raw MONDO edges; SQL ignores either type.
    return """FOR rel IN `DiseaseParentEdge`
    RETURN { "start": rel.start_id, "end": rel.end_id }
    """


def protein_disease_query(last_key: str = None, limit: int = 10000) -> str:
    filter_clause = f'FILTER rel._key > "{last_key}"' if last_key else ""
    return f"""
    FOR rel IN `ProteinDiseaseEdge`
        {filter_clause}
        SORT rel._key
        LIMIT {limit}
        RETURN rel
    """


def tinx_importance_query(last_key: str = None, limit: int = 10000) -> str:
    filter_clause = f'FILTER rel._key > "{last_key}"' if last_key else ""
    return f"""
    FOR rel IN `TINXImportanceEdge`
        {filter_clause}
        SORT rel._key
        LIMIT {limit}
        RETURN rel
    """


def disease_version_query() -> str:
    return """FOR d IN `Disease` LIMIT 1 RETURN d.creation"""


class DiseaseAdapter(PharosArangoAdapter):
    associated_only: bool

    def __init__(self, credentials, database_name: str, associated_only: bool = False):
        self.associated_only = associated_only
        super().__init__(credentials=credentials, database_name=database_name)

    def get_all(self) -> Generator[List[Union[Disease, DiseaseParentEdge]], None, None]:
        db = self.get_db()
        associated_ids = None
        if self.associated_only:
            associated_ids = set()
            for collection_name in ("ProteinDiseaseEdge", "TINXImportanceEdge", "GwasTraitDiseaseEdge"):
                if db.has_collection(collection_name):
                    associated_ids.update(self.runQuery(associated_disease_ids_query(collection_name)))

        diseases = self.runQuery(disease_query())
        disease_map = {}
        rows = []
        for d in diseases:
            disease = Disease.from_dict(d)
            if associated_ids is not None and disease.id not in associated_ids:
                continue
            disease_map[disease.id] = disease
            rows.append(disease)
        yield rows

        if not db.has_collection("DiseaseParentEdge"):
            yield []
            return

        parents = self.runQuery(disease_parent_query())
        yield [
            DiseaseParentEdge(
                start_node=disease_map[rel['start']],
                end_node=disease_map[rel['end']]
            )
            for rel in parents
            if rel['start'] in disease_map and rel['end'] in disease_map
        ]

    def get_version_info_query(self) -> DataSourceDetails:
        raw_version_info = self.runQuery(disease_version_query())[0]
        return DataSourceDetails.parse_tsv(raw_version_info)


class ProteinDiseaseAdapter(PharosArangoAdapter):
    batch_size = 10_000

    @staticmethod
    def _build_details(row: dict) -> List[DiseaseAssociationDetail]:
        details = row.get('details') or []
        return [DiseaseAssociationDetail.from_dict(detail) for detail in details]

    def get_all(self) -> Generator[List[ProteinDiseaseEdge], None, None]:
        last_key = None
        while True:
            rows = list(self.runQuery(protein_disease_query(last_key=last_key, limit=self.batch_size)))
            if not rows:
                break

            yield [
                ProteinDiseaseEdge(
                    start_node=Protein(id=row['start_id']),
                    end_node=Disease(id=row['end_id'], name=row.get('end_node', {}).get('name', '')),
                    details=self._build_details(row)
                )
                for row in rows
            ]
            last_key = rows[-1]['_key']

    def get_version_info_query(self) -> DataSourceDetails:
        raw_version_info = self.runQuery(disease_version_query())[0]
        return DataSourceDetails.parse_tsv(raw_version_info)


class TINXImportanceAdapter(PharosArangoAdapter):
    batch_size = 10_000

    @staticmethod
    def _build_details(row: dict) -> List[DiseaseAssociationDetail]:
        details = row.get('details') or []
        return [DiseaseAssociationDetail.from_dict(detail) for detail in details]

    def get_all(self) -> Generator[List[TINXImportanceEdge], None, None]:
        last_key = None
        while True:
            rows = list(self.runQuery(tinx_importance_query(last_key=last_key, limit=self.batch_size)))
            if not rows:
                break

            yield [
                TINXImportanceEdge(
                    start_node=Protein(id=row['start_id']),
                    end_node=Disease(id=row['end_id'], name=row.get('end_node', {}).get('name', '')),
                    details=self._build_details(row)
                )
                for row in rows
            ]
            last_key = rows[-1]['_key']

    def get_version_info_query(self) -> DataSourceDetails:
        raw_version_info = self.runQuery(disease_version_query())[0]
        return DataSourceDetails.parse_tsv(raw_version_info)
