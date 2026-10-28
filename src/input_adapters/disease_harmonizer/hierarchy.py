"""Published DISEASES MONDO hierarchy projection for graph comparison."""

from __future__ import annotations

import csv

from src.constants import DataSourceName
from src.interfaces.input_adapter import InputAdapter
from src.models.disease import Disease, DiseaseParentEdge
from src.shared.disease_consumer_contract import validated_disease_snapshot_file


HIERARCHY_MEMBER = "app_graph/disease_hierarchy_edges.tsv"
REQUIRED_COLUMNS = {
    "parent_primary_xref", "parent_ncats_disease_id", "parent_label",
    "child_primary_xref", "child_ncats_disease_id", "child_label",
    "relationship_type", "predicate", "source",
}


class DiseaseConsumerHierarchyAdapter(InputAdapter):
    """Keep the derived MONDO subset separate from primary ontology edges."""

    def __init__(self, data_source):
        self.path, _checksum = validated_disease_snapshot_file(data_source, HIERARCHY_MEMBER)
        self.version_info = data_source.version_info()

    def get_datasource_name(self):
        return DataSourceName.DiseaseHarmonizer

    def get_version(self):
        return self.version_info

    def get_all(self):
        batch = []
        seen_pairs = set()
        with self.path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            columns = reader.fieldnames or []
            missing = REQUIRED_COLUMNS - set(columns)
            if missing or len(columns) != len(set(columns)):
                raise ValueError(
                    f"{self.path} invalid hierarchy columns: missing={sorted(missing)}"
                )
            for row_number, row in enumerate(reader, start=2):
                if None in row or any(value is None for value in row.values()):
                    raise ValueError(f"{self.path}:{row_number} has malformed TSV fields")
                child_id = row["child_ncats_disease_id"].strip()
                parent_id = row["parent_ncats_disease_id"].strip()
                child_xref = row["child_primary_xref"].strip()
                parent_xref = row["parent_primary_xref"].strip()
                if not all(value.startswith("IFXDisease:") and value.split(":", 1)[1]
                           for value in (child_id, parent_id)):
                    raise ValueError(f"{self.path}:{row_number} has invalid IFXDisease endpoint")
                if not all(value.startswith("MONDO:") and value.split(":", 1)[1]
                           for value in (child_xref, parent_xref)):
                    raise ValueError(f"{self.path}:{row_number} has invalid MONDO primary xref")
                if (row["source"], row["predicate"], row["relationship_type"]) != (
                    "MONDO", "biolink:subclass_of", "parent_child"
                ):
                    raise ValueError(f"{self.path}:{row_number} has unexpected hierarchy semantics")
                pair = (child_id, parent_id)
                if child_id == parent_id or pair in seen_pairs:
                    raise ValueError(f"{self.path}:{row_number} has repeated or self hierarchy edge {pair}")
                seen_pairs.add(pair)
                batch.append(DiseaseParentEdge(
                    start_node=Disease(id=child_id),
                    end_node=Disease(id=parent_id),
                    source="MONDO",
                    hierarchy_payload=dict(row),
                ))
                if len(batch) >= self.batch_size:
                    yield batch
                    batch = []
        if not seen_pairs:
            raise ValueError(f"{self.path} has no DISEASES hierarchy rows")
        if batch:
            yield batch
