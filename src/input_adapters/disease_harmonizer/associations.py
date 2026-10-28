"""DISEASES pair-level support with its linked source assertions."""

from __future__ import annotations

import csv

from src.constants import DataSourceName
from src.interfaces.input_adapter import InputAdapter
from src.models.disease import Disease, DiseaseAssociationDetail, GeneDiseaseEdge
from src.models.gene import Gene
from src.shared.disease_consumer_contract import validated_disease_snapshot_file


SUMMARY_MEMBER = "app_graph/disease_gene_associations.tsv"
ASSERTION_MEMBER = "app_graph/disease_gene_source_assertions.tsv"
SUMMARY_COLUMNS = {
    "ncats_disease_id", "ncats_gene_id", "support_tier", "support_score_raw",
    "association_ids", "n_supporting_sources", "n_independent_origins",
}
ASSERTION_COLUMNS = {
    "association_id", "ncats_disease_id", "association_source", "negated",
    "subject_type", "ncats_gene_id", "ncats_protein_id", "subject_id_raw",
    "disease_id_raw", "disease_mapping_basis", "predicate", "evidence_type",
    "evidence_label", "evidence_score", "evidence_origins", "publications",
}


def _rows(path, required):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = reader.fieldnames or []
        missing = required - set(columns)
        if missing or len(columns) != len(set(columns)):
            raise ValueError(f"{path} invalid columns: missing={sorted(missing)}")
        for row_number, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"{path}:{row_number} has malformed TSV fields")
            yield row_number, row


class DiseaseConsumerAssociationAdapter(InputAdapter):
    """Emit one gene-disease edge per published support pair."""

    def __init__(self, data_source):
        self.summary_path, _ = validated_disease_snapshot_file(data_source, SUMMARY_MEMBER)
        self.assertion_path, _ = validated_disease_snapshot_file(data_source, ASSERTION_MEMBER)
        self.version_info = data_source.version_info()

    def get_datasource_name(self):
        return DataSourceName.DiseaseHarmonizer

    def get_version(self):
        return self.version_info

    @staticmethod
    def _detail(row):
        raw = {key: value for key, value in row.items() if key != "association_id"}
        source = row["association_source"]
        evidence_type = row["evidence_type"]
        label = row["evidence_label"]
        score = row["evidence_score"]
        zscore = None
        if evidence_type == "jensen_textmining":
            try:
                zscore = float(score)
            except ValueError as exc:
                raise ValueError(f"Invalid Jensen text-mining z-score {score!r}") from exc
        return DiseaseAssociationDetail(
            source=source,
            source_id=row["disease_id_raw"] or None,
            original_sources=[value for value in row["evidence_origins"].split("|") if value],
            evidence_terms=[label] if source == "Orphanet" and label else [],
            zscore=zscore,
            url=(row["publications"] or None) if source == "ClinGen" else None,
            predicate=row["predicate"] or None,
            negated=row["negated"] == "True",
            source_subject_id=row["subject_id_raw"] or None,
            subject_type=row["subject_type"] or None,
            source_protein_id=row["ncats_protein_id"] or None,
            disease_mapping_basis=row["disease_mapping_basis"] or None,
            evidence_type=evidence_type or None,
            evidence_label=label or None,
            evidence_score=score or None,
            publications_raw=row["publications"] or None,
            assertion_payload=raw,
        )

    def get_all(self):
        assertions = {}
        for row_number, row in _rows(self.assertion_path, ASSERTION_COLUMNS):
            association_id = row["association_id"]
            if not association_id or association_id in assertions:
                raise ValueError(f"{self.assertion_path}:{row_number} has duplicate/empty association_id")
            if not row["ncats_disease_id"].startswith("IFXDisease:"):
                raise ValueError(f"{self.assertion_path}:{row_number} has invalid disease endpoint")
            if row["subject_type"] not in {"gene", "protein"}:
                raise ValueError(f"{self.assertion_path}:{row_number} has invalid subject_type")
            expected_id = row["ncats_gene_id"] if row["subject_type"] == "gene" else row["ncats_protein_id"]
            expected_prefix = "IFXGene:" if row["subject_type"] == "gene" else "IFXProtein:"
            if not expected_id.startswith(expected_prefix):
                raise ValueError(f"{self.assertion_path}:{row_number} has invalid subject endpoint")
            if row["negated"] not in {"True", "False"}:
                raise ValueError(f"{self.assertion_path}:{row_number} has invalid negated flag")
            assertions[association_id] = row
        if not assertions:
            raise ValueError(f"{self.assertion_path} has no assertions")

        batch = []
        seen_pairs = set()
        for row_number, row in _rows(self.summary_path, SUMMARY_COLUMNS):
            disease_id = row["ncats_disease_id"]
            gene_id = row["ncats_gene_id"]
            pair = (gene_id, disease_id)
            if not gene_id.startswith("IFXGene:") or not disease_id.startswith("IFXDisease:"):
                raise ValueError(f"{self.summary_path}:{row_number} has invalid pair endpoint")
            if pair in seen_pairs:
                raise ValueError(f"{self.summary_path}:{row_number} repeats pair {pair}")
            seen_pairs.add(pair)
            linked_ids = row["association_ids"].split("|")
            if not linked_ids or any(not value for value in linked_ids):
                raise ValueError(f"{self.summary_path}:{row_number} has no assertion links")
            details = []
            for association_id in linked_ids:
                assertion = assertions.pop(association_id, None)
                if assertion is None:
                    raise ValueError(
                        f"{self.summary_path}:{row_number} has missing/repeated assertion {association_id}"
                    )
                if assertion["ncats_disease_id"] != disease_id or (
                    assertion["ncats_gene_id"] and assertion["ncats_gene_id"] != gene_id
                ):
                    raise ValueError(
                        f"{self.summary_path}:{row_number} assertion {association_id} has different endpoints"
                    )
                details.append(self._detail(assertion))
            batch.append(GeneDiseaseEdge(
                start_node=Gene(id=gene_id),
                end_node=Disease(id=disease_id),
                details=details,
                disease_harmonizer_support={
                    key: value for key, value in row.items() if key != "association_ids"
                },
            ))
            if len(batch) >= self.batch_size:
                yield batch
                batch = []
        if not seen_pairs:
            raise ValueError(f"{self.summary_path} has no associations")
        if assertions:
            raise ValueError(f"{self.assertion_path} has {len(assertions)} unlinked assertions")
        if batch:
            yield batch
