"""TARGETS consumer payload and selected relationship adapters."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

from src.constants import DataSourceName
from src.input_adapters.target_graph.protein_nodes_and_edges import TGProteinFileBase
from src.interfaces.input_adapter import InputAdapter
from src.models.gene import Gene
from src.models.protein import Protein
from src.models.transcript import (
    GeneProteinEdge,
    GeneTranscriptEdge,
    IsoformProteinEdge,
    ParentLinkDetail,
    Transcript,
    TranscriptProteinEdge,
)
from src.shared.targetgraph_parser import (
    TargetGraphGeneParser,
    TargetGraphTranscriptParser,
)


_PAYLOAD_FILES = {
    "gene": "app_graph/target_gene_payload.tsv",
    "transcript": "app_graph/target_transcript_payload.tsv",
    "protein": "app_graph/target_protein_payload.tsv",
}
_ID_PREFIX = {"gene": "IFXGene:", "transcript": "IFXTranscript:", "protein": "IFXProtein:"}
_RELATIONSHIPS = {
    "gene_to_transcript": ("IFXGene:", "IFXTranscript:", GeneTranscriptEdge, Gene, Transcript),
    "transcript_to_protein": ("IFXTranscript:", "IFXProtein:", TranscriptProteinEdge, Transcript, Protein),
    "gene_to_protein_resolved_parent": ("IFXGene:", "IFXProtein:", GeneProteinEdge, Gene, Protein),
    "isoform_of_canonical_protein": ("IFXProtein:", "IFXProtein:", IsoformProteinEdge, Protein, Protein),
}
_DEFERRED_SELECTED = {"gene_to_canonical_protein", "uniprot_entry_to_displayed_protein"}
_RELATIONSHIP_COLUMNS = {
    "subject_id", "predicate", "biolink_predicate", "object_id", "relation_kind",
    "status", "evidence_source", "evidence", "confidence", "source_version",
}


def _checked_reader(path: Path, required: set[str]):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} lacks required columns: {sorted(missing)}")
        for row_number, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(f"{path}:{row_number} has too many fields")
            yield row_number, row


class _TargetConsumerBase(InputAdapter):
    def __init__(self, data_source):
        self.version_info = data_source.version_info()
        self.expected_release = f"v{data_source.version}"
        self.data_source = data_source

    def get_datasource_name(self):
        return DataSourceName.TargetHarmonizer

    def get_version(self):
        return self.version_info


class _PayloadAdapter(_TargetConsumerBase):
    entity_type: str
    required_columns: frozenset[str] = frozenset()

    def __init__(self, data_source):
        super().__init__(data_source)
        file_name = _PAYLOAD_FILES[self.entity_type]
        self.path = Path(data_source.file(file_name))
        if not self.path.is_file():
            raise FileNotFoundError(self.path)

    def _build_node(self, row):
        raise NotImplementedError

    def _include_row(self, row):
        return True

    def get_all(self):
        batch = []
        seen_ids = set()
        for row_number, row in _checked_reader(
            self.path, {"ifx_id", "release_version", "aliases"} | self.required_columns
        ):
            identifier = row["ifx_id"]
            if not identifier.startswith(_ID_PREFIX[self.entity_type]):
                raise ValueError(f"{self.path}:{row_number} has invalid {self.entity_type} ID {identifier!r}")
            if row["release_version"] != self.expected_release:
                raise ValueError(f"{self.path}:{row_number} has release {row['release_version']!r}, expected {self.expected_release!r}")
            if identifier in seen_ids:
                raise ValueError(f"{self.path}:{row_number} duplicates {identifier}")
            seen_ids.add(identifier)
            if not self._include_row(row):
                continue
            node = self._build_node(row)
            node.target_payload = dict(row)
            batch.append(node)
            if len(batch) >= self.batch_size:
                yield batch
                batch = []
        if batch:
            yield batch


class TargetConsumerGeneAdapter(_PayloadAdapter):
    entity_type = "gene"

    def _build_node(self, row):
        gene = Gene(id=row["ifx_id"])
        gene.full_name = TargetGraphGeneParser.get_gene_name(row)
        gene.location = TargetGraphGeneParser.get_gene_location(row)
        gene.created = TargetGraphGeneParser.get_creation_date(row)
        gene.updated = TargetGraphGeneParser.get_updated_time(row)
        gene.type = TargetGraphGeneParser.get_gene_type(row)
        gene.mapping_ratio = TargetGraphGeneParser.get_mapping_ratio(row)
        gene.symbol = TargetGraphGeneParser.get_symbol(row)
        for field in (
            "Description_Provenance", "Location_Provenance", "Ensembl_ID_Provenance",
            "NCBI_ID_Provenance", "HGNC_ID_Provenance", "Symbol_Provenance",
        ):
            target_field = "Name_Provenance" if field == "Description_Provenance" else field
            setattr(gene, target_field, row.get(field) or None)
        return gene


class TargetConsumerTranscriptAdapter(_PayloadAdapter):
    entity_type = "transcript"

    def _build_node(self, row):
        transcript = Transcript(id=row["ifx_id"])
        transcript.created = TargetGraphTranscriptParser.get_creation_date(row)
        transcript.updated = TargetGraphTranscriptParser.get_updated_time(row)
        transcript.location = TargetGraphTranscriptParser.get_transcript_location(row)
        transcript.type = TargetGraphTranscriptParser.get_transcript_type(row)
        transcript.support_level = TargetGraphTranscriptParser.get_transcript_support_level(row)
        transcript.is_canonical = TargetGraphTranscriptParser.get_transcript_is_canonical(row)
        transcript.MANE_select = TargetGraphTranscriptParser.get_mane_select(row)
        transcript.status = TargetGraphTranscriptParser.get_transcript_status(row)
        transcript.ensembl_version = TargetGraphTranscriptParser.get_transcript_version(row)
        transcript.Ensembl_Transcript_ID_Provenance = row.get("Ensembl_Transcript_ID_Provenance") or None
        transcript.RefSeq_Provenance = row.get("RefSeq_Provenance") or None
        return transcript


class TargetConsumerProteinAdapter(_PayloadAdapter):
    entity_type = "protein"
    required_columns = frozenset({"is_canonical"})

    def __init__(self, data_source, canonical_only=False):
        self.canonical_only = canonical_only
        super().__init__(data_source)

    def _include_row(self, row):
        if row["is_canonical"] not in {"True", "False"}:
            raise ValueError(f"{self.path} has invalid is_canonical value {row['is_canonical']!r} for {row['ifx_id']}")
        return not self.canonical_only or row["is_canonical"] == "True"

    def _build_node(self, row):
        return TGProteinFileBase.build_protein_obj({**row, "ncats_protein_id": row["ifx_id"]})


class TargetConsumerRelationshipAdapter(_TargetConsumerBase):
    """Emit only selected structural links for the three-node comparison."""

    def __init__(self, data_source):
        super().__init__(data_source)
        self.path = Path(data_source.file("app_graph/target_relationships.tsv"))
        if not self.path.is_file():
            raise FileNotFoundError(self.path)

    def get_all(self):
        batches = {kind: [] for kind in _RELATIONSHIPS}
        selected_counts = Counter()
        for row_number, row in _checked_reader(self.path, _RELATIONSHIP_COLUMNS):
            if row["status"] not in {"selected", "candidate", "source_claim"}:
                raise ValueError(f"{self.path}:{row_number} has unknown status {row['status']!r}")
            if row["status"] != "selected":
                continue
            kind = row["relation_kind"]
            selected_counts[kind] += 1
            if kind in _DEFERRED_SELECTED:
                continue
            if kind not in _RELATIONSHIPS:
                raise ValueError(f"{self.path}:{row_number} has unmodeled selected kind {kind!r}")
            start_prefix, end_prefix, edge_class, start_class, end_class = _RELATIONSHIPS[kind]
            if not row["subject_id"].startswith(start_prefix) or not row["object_id"].startswith(end_prefix):
                raise ValueError(f"{self.path}:{row_number} has wrong endpoint types for {kind}")
            edge = edge_class(start_node=start_class(id=row["subject_id"]), end_node=end_class(id=row["object_id"]))
            edge.details = [ParentLinkDetail(
                resolution_source=row["evidence_source"],
                **row,
            )]
            batches[kind].append(edge)
            if len(batches[kind]) >= self.batch_size:
                yield batches[kind]
                batches[kind] = []
        for batch in batches.values():
            if batch:
                yield batch
        print(f"TARGETS selected relationship counts: {dict(selected_counts)}")
        for kind in sorted(_DEFERRED_SELECTED):
            print(f"TARGETS selected {kind} rows outside three-node comparison: {selected_counts[kind]}")
