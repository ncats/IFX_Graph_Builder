from typing import Dict, List

import pytest

from src.constants import DataSourceName
from src.interfaces.id_resolver import IdResolver, IdMatch, NoMatchBehavior, MultiMatchBehavior
from src.interfaces.input_adapter import InputAdapter
from src.models.datasource_version_info import DatasourceVersionInfo
from src.models.disease import Disease, GeneDiseaseEdge, ProteinDiseaseEdge
from src.models.expression import ExpressionDetail, GeneTissueExpressionEdge, ProteinTissueExpressionEdge
from src.models.gene import Gene, MeasuredGene, MeasuredGeneEdge
from src.models.mouse_phenotype import GeneMousePhenotypeEdge, MousePhenotype, ProteinMousePhenotypeEdge
from src.models.ortholog import GeneOrthologGeneEdge, OrthologGene, ProteinOrthologGeneEdge
from src.models.pathway import GenePathwayEdge, Pathway, ProteinPathwayEdge
from src.models.protein import Protein, MeasuredProtein, MeasuredProteinEdge
from src.models.tissue import Tissue
from src.models.alliance_genome import AllianceGeneDiseaseEdge


class _RetypeGeneResolver(IdResolver):
    def resolve_internal(self, input_nodes: List[Gene]) -> Dict[str, List[IdMatch]]:
        return {
            node.id: [IdMatch(node.id, f"IFXProtein:{node.id}", equivalent_ids=[f"NCBIGene:{node.id}"])]
            for node in input_nodes
        }


class _IdentityResolver(IdResolver):
    def resolve_internal(self, input_nodes):
        return {
            node.id: [IdMatch(node.id, node.id, equivalent_ids=[node.id])]
            for node in input_nodes
        }


class _SingleBatchAdapter(InputAdapter):
    def __init__(self, entries):
        self._entries = entries

    def get_all(self):
        yield self._entries

    def get_datasource_name(self) -> DataSourceName:
        return DataSourceName.TargetGraph

    def get_version(self) -> DatasourceVersionInfo:
        return DatasourceVersionInfo(version="test", version_date=None, download_date=None)


def _resolver_map():
    return {
        "Gene": _RetypeGeneResolver(
            types=["Gene"],
            no_match_behavior=NoMatchBehavior.Skip,
            multi_match_behavior=MultiMatchBehavior.All,
            canonical_class=Protein,
        )
    }


def test_input_adapter_remaps_gene_tissue_expression_edge_to_protein_edge():
    edge = GeneTissueExpressionEdge(start_node=Gene(id="G1"), end_node=Tissue(id="T1"), details=[])
    adapter = _SingleBatchAdapter([edge])

    batches = list(adapter.get_resolved_and_provenanced_list(_resolver_map()))
    rel = batches[-1][0]

    assert isinstance(rel, ProteinTissueExpressionEdge)
    assert rel.start_node.id == "IFXProtein:G1"


def test_remapped_expression_edges_keep_independent_copies_for_multiple_matches():
    class _MultiMatchResolver(_RetypeGeneResolver):
        def resolve_internal(self, input_nodes):
            return {
                node.id: [
                    IdMatch(node.id, "IFXProtein:P1", equivalent_ids=[f"NCBIGene:{node.id}"]),
                    IdMatch(node.id, "IFXProtein:P2", equivalent_ids=[f"NCBIGene:{node.id}"]),
                ]
                for node in input_nodes
            }

    detail = ExpressionDetail(source="GTEx", number_value=2.0)
    edge = GeneTissueExpressionEdge(
        start_node=Gene(id="G1"), end_node=Tissue(id="UBERON:0001", name="source tissue"), details=[detail]
    )
    adapter = _SingleBatchAdapter([edge])
    resolver = _MultiMatchResolver(
        types=["Gene"],
        no_match_behavior=NoMatchBehavior.Skip,
        multi_match_behavior=MultiMatchBehavior.All,
        canonical_class=Protein,
    )

    tissue_resolver = _IdentityResolver(
        types=["Tissue"],
        no_match_behavior=NoMatchBehavior.Skip,
        multi_match_behavior=MultiMatchBehavior.All,
    )
    batches = list(adapter.get_resolved_and_provenanced_list({"Gene": resolver, "Tissue": tissue_resolver}))
    remapped = batches[-1]

    assert {rel.start_node.id for rel in remapped} == {"IFXProtein:P1", "IFXProtein:P2"}
    assert all(isinstance(rel, ProteinTissueExpressionEdge) for rel in remapped)
    assert all(rel.end_node.name == "source tissue" and rel.end_node.xref == [] for rel in remapped)
    remapped[0].details[0].number_value = 9.0
    remapped[0].start_node.id = "changed"
    assert remapped[1].details[0].number_value == 2.0
    assert remapped[1].start_node.id != "changed"
    assert edge.details[0].number_value == 2.0
    assert edge.start_node.id == "G1"


def test_unmatched_reporter_records_raw_nodes_and_edge_endpoints():
    class _OnlyKnownResolver(IdResolver):
        def resolve_internal(self, input_nodes):
            return {
                node.id: [IdMatch(node.id, node.id)]
                for node in input_nodes if node.id.endswith("KNOWN")
            }

    adapter = _SingleBatchAdapter([
        Gene(id="Gene:MISSING"), Gene(id="Gene:MISSING"),
        GenePathwayEdge(start_node=Gene(id="Gene:MISSING"), end_node=Pathway(id="Pathway:MISSING")),
        GenePathwayEdge(start_node=Gene(id="Gene:MISSING"), end_node=Pathway(id="Pathway:KNOWN")),
    ])
    issues = []
    adapter.set_unmatched_reporter(issues.extend)
    resolvers = {
        node_type: _OnlyKnownResolver(types=[node_type], no_match_behavior=NoMatchBehavior.Skip)
        for node_type in ("Gene", "Pathway")
    }

    list(adapter.get_resolved_and_provenanced_list(resolvers))

    assert issues == [
        {"node_type": "Gene", "source_id": "Gene:MISSING", "endpoint_role": "node", "occurrences": 2, "resolver_outcome": "skipped"},
        {"node_type": "Gene", "source_id": "Gene:MISSING", "endpoint_role": "start", "occurrences": 2, "resolver_outcome": "skipped"},
        {"node_type": "Pathway", "source_id": "Pathway:MISSING", "endpoint_role": "end", "occurrences": 1, "resolver_outcome": "skipped"},
    ]


def test_unmatched_reporter_records_allowed_and_prefix_filtered_ids():
    class _NoMatchesResolver(IdResolver):
        def resolve_internal(self, input_nodes):
            return {}

        def parse_flat_node_list_from_map(self, entity_map):
            filtered = dict(entity_map)
            filtered[self.MatchKeys.unmatched] = {
                source_id: node for source_id, node in entity_map[self.MatchKeys.unmatched].items()
                if source_id.startswith("IFXGene:")
            }
            return super().parse_flat_node_list_from_map(filtered)

    adapter = _SingleBatchAdapter([Gene(id="IFXGene:1"), Gene(id="Other:2")])
    issues = []
    adapter.set_unmatched_reporter(issues.extend)
    resolver = _NoMatchesResolver(types=["Gene"], no_match_behavior=NoMatchBehavior.Allow)

    list(adapter.get_resolved_and_provenanced_list({"Gene": resolver}))

    assert [(issue["source_id"], issue["resolver_outcome"]) for issue in issues] == [
        ("IFXGene:1", "allowed"), ("Other:2", "skipped"),
    ]


def test_allowed_unmatched_endpoint_does_not_imply_edge_was_emitted():
    class _NoMatchesResolver(IdResolver):
        def resolve_internal(self, input_nodes):
            return {}

    edge = GenePathwayEdge(
        start_node=Gene(id="Gene:ALLOWED"),
        end_node=Pathway(id="Pathway:SKIPPED"),
    )
    adapter = _SingleBatchAdapter([edge])
    issues = []
    adapter.set_unmatched_reporter(issues.extend)
    resolvers = {
        "Gene": _NoMatchesResolver(types=["Gene"], no_match_behavior=NoMatchBehavior.Allow),
        "Pathway": _NoMatchesResolver(types=["Pathway"], no_match_behavior=NoMatchBehavior.Skip),
    }

    batches = list(adapter.get_resolved_and_provenanced_list(resolvers))

    assert not any(isinstance(item, GenePathwayEdge) for batch in batches for item in batch)
    assert {(issue["source_id"], issue["resolver_outcome"]) for issue in issues} == {
        ("Gene:ALLOWED", "allowed"), ("Pathway:SKIPPED", "skipped"),
    }


def test_unmatched_reporter_records_error_when_resolver_rejects_miss():
    class _NoMatchesResolver(IdResolver):
        def resolve_internal(self, input_nodes):
            return {}

    adapter = _SingleBatchAdapter([Gene(id="Gene:MISSING")])
    issues = []
    adapter.set_unmatched_reporter(issues.extend)
    resolver = _NoMatchesResolver(types=["Gene"], no_match_behavior=NoMatchBehavior.Error)

    with pytest.raises(Exception, match="unmatched identifiers"):
        list(adapter.get_resolved_and_provenanced_list({"Gene": resolver}))

    assert issues[0]["resolver_outcome"] == "error"


def test_input_adapter_remaps_gene_disease_edge_to_protein_edge():
    edge = GeneDiseaseEdge(start_node=Gene(id="G1"), end_node=Disease(id="D1"), details=[])
    adapter = _SingleBatchAdapter([edge])

    batches = list(adapter.get_resolved_and_provenanced_list(_resolver_map()))
    rel = batches[-1][0]

    assert isinstance(rel, ProteinDiseaseEdge)
    assert rel.start_node.id == "IFXProtein:G1"


def test_input_adapter_remaps_gene_pathway_edge_to_protein_edge():
    edge = GenePathwayEdge(start_node=Gene(id="G1"), end_node=Pathway(id="PW1"), source="x")
    adapter = _SingleBatchAdapter([edge])

    batches = list(adapter.get_resolved_and_provenanced_list(_resolver_map()))
    rel = batches[-1][0]

    assert isinstance(rel, ProteinPathwayEdge)
    assert rel.start_node.id == "IFXProtein:G1"
    assert rel.source == "x"


def test_input_adapter_remaps_gene_ortholog_edge_to_protein_ortholog_edge():
    edge = GeneOrthologGeneEdge(start_node=Gene(id="G1"), end_node=OrthologGene(id="OG1"))
    adapter = _SingleBatchAdapter([edge])

    batches = list(adapter.get_resolved_and_provenanced_list(_resolver_map()))
    rel = batches[-1][0]

    assert isinstance(rel, ProteinOrthologGeneEdge)
    assert rel.start_node.id == "IFXProtein:G1"
    assert rel.end_node.id == "OG1"


def test_input_adapter_remaps_gene_mouse_phenotype_edge_to_protein_mouse_phenotype_edge():
    edge = GeneMousePhenotypeEdge(start_node=Gene(id="G1"), end_node=MousePhenotype(id="MP:1"), details=[])
    adapter = _SingleBatchAdapter([edge])

    batches = list(adapter.get_resolved_and_provenanced_list(_resolver_map()))
    rel = batches[-1][0]

    assert isinstance(rel, ProteinMousePhenotypeEdge)
    assert rel.start_node.id == "IFXProtein:G1"
    assert rel.end_node.id == "MP:1"


def test_input_adapter_leaves_unmapped_gene_edge_classes_unchanged():
    edge = AllianceGeneDiseaseEdge(start_node=Gene(id="G1"), end_node=Disease(id="D1"), evidence_terms=["ECO:1"])
    adapter = _SingleBatchAdapter([edge])

    batches = list(adapter.get_resolved_and_provenanced_list(_resolver_map()))
    rel = batches[-1][0]

    assert isinstance(rel, AllianceGeneDiseaseEdge)
    assert isinstance(rel.start_node, Protein)
    assert rel.start_node.id == "IFXProtein:G1"
    assert rel.evidence_terms == ["ECO:1"]


def test_input_adapter_keeps_measured_gene_edge_start_node_typed_as_measured_gene():
    edge = MeasuredGeneEdge(
        start_node=MeasuredGene(id="Ensembl:G1", symbol="GENE1"),
        end_node=Gene(id="Ensembl:G1"),
    )
    adapter = _SingleBatchAdapter([edge])

    batches = list(adapter.get_resolved_and_provenanced_list({
        "Gene": _IdentityResolver(
            types=["Gene"],
            no_match_behavior=NoMatchBehavior.Skip,
            multi_match_behavior=MultiMatchBehavior.All,
        )
    }))
    rel = batches[-1][0]

    assert isinstance(rel, MeasuredGeneEdge)
    assert isinstance(rel.start_node, MeasuredGene)
    assert rel.start_node.id == "Ensembl:G1"
    assert isinstance(rel.end_node, Gene)
    assert rel.end_node.id == "Ensembl:G1"


def test_input_adapter_keeps_measured_protein_edge_start_node_typed_as_measured_protein():
    edge = MeasuredProteinEdge(
        start_node=MeasuredProtein(id="UniProtKB:P12345", name="Protein 1"),
        end_node=Protein(id="UniProtKB:P12345"),
    )
    adapter = _SingleBatchAdapter([edge])

    batches = list(adapter.get_resolved_and_provenanced_list({
        "Protein": _IdentityResolver(
            types=["Protein"],
            no_match_behavior=NoMatchBehavior.Skip,
            multi_match_behavior=MultiMatchBehavior.All,
        )
    }))
    rel = batches[-1][0]

    assert isinstance(rel, MeasuredProteinEdge)
    assert isinstance(rel.start_node, MeasuredProtein)
    assert rel.start_node.id == "UniProtKB:P12345"
    assert isinstance(rel.end_node, Protein)
    assert rel.end_node.id == "UniProtKB:P12345"


def test_input_adapter_preserves_existing_edge_sources_and_provenance():
    edge = ProteinPathwayEdge(
        start_node=Protein(id="UniProtKB:P1"),
        end_node=Pathway(id="Reactome:R-HSA-1"),
        source="Reactome",
    )
    edge.sources = ["BioPlex\t3.0 (293T)\t2024-01-19\t2026-04-24", "Reactome\t96\t2026-03-24\t2026-04-24"]
    edge.provenance = "BioPlex\t3.0 (293T)\t2024-01-19\t2026-04-24"
    adapter = _SingleBatchAdapter([edge])

    batches = list(adapter.get_resolved_and_provenanced_list({
        "Protein": _IdentityResolver(
            types=["Protein"],
            no_match_behavior=NoMatchBehavior.Skip,
            multi_match_behavior=MultiMatchBehavior.All,
        ),
        "Pathway": _IdentityResolver(
            types=["Pathway"],
            no_match_behavior=NoMatchBehavior.Skip,
            multi_match_behavior=MultiMatchBehavior.All,
        ),
    }))
    rel = batches[-1][0]

    assert rel.sources == ["BioPlex\t3.0 (293T)\t2024-01-19\t2026-04-24", "Reactome\t96\t2026-03-24\t2026-04-24"]
    assert rel.provenance == "BioPlex\t3.0 (293T)\t2024-01-19\t2026-04-24"
