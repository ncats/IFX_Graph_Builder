import gzip
import json
from pathlib import Path

import pytest

from src.id_resolvers.uniprot_resolver import UniProtResolver
from src.models.protein import Protein


@pytest.fixture
def uniprot_fixture_path(tmp_path: Path) -> Path:
    path = tmp_path / "uniprot-human-reviewed.json.gz"
    records = {
        "results": [
            {
                "entryType": "UniProtKB reviewed (Swiss-Prot)",
                "primaryAccession": "P04637",
                "uniProtkbId": "P53_HUMAN",
                "proteinDescription": {
                    "recommendedName": {
                        "fullName": {"value": "Cellular tumor antigen p53"}
                    }
                },
                "genes": [{"geneName": {"value": "TP53"}}],
                "uniProtKBCrossReferences": [
                    {'database': 'GeneID', 'id': '7157'},
                    {'database': 'Ensembl', 'id': 'ENST1.2',
                     'properties': [{'key': 'GeneId', 'value': 'ENSG1.3'}]},
                ],
            },
            {
                "entryType": "UniProtKB reviewed (Swiss-Prot)",
                "primaryAccession": "Q15149",
                "uniProtkbId": "PLEC_HUMAN",
                "proteinDescription": {
                    "recommendedName": {"fullName": {"value": "Plectin"}}
                },
                "genes": [{"geneName": {"value": "PLEC"}}],
            },
        ]
    }
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        json.dump(records, handle)
    return path


def test_uniprot_resolver_initializes_framework_fields(uniprot_fixture_path: Path):
    resolver = UniProtResolver(
        uniprot_json_path=str(uniprot_fixture_path),
        types=["Protein"],
        no_match_behavior="Skip",
    )

    assert resolver.types == ["Protein"]
    assert resolver.no_match_behavior.value == "Skip"
    assert "UniProtKB:P12345" not in resolver.resolve_cache


def test_uniprot_resolver_returns_empty_list_for_unmatched_id(uniprot_fixture_path: Path):
    resolver = UniProtResolver(
        uniprot_json_path=str(uniprot_fixture_path),
        types=["Protein"],
        no_match_behavior="Skip",
    )

    results = resolver.resolve_internal([Protein(id="UniProtKB:NOT_A_REAL_UNIPROT")])

    assert results["UniProtKB:NOT_A_REAL_UNIPROT"] == []


def test_uniprot_resolver_normalizes_isoform_to_base_accession(uniprot_fixture_path: Path):
    resolver = UniProtResolver(
        uniprot_json_path=str(uniprot_fixture_path),
        types=["Protein"],
        no_match_behavior="Skip",
    )

    results = resolver.resolve_internal([Protein(id="UniProtKB:Q15149-9")])

    assert len(results["UniProtKB:Q15149-9"]) == 1
    assert results["UniProtKB:Q15149-9"][0].match == "UniProtKB:Q15149"


def test_uniprot_resolver_resolves_symbol_alias(uniprot_fixture_path: Path):
    resolver = UniProtResolver(
        uniprot_json_path=str(uniprot_fixture_path),
        types=["Protein"],
        no_match_behavior="Skip",
    )

    results = resolver.resolve_internal([Protein(id="TP53")])

    assert len(results["TP53"]) == 1
    assert results["TP53"][0].match == "UniProtKB:P04637"


@pytest.mark.parametrize('identifier', ['NCBIGene:7157', 'Symbol:TP53', 'HGNC.SYMBOL:TP53', 'Ensembl:ENSG1.9'])
def test_gene_identifiers_resolve_to_canonical_protein(uniprot_fixture_path, identifier):
    resolver = UniProtResolver(uniprot_json_path=str(uniprot_fixture_path))
    results = resolver.resolve_internal([Protein(id=identifier)])
    assert [m.match for m in results[identifier]] == ['UniProtKB:P04637']


def test_primary_accession_precedes_lower_priority_alias_match(uniprot_fixture_path):
    from src.interfaces.id_resolver import IdMatch
    resolver = UniProtResolver(uniprot_json_path=str(uniprot_fixture_path))
    resolver.alias_map['UniProtKB:P04637'].insert(
        0, IdMatch(input='UniProtKB:P04637', match='UniProtKB:Q15149', context=['secondary accession']))
    result = resolver.resolve_internal([Protein(id='UniProtKB:P04637')])
    assert [m.match for m in result['UniProtKB:P04637']] == ['UniProtKB:P04637']
