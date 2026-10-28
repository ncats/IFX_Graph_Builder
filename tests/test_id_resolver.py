from typing import List, Dict
import pytest

from src.interfaces.id_resolver import IdResolver, IdMatch, NoMatchBehavior, MultiMatchBehavior
from src.models.gene import Gene, GeneticLocation
from src.models.node import EquivalentId, Node
from src.models.protein import Protein


class TempResolver(IdResolver):
    def resolve_internal(self, input_nodes: List[Node]) -> Dict[str, List[IdMatch]]:
        return {'test': []}



def test_simple_merge():
    testResolver = TempResolver(types=['Node'], no_match_behavior=NoMatchBehavior.Allow, multi_match_behavior=MultiMatchBehavior.All)
    entries = [Node(id='1'), Node(id='2')]
    oneEquivs = ['EC:1','EC:1a']
    oneAEquivs = ['EC:2']
    id_map: Dict[str, List[IdMatch]] = {
        '1': [IdMatch('1', 'Match:1', equivalent_ids=oneEquivs), IdMatch('1', 'Match:1a', equivalent_ids=oneAEquivs)]
    }

    merged_map = testResolver.get_merged_map(entries, id_map)

    assert len(merged_map) == 3
    assert len(merged_map[IdResolver.MatchKeys.matched]) == 1
    assert merged_map[IdResolver.MatchKeys.matched]['1'].id == 'Match:1'
    assert len(merged_map[IdResolver.MatchKeys.matched]['1'].xref) == 3

    assert len(merged_map[IdResolver.MatchKeys.newborns]) == 1
    assert len(merged_map[IdResolver.MatchKeys.newborns]['1']) == 1
    assert merged_map[IdResolver.MatchKeys.newborns]['1'][0].id == 'Match:1a'
    assert len(merged_map[IdResolver.MatchKeys.newborns]['1'][0].xref) == 3

    assert len(merged_map[IdResolver.MatchKeys.unmatched]) == 1
    assert merged_map[IdResolver.MatchKeys.unmatched]['2'].id == '2'


def test_relationship_endpoint_multi_match_error_does_not_select_first():
    resolver = TempResolver(types=['Node'], multi_match_behavior=MultiMatchBehavior.Error)
    entries = [Node(id='source')]
    matches = {
        'source': [
            IdMatch('source', 'Target:one'),
            IdMatch('source', 'Target:two'),
        ]
    }
    entity_map = resolver.get_merged_map(entries, matches)

    with pytest.raises(ValueError, match="resolves to multiple targets"):
        resolver.parse_entity_map(entity_map)


def test_repeated_ids_parse_xrefs_once(monkeypatch):
    resolver = TempResolver(types=['Node'])
    entries = [Node(id='1') for _ in range(100)]
    id_map = {'1': [IdMatch('1', 'Match:1', equivalent_ids=['EC:1', 'EC:2'])]}
    parsed = []
    original_parse = EquivalentId.parse

    def track_parse(value):
        parsed.append(value)
        return original_parse(value)

    monkeypatch.setattr(EquivalentId, 'parse', staticmethod(track_parse))
    merged_map = resolver.get_merged_map(entries, id_map)

    assert parsed == ['EC:1', 'EC:2']
    assert len(merged_map[IdResolver.MatchKeys.matched]['1'].xref) == 2


def test_cross_type_retype_copies_compatible_fields():
    testResolver = TempResolver(
        types=['Gene'],
        no_match_behavior=NoMatchBehavior.Allow,
        multi_match_behavior=MultiMatchBehavior.All,
        canonical_class=Protein,
    )
    entries = [Gene(id='ENSEMBL:ENSG1', calculated_properties={'gtex_tau': 0.5})]
    id_map = {
        'ENSEMBL:ENSG1': [IdMatch('ENSEMBL:ENSG1', 'IFXProtein:ABC', equivalent_ids=['ENSEMBL:ENSG1'])]
    }

    merged_map = testResolver.get_merged_map(entries, id_map, allow_retype=True)
    matched = merged_map[IdResolver.MatchKeys.matched]['ENSEMBL:ENSG1']

    assert isinstance(matched, Protein)
    assert matched.id == 'IFXProtein:ABC'
    assert matched.calculated_properties == {'gtex_tau': 0.5}


def test_cross_type_retype_rejects_incompatible_fields():
    testResolver = TempResolver(
        types=['Gene'],
        no_match_behavior=NoMatchBehavior.Allow,
        multi_match_behavior=MultiMatchBehavior.All,
        canonical_class=Protein,
    )
    entries = [Gene(id='ENSEMBL:ENSG1', location=GeneticLocation(location='1p36'))]
    id_map = {
        'ENSEMBL:ENSG1': [IdMatch('ENSEMBL:ENSG1', 'IFXProtein:ABC', equivalent_ids=['ENSEMBL:ENSG1'])]
    }

    merged_map = testResolver.get_merged_map(entries, id_map, allow_retype=True)

    assert 'ENSEMBL:ENSG1' not in merged_map[IdResolver.MatchKeys.matched]
    assert merged_map[IdResolver.MatchKeys.unmatched]['ENSEMBL:ENSG1'].id == 'ENSEMBL:ENSG1'
