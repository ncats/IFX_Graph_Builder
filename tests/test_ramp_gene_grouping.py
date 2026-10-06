from src.use_cases.ramp.sqlite_gene_grouping import collapse_input_matches


def test_transitive_components_depend_only_on_used_identifiers():
    matches = {
        ('GeneIdentifier', 'G1'): ['UniProtKB:A', 'UniProtKB:B'],
        ('GeneIdentifier', 'G2'): ['UniProtKB:B', 'UniProtKB:C'],
        ('ProteinIdentifier', 'UniProtKB:D'): ['UniProtKB:D'],
        ('GeneIdentifier', 'UNKNOWN'): [],
        ('ProteinIdentifier', 'UNKNOWN'): [],
    }
    result = collapse_input_matches(matches)
    assert result.input_to_group[('GeneIdentifier', 'G1')] == result.input_to_group[('GeneIdentifier', 'G2')]
    assert result.canonical_accessions[('UniProtKB', 'UniProtKB:A')] == ('UniProtKB:A', 'UniProtKB:B', 'UniProtKB:C')
    assert len(result.canonical_accessions) == 4
    assert result.input_to_group[('GeneIdentifier', 'UNKNOWN')] != result.input_to_group[('ProteinIdentifier', 'UNKNOWN')]
    reordered = {key: list(reversed(values)) + values for key, values in reversed(list(matches.items()))}
    assert collapse_input_matches(reordered) == result
    assert matches[('GeneIdentifier', 'G1')] == ['UniProtKB:A', 'UniProtKB:B']


def test_empty_inputs():
    result = collapse_input_matches({})
    assert result.input_to_group == result.canonical_accessions == {}
