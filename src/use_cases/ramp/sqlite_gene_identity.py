"""UniProt-centered identities for the RaMP application projection."""
from collections import Counter
import hashlib
from types import SimpleNamespace
from src.use_cases.ramp.sqlite_gene_grouping import (
    POLICY, POLICY_NAME, POLICY_VERSION, collapse_input_matches,
)

FILE_NAME = 'uniprot-human.json.gz'


def stage_uniprot(metadata):
    datasets = [d for d in metadata['registry_datasets']
                if d['source'] == 'uniprot' and d['dataset'] == 'human']
    if len(datasets) != 1:
        raise ValueError('RaMP gene resolution requires one exact UniProt human snapshot in graph metadata')
    return datasets[0]


class GeneIdentity:
    def __init__(self, resolver, provenance):
        self.resolver = resolver
        self.provenance = provenance
        self.counts = Counter()

    @classmethod
    def from_stage(cls, metadata, registry, *, file_name=FILE_NAME):
        from src.id_resolvers.uniprot_resolver import UniProtResolver
        recorded = stage_uniprot(metadata)
        dataset = registry.resolve(recorded['snapshot_id'])
        path = dataset.file(file_name)
        with path.open('rb') as handle:
            digest = hashlib.file_digest(handle, 'sha256').hexdigest()
        provenance = {'snapshot_id': dataset.snapshot_id, 'file': file_name,
                      'sha256': digest, 'resolver': 'UniProtResolver'}
        identity = cls(None, provenance)
        identity.validate(metadata)
        identity.input_file = path
        identity.resolver = UniProtResolver(data_source=dataset, file_name=file_name)
        return identity

    def validate(self, metadata):
        recorded = stage_uniprot(metadata)
        if self.provenance['snapshot_id'] != recorded['snapshot_id']:
            raise ValueError('Gene resolver snapshot differs from the recorded graph input')
        expected = next((f.get('sha256') for f in recorded.get('files', [])
                         if f.get('path') == self.provenance['file']), None)
        if not expected or expected != self.provenance['sha256']:
            raise ValueError('Gene resolver file checksum differs from the recorded graph input')

    def groups(self, inputs):
        self.counts.clear()
        matches_by_input = {}
        for collection, identifier in sorted(set(inputs)):
            matches = self.resolver.resolve_internal([SimpleNamespace(id=identifier)])[identifier]
            accessions = {m.match for m in matches if m.match.startswith('UniProtKB:')}
            state = 'resolved' if len(accessions) == 1 else 'ambiguous' if accessions else 'unmatched'
            self.counts[f'{collection}.{state}'] += 1
            matches_by_input[(collection, identifier)] = accessions
        groups = collapse_input_matches(matches_by_input)
        self.matches_by_input = matches_by_input
        sizes = [len(values) for values in groups.canonical_accessions.values()]
        self.counts['output_groups'] = len(sizes)
        self.counts['multi_accession_groups'] = sum(size > 1 for size in sizes)
        self.counts['largest_canonical_accession_group'] = max(sizes, default=0)
        return groups

    def manifest(self):
        return {**self.provenance, 'policy': POLICY, 'policy_name': POLICY_NAME, 'policy_version': POLICY_VERSION,
                'counts': dict(sorted(self.counts.items())),
                'ambiguous_match_policy': 'All returned canonical matches are collapsed; counts describe input matches, not singleton output groups',
                'catalyzed_type_policy': 'sorted distinct source protein types joined with semicolon-space; missing type is Unknown'}
