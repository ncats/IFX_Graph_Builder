"""UniProt display annotations and lookup aliases, separate from identity resolution."""
from collections import Counter
import gzip
import hashlib
import ijson


class ProteinAnnotations:
    def __init__(self, names=None, provenance=None, *, aliases=None):
        self.names = names or {}
        self.provenance = provenance or {}
        self.aliases = {key: tuple(sorted(set(values))) for key, values in (aliases or {}).items()}

    @staticmethod
    def lookup_identifiers(record):
        """Explicit UniProt assertions only; never feed this expansion into grouping."""
        from src.shared.uniprot_parser import UniProtParser
        identifiers = set()

        def add(prefix, value):
            if isinstance(value, str) and value.strip() and value.strip() != '-':
                identifiers.add(prefix + ':' + value.strip())

        add('UniProtKB', record.get('primaryAccession'))
        for value in record.get('secondaryAccessions') or []:
            add('UniProtKB', value)
        for gene in record.get('genes') or []:
            add('HGNC.SYMBOL', (gene.get('geneName') or {}).get('value'))
            for synonym in gene.get('synonyms') or []:
                add('HGNC.SYMBOL', synonym.get('value'))
        for xref in UniProtParser.find_cross_refs(record, 'GeneID'):
            add('NCBIGene', xref.get('id'))
        for xref in UniProtParser.find_cross_refs(record, 'HGNC'):
            value = xref.get('id')
            if isinstance(value, str):
                value = value.strip()
                if value.upper().startswith('HGNC:'):
                    value = value[5:]
                add('hgnc', value)
        for xref in UniProtParser.find_cross_refs(record, 'Ensembl'):
            values = [xref.get('id')]
            values.extend(p.get('value') for p in xref.get('properties') or []
                          if p.get('key') in ('GeneId', 'ProteinId'))
            for value in values:
                if isinstance(value, str):
                    add('Ensembl', UniProtParser.trim_version(value.strip()))
        return tuple(sorted(identifiers))

    @classmethod
    def from_file(cls, path, provenance):
        from src.shared.uniprot_parser import UniProtParser
        with path.open('rb') as handle:
            if hashlib.file_digest(handle, 'sha256').hexdigest() != provenance['sha256']:
                raise ValueError('Protein annotation file differs from verified UniProt resolver input')
        names, aliases = {}, {}
        with gzip.open(path, 'rb') as handle:
            for record in ijson.items(handle, 'results.item'):
                name = UniProtParser.get_full_name(record)
                if name and record.get('primaryAccession'):
                    names['UniProtKB:' + record['primaryAccession']] = name
                if record.get('primaryAccession'):
                    aliases['UniProtKB:' + record['primaryAccession']] = cls.lookup_identifiers(record)
        return cls(names, {k: provenance[k] for k in ('snapshot_id', 'file', 'sha256')}, aliases=aliases)

    def name(self, accession):
        return self.names.get(accession)

    def manifest(self):
        counts = Counter(identifier.split(':', 1)[0] for values in self.aliases.values() for identifier in values)
        return {**self.provenance, 'named_accessions': len(self.names),
                'name_policy': 'UniProt recommended full name, then submitted, then alternative; primary accession only',
                'lookup_alias_policy_version': 1,
                'lookup_alias_policy': 'Primary/secondary accessions, gene names/synonyms, GeneID, HGNC and Ensembl transcript/GeneId/ProteinId; Ensembl versions stripped; only already matched records; dataSource=uniprot; no regrouping; shared aliases retain every group',
                'file_accession_alias_counts_by_namespace': dict(sorted(counts.items()))}
