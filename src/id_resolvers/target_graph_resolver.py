from typing import List, Any, Generator, Optional

from src.constants import Prefix
from src.id_resolvers.sqlite_cache_resolver import SqliteCacheResolver, MatchingPair
from src.interfaces.id_resolver import IdMatch
from src.models.node import EquivalentId
from src.shared.targetgraph_parser import TargetGraphGeneParser, TargetGraphTranscriptParser, TargetGraphProteinParser, \
    TargetGraphParser

def _resolve_canonical_class(type_name: Optional[str]):
    if type_name is None:
        return None
    from src.models.protein import Protein
    from src.models.gene import Gene
    registry = {
        "Protein": Protein,
        "Gene": Gene,
    }
    cls = registry.get(type_name)
    if cls is None:
        raise ValueError(f"Unknown canonical_type '{type_name}'. Known types: {list(registry.keys())}")
    return cls

scores = {
    "exact": 0,
    "isoform": 0.5,
    "primary accession": 1,
    "secondary accession": 2,
    "NCBI Gene ID": 2.5,
    "uniprot kb": 3,
    "symbol": 4,
    "full name": 5,
    "Ensembl": 6,
    "STRING": 7,
    "RefSeq": 8,
    "short name": 9,
    "synonym": 10
}


class TargetGraphResolver(SqliteCacheResolver):
    name = "TargetGraph Resolver"

    parsers: List[TargetGraphParser]

    @staticmethod
    def sort_matches(match_list: List[IdMatch]):
        match_list.sort(key=lambda x: scores.get(x.context[0], float('inf')))
        return match_list


    def get_version_info(self) -> str:
        version_info = []
        for parser in self.parsers:
            version_info.append(parser.get_version_info())
        return '\t'.join(version_info)


    def matching_ids(self) -> Generator[MatchingPair, Any, None]:
        for parser in self.parsers:
            yield from self.get_one_match(parser)

    def get_one_match(self, parser):
        for line in parser.all_rows():
            id = parser.get_id(line)
            equiv_ids = parser.get_equivalent_ids(line)
            yield MatchingPair(id=id, match=id, type='exact')
            for equiv_ids in equiv_ids:
                yield MatchingPair(id=id, match=equiv_ids.id_str(), type=equiv_ids.type.value)


class TargetGraphProteinResolver(TargetGraphResolver):
    name = "TargetGraph Protein Resolver"
    parsers: List[TargetGraphProteinParser]

    def __init__(self, data_source, **kwargs):
        file_paths = [str(data_source.file("protein_ids.tsv"))]
        additional_ids = str(data_source.file("uniprot_mapping.csv"))
        self.parsers = [
            TargetGraphProteinParser(file_path=path, additional_id_file_path=additional_ids)
            for path in file_paths]
        TargetGraphResolver.__init__(self, **kwargs)


class TargetGraphGeneResolver(TargetGraphResolver):
    name = "TargetGraph Protein Resolver"
    parsers: List[TargetGraphGeneParser]

    def __init__(self, data_source, **kwargs):
        file_path = str(data_source.file("gene_ids.tsv"))
        self.parsers = [TargetGraphGeneParser(file_path=file_path)]
        TargetGraphResolver.__init__(self, **kwargs)


class TargetGraphTranscriptResolver(TargetGraphResolver):
    name = "TargetGraph Transcript Resolver"
    parsers: List[TargetGraphTranscriptParser]

    def __init__(self, data_source, **kwargs):
        file_path = str(data_source.file("transcript_ids.tsv"))
        self.parsers = [TargetGraphTranscriptParser(file_path=file_path)]
        TargetGraphResolver.__init__(self, **kwargs)

class TCRDTargetResolver(TargetGraphResolver):
    name = "TCRD Target Resolver"
    protein_parsers: List[TargetGraphProteinParser]
    gene_parser: TargetGraphGeneParser
    transcript_parser: TargetGraphTranscriptParser
    reviewed_only: bool
    collapse_to_canonical: bool

    def get_version_info(self) -> str:
        version_info = [
            f"reviewed only: {self.reviewed_only}",
            f"collapse to canonical: {self.collapse_to_canonical}",
        ]
        version_info.append(self.gene_parser.get_version_info())
        version_info.append(self.transcript_parser.get_version_info())
        for parser in self.protein_parsers:
            version_info.append(parser.get_version_info())
        return '\t'.join(version_info)

    def __init__(
        self,
        gene_data_source,
        transcript_data_source,
        protein_data_source,
        reviewed_only: bool = False,
        collapse_to_canonical: bool = False,
        canonical_type: str | None = None,
        **kwargs,
    ):
        gene_file_path = str(gene_data_source.file("gene_ids.tsv"))
        transcript_file_path = str(transcript_data_source.file("transcript_ids.tsv"))
        protein_file_paths = [str(protein_data_source.file("protein_ids.tsv"))]
        additional_ids = str(protein_data_source.file("uniprot_mapping.csv"))

        self.parsers = []
        self.protein_parsers = [
            TargetGraphProteinParser(file_path=path, additional_id_file_path=additional_ids)
            for path in protein_file_paths]
        self.gene_parser = TargetGraphGeneParser(file_path=gene_file_path)
        self.transcript_parser = TargetGraphTranscriptParser(file_path=transcript_file_path)
        self.collapse_to_canonical = collapse_to_canonical
        self.reviewed_only = reviewed_only

        TargetGraphResolver.__init__(self, canonical_class=_resolve_canonical_class(canonical_type), **kwargs)

    def matching_ids(self) -> Generator[MatchingPair, Any, None]:
        transcript_ids, transcript_id_idx, _ = self.get_transcript_ids()

        gene_ids, gene_ids_idx = self.get_gene_ids()
        protein_ids, protein_transcript_map, protein_gene_map = self.get_protein_ids(
            self.reviewed_only, self.collapse_to_canonical
        )

        missing_transcripts = set()
        missing_genes = set()

        for protein_ifx_id in protein_ids.keys():
            if protein_ifx_id in protein_transcript_map and len(protein_transcript_map[protein_ifx_id]) > 0:
                for match in protein_transcript_map[protein_ifx_id]:
                    if match not in transcript_id_idx:
                        missing_transcripts.add(match)
                    else:
                        transcript_ifx_ids = transcript_id_idx[match]
                        for transcript_ifx_id in transcript_ifx_ids:
                            equivalent_ids = transcript_ids[transcript_ifx_id]
                            for transcript_alias in equivalent_ids:
                                protein_ids[protein_ifx_id].add(MatchingPair(id=protein_ifx_id, match=transcript_alias.match, type=transcript_alias.type))

            if protein_ifx_id in protein_gene_map and len(protein_gene_map[protein_ifx_id]) > 0:
                for match in protein_gene_map[protein_ifx_id]:
                    if match not in gene_ids_idx:
                        missing_genes.add(match)
                    else:
                        gene_ifx_ids = gene_ids_idx[match]
                        for gene_ifx_id in gene_ifx_ids:
                            equivalent_ids = gene_ids[gene_ifx_id]
                            for gene_alias in equivalent_ids:
                                protein_ids[protein_ifx_id].add(MatchingPair(id=protein_ifx_id, match=gene_alias.match, type=gene_alias.type))

        print('missing transcripts')
        for tx in missing_transcripts:
            print(f"\t{tx}")

        print('missing genes')
        for gene in missing_genes:
            print(f"\t{gene}")

        for p in protein_ids.keys():
            for match in protein_ids[p]:
                yield match

    def get_gene_ids(self):
        gene_map = {}
        gene_id_idx = {}

        for line in self.gene_parser.all_rows():
            gene_id = self.gene_parser.get_id(line)

            equiv_ids = self.gene_parser.get_equivalent_ids(line)
            ids = [MatchingPair(id=gene_id, match=gene_id, type='exact')]
            for equiv_ids in equiv_ids:
                equiv_id_str = equiv_ids.id_str()
                ids.append(MatchingPair(id=gene_id, match=equiv_id_str, type=equiv_ids.type.value))
                if equiv_id_str not in gene_id_idx:
                    gene_id_idx[equiv_id_str] = set()
                gene_id_idx[equiv_id_str].add(gene_id)
            gene_map[gene_id] = set(ids)
            gene_id_idx.setdefault(gene_id, set()).add(gene_id)

        return gene_map, gene_id_idx

    def get_protein_ids(self, reviewed_only: bool, collapse_to_canonical: bool):
        protein_ids = {}
        protein_transcript_map = {}
        protein_gene_map = {}
        protein_rows = []
        canonical_protein_ids = set()

        for parser in self.protein_parsers:
            for line in parser.all_rows():
                if not parser.get_uniprot_id(line):
                    continue
                if reviewed_only and not parser.get_uniprot_reviewed(line):
                    continue
                protein_id = parser.get_id(line)
                is_canonical = parser.get_is_canonical(line) is True
                if is_canonical or not collapse_to_canonical:
                    canonical_protein_ids.add(protein_id)
                protein_rows.append((parser, line, protein_id, is_canonical))

        for parser, line, protein_id, is_canonical in protein_rows:
            target_protein_id = protein_id
            if collapse_to_canonical and not is_canonical:
                canonical_protein_id = parser.get_isoform_id(line)
                if canonical_protein_id not in canonical_protein_ids:
                    continue
                target_protein_id = canonical_protein_id

            if target_protein_id not in protein_transcript_map:
                protein_transcript_map[target_protein_id] = set()
            if target_protein_id not in protein_gene_map:
                protein_gene_map[target_protein_id] = set()

            if target_protein_id not in protein_ids:
                protein_ids[target_protein_id] = set()

            match_type = 'exact' if target_protein_id == protein_id else 'isoform'
            protein_ids[target_protein_id].add(
                MatchingPair(id=target_protein_id, match=protein_id, type=match_type)
            )

            for equiv_id in parser.get_equivalent_ids(line):
                protein_ids[target_protein_id].add(
                    MatchingPair(id=target_protein_id, match=equiv_id.id_str(), type=equiv_id.type.value)
                )

            transcript_ids = parser.get_transcript_ids(line)
            gene_id = parser.get_parent_gene_id(line)

            for transcript_id in transcript_ids:
                transcript_id_to_use = EquivalentId(id=transcript_id, type=Prefix.ENSEMBL).id_str()
                protein_transcript_map[target_protein_id].add(transcript_id_to_use)

            # Collapsing an isoform does not make its gene-parent assertion a
            # parent assertion about the canonical protein.
            if collapse_to_canonical and not is_canonical:
                continue
            if gene_id is not None:
                gene_id_to_use = gene_id
                protein_gene_map[target_protein_id].add(gene_id_to_use)
            elif 'parent_ifx_gene_id' not in parser.fieldnames:
                for legacy_gene_id in (parser.get_gene_id(line) or '').split('|'):
                    if legacy_gene_id:
                        protein_gene_map[target_protein_id].add(
                            EquivalentId(id=legacy_gene_id, type=Prefix.NCBIGene).id_str()
                        )

        return protein_ids, protein_transcript_map, protein_gene_map

    def get_transcript_ids(self):
        transcript_gene_map = {}
        transcript_id_idx = {}
        transcript_ids = {}
        for line in self.transcript_parser.all_rows():
            transcript_id = self.transcript_parser.get_id(line)
            if transcript_id not in transcript_gene_map:
                transcript_gene_map[transcript_id] = set()

            if 'parent_ifx_gene_id' in self.transcript_parser.fieldnames:
                parent_id = self.transcript_parser.get_parent_gene_id(line)
                if not parent_id:
                    raise ValueError(f"Transcript {transcript_id} has no parent_ifx_gene_id")
                transcript_gene_map[transcript_id].add(parent_id)
            else:
                ensg_id = self.transcript_parser.get_associated_ensg_id(line)
                ncbi_id = self.transcript_parser.get_associated_ncbi_id(line)
                if ensg_id:
                    transcript_gene_map[transcript_id].add(
                        EquivalentId(id=ensg_id, type=Prefix.ENSEMBL).id_str()
                    )
                if ncbi_id:
                    transcript_gene_map[transcript_id].add(
                        EquivalentId(id=ncbi_id, type=Prefix.NCBIGene).id_str()
                    )

            equiv_ids = self.transcript_parser.get_equivalent_ids(line)
            ids = [MatchingPair(id=transcript_id, match=transcript_id, type='exact')]
            for equiv_ids in equiv_ids:
                transcript_id_str = equiv_ids.id_str()
                ids.append(MatchingPair(id=transcript_id, match=transcript_id_str, type=equiv_ids.type.value))
                if transcript_id_str not in transcript_id_idx:
                    transcript_id_idx[transcript_id_str] = set()
                transcript_id_idx[transcript_id_str].add(transcript_id)

            transcript_ids[transcript_id] = set(ids)
        return transcript_ids, transcript_id_idx, transcript_gene_map
