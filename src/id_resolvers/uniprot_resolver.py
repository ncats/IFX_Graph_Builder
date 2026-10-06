import gzip
import os
from typing import List, Dict, Optional

import ijson

from src.constants import Prefix
from src.models.node import Node
from src.shared.uniprot_parser import UniProtParser
from src.interfaces.id_resolver import IdResolver, IdMatch

scores = {
    "exact": 0,
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


class UniProtResolver(IdResolver):
    name = "UniProt Resolver"
    alias_map: dict

    @staticmethod
    def sort_list(context_list: list):
        context_list.sort(key=lambda context: scores.get(context, float('inf')))
        return context_list

    @staticmethod
    def sort_matches(match_list: List[IdMatch]):
        match_list.sort(key=lambda x: scores.get(x.context[0], float('inf')))
        return match_list

    def __init__(self,
                 uniprot_json_path: Optional[str] = None,
                 data_source=None,
                 file_name: str = "uniprot-human-reviewed.json.gz",
                 types: Optional[List[str]] = None,
                 no_match_behavior="Allow",
                 multi_match_behavior="All"):
        self.alias_map = {}
        if data_source is not None:
            uniprot_json_path = str(data_source.file(file_name))
        if uniprot_json_path is None:
            raise ValueError("UniProtResolver requires uniprot_json_path or data_source")
        self.uniprot_json_path = uniprot_json_path
        IdResolver.__init__(
            self,
            types=types or ["Protein"],
            no_match_behavior=no_match_behavior,
            multi_match_behavior=multi_match_behavior,
        )
        path = os.path.expanduser(self.uniprot_json_path)
        print(f"reading file at {path}")
        with gzip.open(path, 'rb') as gzip_file:
            for entry in ijson.items(gzip_file, 'results.item'):
                primary_accession = UniProtParser.get_primary_accession(entry)
                aliases = UniProtParser.parse_aliases(entry)
                for alias in aliases:
                    if alias.term not in self.alias_map:
                        self.alias_map[alias.term] = [
                            IdMatch(input=alias.term, match=primary_accession, context=[alias.type])
                        ]
                    else:
                        list = self.alias_map[alias.term]
                        same_id_match = [matching_entry for matching_entry in list if
                                         matching_entry.match == primary_accession]
                        if len(same_id_match) > 0:
                            same_id_match[0].context.append(alias.type)
                        else:
                            list.append(
                                IdMatch(input=alias.term, match=primary_accession, context=[alias.type])
                            )
                for go_term_json in UniProtParser.find_cross_refs(entry, 'GO'):
                    go_id, go_type, go_term, eco_term, eco_assigned_by = UniProtParser.parse_go_term(go_term_json)
                    self.alias_map[go_id] = [
                        IdMatch(input=go_id, match=go_id, context=['exact'])
                    ]

    def get_matches_for_merged_list(self, input_ids: List[str]):
        id_list = list(set(input_ids))
        unsorted_matches = []
        for input_id in id_list:
            individual_matches = self.alias_map.get(input_id)
            if individual_matches:
                unsorted_matches.extend(self.alias_map.get(input_id))
        return self.return_matches(unsorted_matches)

    def return_matches(self, unsorted_matches):
        for match in unsorted_matches:
            match.context = UniProtResolver.sort_list(match.context)
        sorted_matches = UniProtResolver.sort_matches(unsorted_matches)
        if len(sorted_matches) == 0:
            return []
        best_context = sorted_matches[0].context[0]
        best_matches = [match for match in sorted_matches if match.context[0] == best_context]
        return best_matches

    def _resolve(self, input_id):
        if input_id not in self.alias_map:
            return []
        unsorted_matches = self.alias_map.get(input_id)
        for match in unsorted_matches:
            match.context = UniProtResolver.sort_list(match.context)
        return self.return_matches(unsorted_matches)

    def resolve_internal(self, input_nodes: List[Node]) -> Dict[str, List[IdMatch]]:
        result_list = {}
        for node in input_nodes:
            lookup_id = node.id
            if lookup_id.startswith(('Symbol:', 'HGNC.SYMBOL:')):
                # Restrict explicit symbols to symbol/synonym aliases, not protein names.
                aliases = self.alias_map.get(lookup_id.split(':', 1)[1], [])
                best_matches = self.return_matches([
                    IdMatch(input=node.id, match=m.match,
                            context=[c for c in m.context if c in ('symbol', 'synonym')])
                    for m in aliases if any(c in ('symbol', 'synonym') for c in m.context)
                ])
            else:
                if lookup_id.split(':', 1)[0].casefold() == str(Prefix.ENSEMBL).casefold():
                    lookup_id = self.clean_id(lookup_id, 'ensembl')
                best_matches = self.get_matches_for_merged_list([lookup_id])
            if not best_matches and isinstance(node.id, str) and node.id.startswith(f"{Prefix.UniProtKB}:"):
                normalized_id = self.clean_id(node.id, 'uniprot')
                if normalized_id != node.id:
                    best_matches = self.get_matches_for_merged_list([normalized_id])
            result_list[node.id] = best_matches
        return result_list

    def clean_id(self, input_id: str, id_type: str):
        def remove_prefix(pre_id: str):
            pieces = pre_id.split(':')
            if len(pieces) > 1:
                return pieces[1]
            return pre_id

        def remove_isoform(pre_id: str):
            return pre_id.split('-')[0]

        def remove_version(pre_id: str):
            return pre_id.split('.')[0]

        id_type = id_type.lower()
        input_id = remove_prefix(input_id)
        if id_type == 'uniprot':
            input_id = remove_isoform(input_id)
            input_id = f"{Prefix.UniProtKB}:{input_id}"
        if id_type == 'ensembl':
            input_id = remove_version(input_id)
            input_id = f"{Prefix.ENSEMBL}:{input_id}"
        if id_type == 'brenda':
            input_id = f"{Prefix.BRENDA}:EC{input_id}"
        if id_type == 'chebi':
            input_id = f"{Prefix.CHEBI}:{input_id}"
        if id_type == 'entrez':
            input_id = f"{Prefix.NCBIGene}:{input_id}"
        if id_type == 'ncbiprotein':
            if id_type.isdigit():
                input_id = f"{Prefix.NCBIGene}:{input_id}"
            elif len(input_id) > 2 and input_id[0] == 'N' and input_id[2] == '_':
                input_id = remove_version(input_id)
                input_id = f"{Prefix.RefSeq}:{input_id}"
        if id_type == 'wikidata':
            input_id = f"{Prefix.Wikidata}:{input_id}"
        if id_type == 'hmdb':
            input_id = f"{Prefix.HMDB}:{input_id}"

        return input_id
