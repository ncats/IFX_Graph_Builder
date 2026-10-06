from dataclasses import dataclass
import re

from src.constants import Prefix
from src.models.go_term import GoTerm, GoType, ProteinGoTermEdge, GoEvidence
from src.models.keyword import Keyword
from src.models.pathway import Pathway
from src.models.protein import Protein
from src.models.disease import Disease, DiseaseAssociationDetail
from src.models.node import EquivalentId


@dataclass
class Alias:
    type: str
    term: str


class UniProtParser:
    PATHWAY_STEP_SUFFIX_PATTERN = re.compile(r'\s*:\s*step\s+\d+\s*/\s*\d+\.?\s*$', re.IGNORECASE)
    @staticmethod
    def find_matches(uniprot_obj, list_field, data_field, match_string):
        if list_field in uniprot_obj:
            return filter(lambda match: match[data_field] == match_string, uniprot_obj[list_field])
        return iter([])

    @staticmethod
    def find_comments(uniprot_obj, comment_type):
        return UniProtParser.find_matches(uniprot_obj, 'comments', 'commentType', comment_type)

    @staticmethod
    def find_first_comment(uniprot_obj, comment_type):
        first = next(UniProtParser.find_comments(uniprot_obj, comment_type), None)
        if first is not None and len(first) > 0:
            return first['texts'][0]['value']
        return None

    @staticmethod
    def find_cross_refs(uniprot_obj, ref_type):
        return UniProtParser.find_matches(uniprot_obj, 'uniProtKBCrossReferences', 'database', ref_type)

    @staticmethod
    def get_uniprot_id(uniprot_obj):
        return uniprot_obj['primaryAccession']

    @staticmethod
    def get_primary_accession(uniprot_obj):
        return f"{Prefix.UniProtKB}:" + UniProtParser.get_uniprot_id(uniprot_obj)

    @staticmethod
    def get_gene_name(uniprot_obj):
        return uniprot_obj['uniProtkbId']

    @staticmethod
    def get_full_name(uniprot_obj):
        protein_description = uniprot_obj.get('proteinDescription') or {}

        recommended_name = protein_description.get('recommendedName') or {}
        recommended_full_name = ((recommended_name.get('fullName') or {}).get('value'))
        if recommended_full_name:
            return recommended_full_name

        for submission_name in protein_description.get('submissionNames', []):
            submission_full_name = (submission_name.get('fullName') or {}).get('value')
            if submission_full_name:
                return submission_full_name

        for alternative_name in protein_description.get('alternativeNames', []):
            alternative_full_name = (alternative_name.get('fullName') or {}).get('value')
            if alternative_full_name:
                return alternative_full_name

        return None

    @staticmethod
    def get_description(uniprot_obj):
        return UniProtParser.find_first_comment(uniprot_obj, "FUNCTION")

    @staticmethod
    def get_similarity(uniprot_obj):
        return UniProtParser.find_first_comment(uniprot_obj, "SIMILARITY")

    @staticmethod
    def get_uniprot_entry_type(uniprot_obj):
        return uniprot_obj.get('entryType')

    @staticmethod
    def get_uniprot_reviewed(uniprot_obj):
        entry_type = UniProtParser.get_uniprot_entry_type(uniprot_obj)
        if entry_type is None:
            return None
        return entry_type == 'UniProtKB reviewed (Swiss-Prot)'

    @staticmethod
    def get_uniprot_annotation_score(uniprot_obj):
        score = uniprot_obj.get('annotationScore')
        if score is None:
            return None
        if isinstance(score, dict):
            return score.get('value')
        return score

    @staticmethod
    def get_symbols(uniprot_obj):
        if 'genes' not in uniprot_obj:
            return None
        symbols = []
        for g in uniprot_obj['genes']:
            if 'geneName' in g:
                symbols.append(g['geneName']['value'])
        if len(symbols) == 0:
            return None
        return symbols

    @staticmethod
    def get_primary_symbol(uniprot_obj):
        symbols = UniProtParser.get_symbols(uniprot_obj)
        if not symbols:
            return None
        return symbols[0]

    @staticmethod
    def get_sequence(uniprot_obj):
        return uniprot_obj['sequence']['value']

    @staticmethod
    def get_secondary_accessions(uniprot_obj):
        if 'secondaryAccessions' in uniprot_obj:
            return uniprot_obj['secondaryAccessions']
        return None

    @staticmethod
    def get_keywords(uniprot_obj):
        keywords = {}
        for keyword_obj in uniprot_obj.get('keywords', []):
            source_id = keyword_obj.get('id')
            value = keyword_obj.get('name')
            if not value:
                continue
            category = keyword_obj.get('category')
            keyword_id = f"keyword:uniprot:{(category or 'unknown').strip().lower()}:{value.strip().lower()}"
            if keyword_id not in keywords:
                keywords[keyword_id] = Keyword(
                    id=keyword_id,
                    source_id=source_id,
                    category=category,
                    source='UniProt',
                    value=value,
                )
        if len(keywords) == 0:
            return None
        return keywords

    @staticmethod
    def normalize_pathway_name(pathway_text):
        if pathway_text is None:
            return None
        normalized = UniProtParser.PATHWAY_STEP_SUFFIX_PATTERN.sub('', pathway_text).strip()
        if normalized.endswith('.'):
            normalized = normalized[:-1].strip()
        return normalized or None

    @staticmethod
    def get_pathways(uniprot_obj):
        pathways = {}
        for comment in UniProtParser.find_comments(uniprot_obj, 'PATHWAY'):
            for text_obj in comment.get('texts', []):
                pathway_text = text_obj.get('value')
                normalized_text = UniProtParser.normalize_pathway_name(pathway_text)
                if normalized_text is None:
                    continue
                for pathway_name in [piece.strip() for piece in normalized_text.split(';') if piece.strip()]:
                    pathway_id = f"Pathway:uniprot:{pathway_name.lower()}"
                    if pathway_id not in pathways:
                        pathways[pathway_id] = Pathway(id=pathway_id, source_id=pathway_name, type='UniProt', name=pathway_name)
        if len(pathways) == 0:
            return None
        return pathways

    @staticmethod
    def get_diseases(uniprot_obj):
        diseases = {}
        for comment in UniProtParser.find_comments(uniprot_obj, 'DISEASE'):
            disease_obj = comment.get('disease')
            if disease_obj is None:
                continue

            cross_ref = disease_obj.get('diseaseCrossReference') or {}
            database = cross_ref.get('database')
            identifier = cross_ref.get('id')
            if not database or not identifier:
                continue

            normalized_database = 'OMIM' if str(database).upper() == 'MIM' else database
            prefix = Prefix.parse(normalized_database)
            if prefix is None:
                disease_id = f"{normalized_database}:{identifier}"
            else:
                disease_id = EquivalentId(id=identifier, type=prefix).id_str()

            disease_name = disease_obj.get('diseaseId')
            description = disease_obj.get('description')
            evidences = disease_obj.get('evidences') or []
            evidence_codes = []
            pmids = []
            for evidence_obj in evidences:
                evidence_code = evidence_obj.get('evidenceCode')
                if evidence_code:
                    evidence_codes.append(evidence_code)
                if evidence_obj.get('source') == 'PubMed' and evidence_obj.get('id'):
                    pmids.append(str(evidence_obj['id']))
            if disease_id not in diseases:
                diseases[disease_id] = {
                    'disease': Disease(
                        id=disease_id,
                        name=disease_name,
                        uniprot_description=description,
                    ),
                    'detail': DiseaseAssociationDetail(
                        source='UniProt',
                        source_id=disease_id,
                        pmids=list(dict.fromkeys(pmids)),
                        evidence_codes=list(dict.fromkeys(evidence_codes)),
                    )
                }

        if len(diseases) == 0:
            return None
        return diseases

    @staticmethod
    def parse_aliases(uniprot_obj):
        aliases = []
        UniProtParser.append_to_list(aliases, Alias(type='primary accession', term=UniProtParser.get_primary_accession(uniprot_obj)))
        if 'secondaryAccessions' in uniprot_obj:
            for id_val in uniprot_obj['secondaryAccessions']:
                UniProtParser.append_to_list(aliases, Alias('secondary accession', f"{Prefix.UniProtKB}:" + id_val))
        UniProtParser.append_to_list(aliases, Alias('uniprot kb', uniprot_obj['uniProtkbId']))
        full_name = UniProtParser.get_full_name(uniprot_obj)
        if full_name:
            UniProtParser.append_to_list(aliases, Alias('full name', full_name))

        recommended_name = (uniprot_obj.get('proteinDescription') or {}).get('recommendedName') or {}
        if 'shortNames' in recommended_name:
            for obj in recommended_name['shortNames']:
                UniProtParser.append_to_list(aliases, Alias('short name', obj['value']))
        if 'genes' in uniprot_obj and len(uniprot_obj['genes']) > 0:
            for gene in uniprot_obj['genes']:
                if 'geneName' in gene:
                    UniProtParser.append_to_list(aliases, Alias('symbol', gene['geneName']['value']))
                if 'synonyms' in gene and len(gene['synonyms']) > 0:
                    for synonym in gene['synonyms']:
                        UniProtParser.append_to_list(aliases, Alias('synonym', synonym['value']))
        for match in UniProtParser.find_cross_refs(uniprot_obj, 'GeneID'):
            UniProtParser.append_to_list(aliases, Alias('NCBI Gene ID', f"{Prefix.NCBIGene}:{match['id']}"))
        ensembl_objs = UniProtParser.find_matches(uniprot_obj, 'uniProtKBCrossReferences', 'database', 'Ensembl')
        for match in ensembl_objs:
            UniProtParser.append_to_list(aliases, Alias('Ensembl', Prefix.ENSEMBL + ":" + UniProtParser.trim_version(match['id'])))
            if 'properties' in match:
                for prop in match['properties']:
                    UniProtParser.append_to_list(aliases, Alias('Ensembl', Prefix.ENSEMBL + ":" + UniProtParser.trim_version(prop['value'])))
        string_objs = UniProtParser.find_matches(uniprot_obj, 'uniProtKBCrossReferences', 'database', 'STRING')
        for match in string_objs:
            UniProtParser.append_to_list(aliases, Alias('STRING', UniProtParser.trim_species(match['id'])))
        refseq_objs = UniProtParser.find_matches(uniprot_obj, Prefix.ENSEMBL + ":" + 'uniProtKBCrossReferences', 'database', 'RefSeq')
        for match in refseq_objs:
            UniProtParser.append_to_list(aliases, Alias('RefSeq', Prefix.RefSeq + ":" + UniProtParser.trim_version(match['id'])))
        return aliases

    @staticmethod
    def append_to_list(aliases, new_alias):
        found = list(filter(lambda each: new_alias.type == each.type and new_alias.term == each.term, aliases))
        if len(found) == 0:
            aliases.append(new_alias)

    @staticmethod
    def trim_version(ensembl_id):
        return ensembl_id.split('.')[0]

    @staticmethod
    def trim_species(string_id):
        return string_id.split('.')[1]

    @staticmethod
    def parse_go_term(go_term_json):
        go_id = go_term_json['id']
        term = next(UniProtParser.find_matches(go_term_json, 'properties', 'key', 'GoTerm'))
        evidence = next(UniProtParser.find_matches(go_term_json, 'properties', 'key', 'GoEvidenceType'))
        pieces = term['value'].split(':')
        go_type = pieces[0]
        go_term = pieces[1]
        pieces = evidence['value'].split(':')
        eco_term = pieces[0]
        eco_assigned_by = pieces[1]

        return go_id, go_type, go_term, eco_term, eco_assigned_by

    @staticmethod
    def get_go_term_associations(uniprot_obj, protein: Protein):
        go_associations = []
        for go_term_json in UniProtParser.find_cross_refs(uniprot_obj, 'GO'):
            go_id, go_type, go_term, eco_term, eco_assigned_by = UniProtParser.parse_go_term(go_term_json)
            go_term = GoTerm(
                id= go_id,
                type= GoType.parse(go_type),
                term= go_term)
            go_associations.append(
                ProteinGoTermEdge(
                    start_node=protein,
                    end_node=go_term,
                    evidence=[GoEvidence(code=eco_term, assigned_by=eco_assigned_by)]
                )
            )
        return go_associations
