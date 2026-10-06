"""Explicit RaMP row mappings using supplied identity groups; no graph mutations."""
from collections import defaultdict
from functools import lru_cache
from src.use_cases.ramp.sqlite_source_evidence import SourceRows, evidence
import json
import math


PREFIXES = {
    "HMDB": "hmdb", "CHEBI": "chebi", "PUBCHEM.COMPOUND": "pubchem",
    "PubChem": "pubchem", "KEGG.COMPOUND": "kegg", "KEGG.GLYCAN": "kegg_glycan",
    "UniProtKB": "uniprot", "NCBIGene": "entrez", "NCBI.GENE": "entrez",
    "Entrez": "entrez", "Ensembl": "ensembl", "HGNC.SYMBOL": "gene_symbol",
    "Wikidata": "wikidata", "ChemSpider": "chemspider", "RefMet": "refmet", "REFMET": "refmet",
    "Symbol": "gene_symbol", "PlantFA": "plantfa",
    "RHEA": "rhea", "RHEA.COMP": "rhea-comp", "RHEA.POLYMER": "polymer",
    "LIPIDMAPS": "LIPIDMAPS", "CAS": "CAS", "SwissLipids": "swisslipids",
    "LipidBank": "lipidbank",
}
PROVIDERS = {"wikipathways": "wiki", "lipidmaps": "lipidmaps", "chebi": "chebi",
             "uniprotkb": "uniprot", "expasy": "expasy"}
NAME_PRIORITY = {"hmdb": 0, "chebi": 1, "lipidmaps": 2, "refmet": 3, "pubchem": 4}


def source_id(identifier):
    prefix, sep, accession = identifier.partition(":")
    if not sep or not accession:
        raise ValueError(f"Expected source CURIE: {identifier!r}")
    return f"{PREFIXES.get(prefix, prefix)}:{accession}"


def provider(value):
    value = value.split("\t", 1)[0].lower()
    return PROVIDERS.get(value, value)


def sources(doc):
    return sorted({provider(s) for s in doc.get("sources", []) if not s.startswith("Manual Curation\t")})


def number(value):
    if value in (None, ""):
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Nonfinite chemical property: {value}")
    return result


def name_candidates(doc):
    for name in doc.get("names") or []:
        if isinstance(name, str):
            yield (9, "", name)
        elif name.get("value"):
            src = provider(name.get("source") or "")
            rank = NAME_PRIORITY.get(src, 8)
            if src == 'lipidmaps':
                # The source SDF has no COMMON_NAME in this pin. Its concise
                # abbreviation is preferable to NAME's systematic-style text.
                rank += {'ABBREVIATION': 0, 'NAME': 0.1,
                         'COMMON_NAME': 0.2, 'SYSTEMATIC_NAME': 0.3}.get(
                             name.get('source_field'), 0.4)
            yield (rank, src, name["value"])
    for props in doc.get("chem_props") or []:
        if props.get("common_name"):
            src = provider(props["source"])
            rank = NAME_PRIORITY.get(src, 8) + (0.1 if src == 'lipidmaps' else 0)
            yield (rank, src, props["common_name"])
    if doc.get("name"):
        yield (9, "", doc["name"])
    yield (99, "", source_id(doc["id"]))


def hierarchy(reader, term_collection, edge_collection, *, parent_first=True):
    terms = {d["id"]: d for d in reader.records(term_collection)}
    parents = defaultdict(set)
    for e in reader.records(edge_collection):
        parent, child = (e["start_id"], e["end_id"]) if parent_first else (e["end_id"], e["start_id"])
        if parent not in terms or child not in terms:
            raise ValueError(f"Missing hierarchy term: {e}")
        parents[child].add(parent)

    @lru_cache(None)
    def ancestors(identifier):
        found, visiting = set(), set()

        def visit(node):
            if node in visiting:
                raise ValueError(f"Cycle in {term_collection}: {node}")
            if node not in terms:
                raise ValueError(f"Unknown {term_collection}: {node}")
            if node in found:
                return
            visiting.add(node)
            for p in sorted(parents[node]):
                visit(p)
            visiting.remove(node)
            found.add(node)
        visit(identifier)
        return tuple(sorted(found))
    return terms, ancestors


class Projection:
    def __init__(self, reader, writer, ontology_policy, gene_identity, progress=print, *, protein_annotations=None):
        from src.use_cases.ramp.sqlite_protein_annotations import ProteinAnnotations
        self.protein_annotations = protein_annotations or ProteinAnnotations()
        self.hmdb_status_manifest = {}
        self.reader, self.writer = reader, writer
        self.ontology_policy = ontology_policy
        self.progress = progress
        self.gene_identity = gene_identity
        self.metabolites = {}
        self.genes = {}
        self.proteins = {}
        self.pathways = {}
        self.pathway_sources = {}
        self.skipped = defaultdict(int)
        self.source_rows = SourceRows(writer)
        self.metabolite_source_groups = {}
        self.statuses = {}

    def run(self, *, source_only=False):
        if source_only:
            self.progress('Writing metabolite lookup identities')
            self.write_metabolites()
            self.progress('Writing gene/protein lookup identities')
            self.write_genes()
            self.progress('Writing association lookup identities')
            self.write_source_associations()
            self.progress('Finalizing source lookup rows')
            self.source_rows.finish()
            return
        for name, step in [("metabolites and chemistry", self.write_metabolites),
                           ("gene and protein identifiers", self.write_genes),
                           ("pathways and associations", self.write_pathways),
                           ("classifications and ontology", self.write_classifications),
                           ("reactions", self.write_reactions)]:
            self.progress(f"Writing {name}")
            step()
        self.source_rows.finish()

    def write_source_associations(self):
        """Read only the association evidence needed for source lookup rows."""
        ontology_terms, ontology_ancestors = hierarchy(
            self.reader, 'HmdbOntologyTerm', 'HmdbOntologyParentEdge')
        retained_ontology_terms = {tid for tid, term in ontology_terms.items()
                                   if term.get('term_type') == 'child'
                                   and term['name'] not in self.ontology_policy.get(term.get('ontology_type'), [])}
        for collection, kind, node_type in (
                ('MetabolitePathwayEdge', 'compound', None),
                ('HmdbMetaboliteOntologyEdge', 'compound', None),
                ('MetaboliteClassificationEdge', 'compound', None),
                ('HmdbMetaboliteProteinAssociationEdge', 'compound', None),
                ('RheaMetaboliteReactionEdge', 'compound', None),
                ('GenePathwayEdge', 'gene', 'GeneIdentifier'),
                ('ProteinPathwayEdge', 'gene', 'ProteinIdentifier'),
                ('RheaProteinReactionEdge', 'gene', 'ProteinIdentifier')):
            for edge in self.reader.records(collection):
                rid = (self.metabolites.get(edge['start_id']) if kind == 'compound'
                       else self.genes.get((node_type, edge['start_id'])))
                if rid is None:
                    if kind == 'gene':
                        raise ValueError(f'Missing {node_type}: {edge["start_id"]}')
                    self.skipped[collection] += 1
                    continue
                if collection == 'HmdbMetaboliteOntologyEdge' and not any(
                        tid in retained_ontology_terms for tid in ontology_ancestors(edge['end_id'])):
                    continue
                assertions = self.register_edge(collection, edge, rid, kind)
                if collection == 'HmdbMetaboliteProteinAssociationEdge':
                    gene = self.genes[('ProteinIdentifier', edge['end_id'])]
                    for src, _, detail in assertions:
                        accession = detail.get('hmdb_protein_accession')
                        if not accession:
                            raise ValueError('HMDB protein association lacks hmdb_protein_accession')
                        self.register_source('HMDB:' + accession, gene, 'gene', src)

    def register_source(self, raw, rid, kind, src, name=None, name_rank=0):
        identifier = source_id(raw)
        src = provider(src)
        known = self.metabolite_source_groups.get(identifier) if kind == 'compound' else None
        if known is not None and known != rid:
            raise ValueError(f'Source evidence {src}:{raw} contradicts metabolite groups {known} and {rid}')
        attribution = {'hmdb': 'hmdb_kegg', 'wiki': 'wikipathways_kegg'}.get(src, src) if identifier.startswith('kegg:') else src
        self.source_rows.add(identifier, rid, kind, attribution, name,
                             self.statuses.get(rid, 'no_HMDB_status'), name_rank)

    def register_edge(self, collection, edge, rid, kind):
        assertions = list(evidence(collection, edge))
        for src, raw, detail in assertions:
            self.register_source(raw, rid, kind, src,
                                 detail.get('name') if collection == 'RheaMetaboliteReactionEdge' else None)
            if kind == 'gene' and provider(src) == 'hmdb' and detail.get('hmdb_protein_accession'):
                self.register_source('HMDB:' + detail['hmdb_protein_accession'], rid, kind, src)
        return assertions

    def write_metabolites(self):
        groups = self.reader.groups()
        for i, members in enumerate(groups, 1):
            for member in members:
                self.metabolites[member] = f"RAMP_C_{i:09d}"
                normalized = source_id(member)
                previous = self.metabolite_source_groups.setdefault(normalized, self.metabolites[member])
                if previous != self.metabolites[member]:
                    raise ValueError(f'Normalized source identifier {normalized} spans metabolite groups')
        names, found = {}, set()
        missing_chem, groups_with_chem = set(), set()
        from src.use_cases.ramp.sqlite_hmdb_status import group_statuses
        statuses, self.hmdb_status_manifest = group_statuses(
            self.reader.records("MetaboliteIdentifier"), self.metabolites)
        self.statuses = statuses
        status_schema = 'hmdb_status' in getattr(self.reader, 'schemas', {}).get('MetaboliteIdentifier', {}).get('fields', {})
        self.hmdb_status_manifest['graph_schema_supports_status'] = status_schema
        if not status_schema:
            raise ValueError('MetaboliteIdentifier graph schema lacks hmdb_status; rebuild the graph and harmonization stage before export')
        for doc in self.reader.records("MetaboliteIdentifier"):
            identifier = doc["id"]
            if identifier not in self.metabolites:
                continue
            found.add(identifier)
            rid = self.metabolites[identifier]
            candidate = min(name_candidates(doc))
            if rid not in names or candidate < names[rid]:
                names[rid] = candidate
            self.write_identifier(doc, rid, "compound", candidate[2])
            if not doc.get("chem_props"):
                missing_chem.add(identifier)
            for props in doc.get("chem_props") or []:
                chem_id = props["source_id"]
                # Chemistry belongs to its reported identity, never another split group.
                if self.metabolites.get(chem_id) != rid:
                    raise ValueError(f"Chemistry identity {chem_id} disagrees with member {identifier}")
                self.writer.add("chem_props", ramp_id=rid,
                    chem_data_source=provider(props["source"]), chem_source_id=source_id(chem_id),
                    iso_smiles=props.get("iso_smiles"), inchi_key_prefix=props.get("inchi_key_prefix"),
                    inchi_key=props.get("inchi_key"), inchi=props.get("inchi"),
                    mw=number(props.get("mw")), monoisotop_mass=number(props.get("monoisotopic_mass")),
                    common_name=props.get("common_name"), mol_formula=props.get("molecular_formula"))
                self.register_source(chem_id, rid, "compound", props["source"],
                                     props.get("common_name"), name_rank=1)
                groups_with_chem.add(rid)
        if found != set(self.metabolites):
            raise ValueError(f"Missing active metabolite records: {sorted(set(self.metabolites)-found)[:5]}")
        self.write_chebi_bridge_sources_and_chemistry(missing_chem, groups_with_chem)
        for rid, (_, _, name) in sorted(names.items()):
            self.writer.add("analyte", rampId=rid, type="compound", common_name=name)

    def write_chebi_bridge_sources_and_chemistry(self, missing_chem, groups_with_chem):
        """Register bridged ChEBI IDs; use linked chemistry when the ID has none."""
        linked = {}
        active_bridges = set()
        for edge in self.reader.records("ChebiChemicalEntityMetaboliteIdentifierEdge"):
            identifier = edge["end_id"]
            if identifier not in self.metabolites:
                continue
            chemical_id = edge["start_id"]
            if chemical_id != identifier:
                raise ValueError(f"ChEBI chemistry link changes identifier: {chemical_id} -> {identifier}")
            if identifier in active_bridges:
                raise ValueError(f"Multiple ChEBI chemistry links for {identifier}")
            active_bridges.add(identifier)
            rid = self.metabolites[identifier]
            if identifier in missing_chem:
                linked[identifier] = chemical_id
        first_chem_groups = set()
        rows = 0
        for chemical in self.reader.records("ChemicalEntity"):
            identifier = chemical["id"]
            if linked.get(identifier) != identifier:
                continue
            if not any(chemical.get(field) not in (None, "") for field in
                       ("smiles", "inchi", "inchi_key", "formula", "mass", "monoisotopic_mass")):
                continue
            rid = self.metabolites[identifier]
            inchi_key = chemical.get("inchi_key")
            self.writer.add("chem_props", ramp_id=rid, chem_data_source="chebi",
                            chem_source_id=source_id(identifier), iso_smiles=chemical.get("smiles"),
                            inchi_key_prefix=inchi_key.split("-", 1)[0] if inchi_key else None,
                            inchi_key=inchi_key, inchi=chemical.get("inchi"),
                            mw=number(chemical.get("mass")),
                            monoisotop_mass=number(chemical.get("monoisotopic_mass")),
                            common_name=chemical.get("name"), mol_formula=chemical.get("formula"))
            self.register_source(identifier, rid, "compound", "chebi",
                                 chemical.get("name"), name_rank=1)
            rows += 1
            if rid not in groups_with_chem:
                first_chem_groups.add(rid)
        self.chebi_chemistry_fallback = {
            "policy": "Use same-CURIE ChemicalEntity fields only when active MetaboliteIdentifier.chem_props is empty",
            "source_lookup_policy": "Register a ChEBI bridge as ChEBI only when it supplies retained fallback chemistry",
            "active_bridge_identifiers": len(active_bridges),
            "added_rows": rows,
            "groups_gaining_first_chemistry": len(first_chem_groups),
        }

    def write_identifier(self, doc, rid, kind, name):
        for src in sources(doc):
            candidates = [n for n in name_candidates(doc) if n[1] == src or
                          (kind == 'gene' and not n[1] and n[0] < 99)]
            source_name = min(candidates)[2] if candidates else None
            # Preserve the legacy distinction between KEGG IDs reported by HMDB
            # and WikiPathways without pretending KEGG supplied those records.
            self.register_source(doc['id'], rid, kind, src, source_name)
        for item in (doc.get("names") or []) + (doc.get("synonyms") or []):
            value = item.get("value") if isinstance(item, dict) else item
            attribution = [provider(item["source"])] if isinstance(item, dict) and item.get("source") else sources(doc)
            if value:
                for src in attribution:
                    self.writer.add("analytesynonym", Synonym=value, rampId=rid, geneOrCompound=kind, source=src)

    def write_genes(self):
        groups = defaultdict(list)
        documents = []
        for collection in ("GeneIdentifier", "ProteinIdentifier"):
            for doc in self.reader.records(collection):
                documents.append((collection, doc))
        identities = self.gene_identity.groups((collection, doc['id']) for collection, doc in documents)
        for collection, doc in documents:
            groups[identities.input_to_group[(collection, doc['id'])]].append((collection, doc))
        for index, (key, members) in enumerate(sorted(groups.items()), 1):
            rid = f"RAMP_G_{index:09d}"
            names = []
            for collection, doc in sorted(members, key=lambda item: (item[0], item[1]['id'])):
                protein_name = self.protein_annotations.name(doc['id']) if collection == 'ProteinIdentifier' else None
                if collection == 'ProteinIdentifier' and not protein_name:
                    matches = self.gene_identity.matches_by_input[(collection, doc['id'])]
                    if len(matches) == 1:
                        protein_name = self.protein_annotations.name(next(iter(matches)))
                self.genes[(collection, doc["id"])] = rid
                name = doc.get("gene_name") or protein_name or min(name_candidates(doc))[2]
                names.append((not bool(doc.get('gene_name')), name, collection, doc['id']))
                self.write_identifier(doc, rid, "gene", name)
                if collection == "ProteinIdentifier":
                    if doc.get('hmdb_accession') and 'hmdb' in sources(doc):
                        self.register_source('HMDB:' + doc['hmdb_accession'], rid, 'gene', 'hmdb', doc.get('name'))
                    self.proteins[doc["id"]] = {k: doc.get(k) for k in ("name", "gene_name", "protein_type", "is_reviewed")}
                    if protein_name:
                        self.proteins[doc['id']]['name'] = protein_name
            self.writer.add("analyte", rampId=rid, type="gene", common_name=min(names)[1])
            for accession in identities.canonical_accessions[key]:
                # Preserve every matched accession, including gene-only groups.
                self.write_identifier({'id': accession, 'sources': ['uniprot'],
                                       'name': self.protein_annotations.name(accession)}, rid, 'gene', accession)
                for alias in self.protein_annotations.aliases.get(accession, ()):
                    self.register_source(alias, rid, 'gene', 'uniprot', self.protein_annotations.name(accession))

    def write_pathways(self):
        index = 0
        for doc in self.reader.records("PathwayIdentifier"):
            for src in sources(doc):
                index += 1
                rid = f"RAMP_P_{index:09d}"
                self.pathways[(doc["id"], src)] = rid
                attribution = 'kegg' if src == 'hmdb' and doc.get('category') == 'kegg' else src
                self.pathway_sources[rid] = attribution
                names = sorted(n["value"] for n in doc.get("names", []) if provider(n["source"]) == src and n.get("value"))
                self.writer.add("pathway", pathwayRampId=rid, sourceId=doc["id"].split(":", 1)[1],
                                type=attribution, pathwayCategory=doc.get("category"), pathwayName=names[0] if names else doc["id"])
        for collection, node_type in (("MetabolitePathwayEdge", None), ("GenePathwayEdge", "GeneIdentifier"),
                                      ("ProteinPathwayEdge", "ProteinIdentifier")):
            for edge in self.reader.records(collection):
                sid, pid = edge["start_id"], edge["end_id"]
                rid = self.genes.get((node_type, sid)) if node_type else self.metabolites.get(sid)
                if rid is None:
                    if node_type:
                        raise ValueError(f"Missing {node_type}: {sid}")
                    self.skipped[collection] += 1
                    continue
                assertions = self.register_edge(collection, edge, rid, 'gene' if node_type else 'compound')
                for src in sorted({provider(src) for src, _, _ in assertions}):
                    pathway = self.pathways.get((pid, src))
                    if not pathway:
                        raise ValueError(f"Missing provider pathway: {pid}, {src}")
                    self.writer.add("analytehaspathway", rampId=rid, pathwayRampId=pathway,
                                    pathwaySource=self.pathway_sources[pathway])
        catalyzed_types = defaultdict(set)
        for edge in self.reader.records("HmdbMetaboliteProteinAssociationEdge"):
            met = self.metabolites.get(edge["start_id"])
            if not met:
                self.skipped["HmdbMetaboliteProteinAssociationEdge"] += 1
                continue
            pid = edge["end_id"]
            gene = self.genes[("ProteinIdentifier", pid)]
            for src, raw, detail in self.register_edge('HmdbMetaboliteProteinAssociationEdge', edge, met, 'compound'):
                accession = detail.get('hmdb_protein_accession')
                if not accession:
                    raise ValueError('HMDB protein association lacks hmdb_protein_accession')
                self.register_source('HMDB:' + accession, gene, 'gene', src)
            catalyzed_types[(met, gene)].add(self.proteins[pid].get("protein_type") or "Unknown")
        for (met, gene), types in sorted(catalyzed_types.items()):
            self.writer.add("catalyzed", rampCompoundId=met, rampGeneId=gene,
                            proteinType='; '.join(sorted(types)))

    def write_classifications(self):
        terms, ancestors = hierarchy(self.reader, "MetaboliteClassificationTerm", "MetaboliteClassificationParentEdge")
        for edge in self.reader.records("MetaboliteClassificationEdge"):
            rid = self.metabolites.get(edge["start_id"])
            if not rid:
                self.skipped["MetaboliteClassificationEdge"] += 1
                continue
            for src, raw, detail in self.register_edge('MetaboliteClassificationEdge', edge, rid, 'compound'):
                for tid in ancestors(edge["end_id"]):
                    term = terms[tid]
                    self.writer.add("metabolite_class", ramp_id=rid, class_source_id=source_id(raw),
                        class_level_name=term["level_name"], class_name=term["name"], source=provider(src))
        terms, ancestors = hierarchy(self.reader, "HmdbOntologyTerm", "HmdbOntologyParentEdge")
        allowed = {tid: t for tid, t in terms.items() if t.get("term_type") == "child"
                   and t["name"] not in self.ontology_policy.get(t.get("ontology_type"), [])}
        ids = {tid: f"RAMP_OL_{i:09d}" for i, tid in enumerate(sorted(allowed), 1)}
        for tid, term in sorted(allowed.items()):
            self.writer.add("ontology", rampOntologyId=ids[tid], commonName=term["name"],
                            HMDBOntologyType=term["ontology_type"], metCount=None)
        for edge in self.reader.records("HmdbMetaboliteOntologyEdge"):
            rid = self.metabolites.get(edge["start_id"])
            if not rid:
                self.skipped["HmdbMetaboliteOntologyEdge"] += 1
                continue
            retained_ids = [tid for tid in ancestors(edge["end_id"]) if tid in ids]
            if not retained_ids:
                continue
            self.register_edge('HmdbMetaboliteOntologyEdge', edge, rid, 'compound')
            for tid in retained_ids:
                self.writer.add("analytehasontology", rampCompoundId=rid, rampOntologyId=ids[tid])

    def write_reactions(self):
        reactions = {d["id"]: d for d in self.reader.records("RheaReaction")}
        ids = {rid: f"RAMP_R_{i:09d}" for i, rid in enumerate(sorted(reactions), 1)}
        classes, ancestors = hierarchy(self.reader, "RheaReactionClass", "RheaReactionClassParentEdge", parent_first=False)
        ec = defaultdict(set)
        for edge in self.reader.records("RheaReactionReactionClassEdge"):
            rxn, term_id = edge["start_id"], edge["end_id"]
            ec[rxn].add(term_id.split(":", 1)[1])
            for tid in ancestors(term_id):
                term = classes[tid]
                lineage = sorted((classes[a] for a in ancestors(tid)), key=lambda d: (d["ec_level"], d["id"]))
                self.writer.add("reaction_ec_class", ramp_rxn_id=ids[rxn], rxn_source_id=source_id(rxn),
                    rxn_class_ec=tid.split(":", 1)[1], ec_level=term["ec_level"], rxn_class=term["name"],
                    rxn_class_hierarchy=" | ".join(t["name"] for t in lineage))
        for rid, doc in sorted(reactions.items()):
            status = {"Approved": 1, "Obsolete": -1, "Preliminary": 0}[doc["status"]]
            self.writer.add("reaction", ramp_rxn_id=ids[rid], rxn_source_id=source_id(rid), status=status,
                is_transport=int(doc["is_transport"]), direction=doc["direction"], label=doc.get("label") or "",
                equation=doc.get("equation") or "", html_equation=doc.get("html_equation") or "",
                ec_num="; ".join(sorted(ec[rid])), has_human_prot=-1, only_human_mets=-1)
        for edge in self.reader.records("RheaMetaboliteReactionEdge"):
            rid = self.metabolites.get(edge["start_id"])
            if not rid:
                self.skipped["RheaMetaboliteReactionEdge"] += 1
                continue
            rxn = edge["end_id"]
            _, raw, _ = self.register_edge('RheaMetaboliteReactionEdge', edge, rid, 'compound')[0]
            self.writer.add("reaction2met", ramp_rxn_id=ids[rxn], rxn_source_id=source_id(rxn), ramp_cmpd_id=rid,
                substrate_product={"left": 0, "right": 1}[edge["side"]], met_source_id=source_id(raw),
                met_name=edge.get("name"), is_cofactor=int(edge["is_cofactor"]) if edge.get("is_cofactor") is not None else -1)
        for edge in self.reader.records("RheaProteinReactionEdge"):
            pid, rxn = edge["start_id"], edge["end_id"]
            doc = self.proteins[pid]
            _, raw, _ = self.register_edge('RheaProteinReactionEdge', edge, self.genes[("ProteinIdentifier", pid)], 'gene')[0]
            self.writer.add("reaction2protein", ramp_rxn_id=ids[rxn], rxn_source_id=source_id(rxn),
                ramp_gene_id=self.genes[("ProteinIdentifier", pid)], uniprot=source_id(raw),
                protein_name=doc.get("name") or doc.get("gene_name") or source_id(pid),
                is_reviewed=int(doc["is_reviewed"]) if doc.get("is_reviewed") is not None else -1)


def write_versions(writer, reader, release_version, timestamp, *, kegg_via_hmdb=False):
    from src.use_cases.ramp.sqlite_version_metadata import display_metadata
    writer.add("db_version", ramp_version=release_version, load_timestamp=timestamp,
        version_notes="Base tables only; post-processing pending. See ramp_export_metadata.",
        met_intersects_json=None, gene_intersects_json=None, met_intersects_json_pw_mapped=None,
        gene_intersects_json_pw_mapped=None, db_sql_url=None)
    by_source = defaultdict(list)
    for dataset in reader.metadata["registry_datasets"]:
        by_source[dataset["source"]].append(dataset)
    for src, datasets in sorted(by_source.items()):
        writer.add("version_info", ramp_db_version=release_version, db_mod_date=timestamp,
            status="current", data_source_id=provider(src), **display_metadata(src, datasets),
            data_source_snapshot_ids=json.dumps(sorted({d["snapshot_id"] for d in datasets})))
    if kegg_via_hmdb:
        datasets = [d for d in by_source['hmdb'] if d['dataset'] == 'metabolites_xml']
        if not datasets:
            raise ValueError('KEGG pathways via HMDB require recorded HMDB metabolite input metadata')
        hmdb_version = display_metadata('hmdb', datasets)['data_source_version']
        writer.add('version_info', ramp_db_version=release_version, db_mod_date=timestamp,
                   status='current', data_source_id='kegg', data_source_name='KEGG',
                   data_source_url='https://www.genome.jp/kegg/',
                   data_source_version=f'From HMDB ({hmdb_version})',
                   data_source_snapshot_ids=json.dumps(sorted({d['snapshot_id'] for d in datasets})))
