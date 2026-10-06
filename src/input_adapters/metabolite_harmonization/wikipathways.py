from pathlib import Path
from typing import Dict, Generator, Iterable, List, Optional, Set, Tuple, Union
import re
from urllib.parse import unquote
import zipfile

from rdflib import Graph, URIRef, RDF, RDFS
from rdflib.namespace import DC

from src.constants import DataSourceName
from src.interfaces.input_adapter import InputAdapter
from src.models.datasource_version_info import DatasourceVersionInfo
from src.models.metabolite_harmonization import (
    GeneIdentifier,
    GenePathwayDetail,
    GenePathwayEdge,
    MetaboliteIdentifier,
    MetaboliteIdentifierMappingDetail,
    MetaboliteIdentifierMappingEdge,
    MetaboliteName,
    MetabolitePathwayDetail,
    MetabolitePathwayEdge,
    PathwayIdentifier,
    PathwayName,
    ProteinIdentifier,
    ProteinPathwayDetail,
    ProteinPathwayEdge,
)


WP = "http://vocabularies.wikipathways.org/wp#"
PUBCHEM_COMPOUND_RDF_PREFIX = "http://rdf.ncbi.nlm.nih.gov/pubchem/compound/CID"
WIKIDATA_ENTITY_PREFIX = "http://www.wikidata.org/entity/"
IDENTIFIERS_ORG_PREFIX = "https://identifiers.org/"

IDENTIFIERS_ORG_PREFIXES = {
    "cas": "CAS",
    "chebi": "CHEBI",
    "chembl.compound": "ChEMBL.COMPOUND",
    "chemspider": "ChemSpider",
    "drugbank": "DRUGBANK",
    "hmdb": "HMDB",
    "kegg.compound": "KEGG.COMPOUND",
    "kegg.drug": "KEGG.DRUG",
    "kegg.glycan": "KEGG.GLYCAN",
    "knapsack": "KNApSAcK",
    "lipidmaps": "LIPIDMAPS",
    "pharmgkb.drug": "PharmGKB.DRUG",
    "pid.pathway": "PID.PATHWAY",
    "pubchem.compound": "PUBCHEM.COMPOUND",
    "pubchem.substance": "PUBCHEM.SUBSTANCE",
    "reactome": "Reactome",
    "wikidata": "Wikidata",
}

GENE_IDENTIFIER_PREFIXES = {
    "ncbigene": "NCBIGene",
    "ensembl": "Ensembl",
    "hgnc.symbol": "Symbol",
    "wikidata": "Wikidata",
}

PROTEIN_IDENTIFIER_PREFIXES = {
    "uniprot": "UniProtKB",
}


class WikiPathwaysMetaboliteEquivalenceAdapter(InputAdapter):
    def __init__(
        self,
        data_source=None,
        rdf_zip_file: Optional[str] = None,
        max_files: Optional[int] = None,
        max_records: Optional[int] = None,
    ):
        if data_source is not None:
            rdf_zip_file = str(data_source.file())
            self.version_info = data_source.version_info()
        else:
            self.version_info = DatasourceVersionInfo(
                version=None,
                version_date=None,
                download_date=None,
            )
        if rdf_zip_file is None:
            raise ValueError("WikiPathwaysMetaboliteEquivalenceAdapter requires data_source or rdf_zip_file")
        self.rdf_zip_file = Path(rdf_zip_file)
        self.max_files = max_files
        self.max_records = max_records

    def get_datasource_name(self) -> DataSourceName:
        return DataSourceName.WikiPathways

    def get_version(self) -> DatasourceVersionInfo:
        return self.version_info

    def get_all(self) -> Generator[
        List[Union[MetaboliteIdentifier, MetaboliteIdentifierMappingEdge]],
        None,
        None,
    ]:
        yield from self._iter_node_batches()
        yield from self._iter_edge_batches()

    def _iter_node_batches(self) -> Generator[List[MetaboliteIdentifier], None, None]:
        batch: List[MetaboliteIdentifier] = []
        emitted_ids: Set[str] = set()
        emitted_labeled_records: Set[Tuple[str, Tuple[str, ...]]] = set()
        emitted_xref_labels: Set[Tuple[str, Tuple[str, ...]]] = set()

        for record in self._iter_metabolite_records():
            source_id = record["source_id"]
            label_key = (source_id, tuple(record["labels"]))
            if label_key not in emitted_labeled_records:
                batch.append(self._source_node(record))
                emitted_labeled_records.add(label_key)
                emitted_ids.add(source_id)

            for _source_field, node_id in record["xrefs"]:
                if node_id == source_id:
                    continue
                labels = tuple(record["labels"])
                if labels:
                    xref_key = (node_id, labels)
                    if xref_key in emitted_xref_labels:
                        continue
                    emitted_xref_labels.add(xref_key)
                    names = [MetaboliteName(value=label, source="WikiPathways",
                                            source_field="rdfs:label") for label in labels]
                else:
                    if node_id in emitted_ids:
                        continue
                    names = []
                batch.append(MetaboliteIdentifier(id=node_id, names=names))
                emitted_ids.add(node_id)
                if len(batch) >= self.batch_size:
                    yield batch
                    batch = []

            if len(batch) >= self.batch_size:
                yield batch
                batch = []

        if batch:
            yield batch

    def _iter_edge_batches(self) -> Generator[List[MetaboliteIdentifierMappingEdge], None, None]:
        batch: List[MetaboliteIdentifierMappingEdge] = []
        emitted_details: Set[Tuple[str, str, str, str, str]] = set()

        for record in self._iter_metabolite_records():
            source_id = record["source_id"]
            for source_field, end_id in record["xrefs"]:
                if end_id == source_id:
                    continue
                detail_key = (source_id, end_id, "WikiPathways", source_field, source_id)
                if detail_key in emitted_details:
                    continue
                emitted_details.add(detail_key)
                batch.append(
                    MetaboliteIdentifierMappingEdge(
                        start_node=MetaboliteIdentifier(id=source_id),
                        end_node=MetaboliteIdentifier(id=end_id),
                        details=[
                            MetaboliteIdentifierMappingDetail(
                                source="WikiPathways",
                                source_field=source_field,
                                source_id=source_id,
                            )
                        ],
                    )
                )
                if len(batch) >= self.batch_size:
                    yield batch
                    batch = []

        if batch:
            yield batch

    def _iter_metabolite_records(self) -> Iterable[Dict]:
        count = 0
        with zipfile.ZipFile(self.rdf_zip_file) as archive:
            ttl_members = [
                info.filename
                for info in archive.infolist()
                if info.filename.endswith(".ttl")
            ]
            for file_idx, member in enumerate(ttl_members):
                if self.max_files is not None and file_idx >= self.max_files:
                    break
                graph = Graph()
                graph.parse(data=archive.read(member).decode("utf-8", "ignore"), format="turtle")
                if not self._is_human_pathway(graph):
                    continue
                for subject in graph.subjects(RDF.type, URIRef(WP + "Metabolite")):
                    record = self._parse_metabolite(graph, subject)
                    if record is None:
                        continue
                    yield record
                    count += 1
                    if self.max_records is not None and count >= self.max_records:
                        return

    @classmethod
    def _parse_metabolite(cls, graph: Graph, subject) -> Optional[Dict]:
        source_id = cls._normalize_uri(str(subject))
        if source_id is None:
            return None

        labels = cls._unique_clean_values(str(label) for label in graph.objects(subject, RDFS.label))
        xrefs = [("subject", source_id)]

        for predicate, obj in graph.predicate_objects(subject):
            predicate_uri = str(predicate)
            if not predicate_uri.startswith(WP + "bdb"):
                continue
            node_id = cls._normalize_uri(str(obj))
            if node_id is None:
                continue
            source_field = predicate_uri.rsplit("#", 1)[-1]
            xrefs.append((source_field, node_id))

        return {
            "source_id": source_id,
            "labels": labels,
            "xrefs": cls._unique_xrefs(xrefs),
        }

    @staticmethod
    def _is_human_pathway(graph: Graph) -> bool:
        organism = URIRef(WP + "organism")
        return any(str(obj).endswith("NCBITaxon_9606") for obj in graph.objects(None, organism))

    @staticmethod
    def _source_node(record: Dict) -> MetaboliteIdentifier:
        names = [
            MetaboliteName(
                value=label,
                source="WikiPathways",
                source_field="rdfs:label",
            )
            for label in record["labels"]
        ]
        return MetaboliteIdentifier(id=record["source_id"], names=names)

    @classmethod
    def _normalize_uri(cls, uri: str) -> Optional[str]:
        uri = uri.strip()
        if uri.startswith(PUBCHEM_COMPOUND_RDF_PREFIX):
            return cls._prefixed_id("PUBCHEM.COMPOUND", uri[len(PUBCHEM_COMPOUND_RDF_PREFIX):])
        if uri.startswith(WIKIDATA_ENTITY_PREFIX):
            return cls._prefixed_id("Wikidata", uri[len(WIKIDATA_ENTITY_PREFIX):])
        if uri.startswith(IDENTIFIERS_ORG_PREFIX):
            suffix = uri[len(IDENTIFIERS_ORG_PREFIX):]
            if "/" not in suffix:
                return None
            family, value = suffix.split("/", 1)
            prefix = IDENTIFIERS_ORG_PREFIXES.get(family)
            if prefix is None:
                return None
            return cls._prefixed_id(prefix, value)
        return None

    @staticmethod
    def _prefixed_id(prefix: str, value: str) -> Optional[str]:
        value = unquote(value).replace("\u00a0", " ").strip().rstrip("\u00c2").strip()
        if not value:
            return None
        if prefix == "CHEBI" and value.upper().startswith("CHEBI:"):
            value = value.split(":", 1)[1]
        if prefix == "PUBCHEM.COMPOUND":
            value = re.sub(r"^CID", "", value, flags=re.IGNORECASE)
        return f"{prefix}:{value}"

    @staticmethod
    def _unique_clean_values(values: Iterable[Optional[str]]) -> List[str]:
        result = []
        seen = set()
        for value in values:
            if value is None:
                continue
            value = value.strip()
            if not value or value in seen:
                continue
            seen.add(value)
            result.append(value)
        return result

    @staticmethod
    def _unique_xrefs(values: Iterable[Tuple[str, str]]) -> List[Tuple[str, str]]:
        result = []
        seen = set()
        for source_field, node_id in values:
            key = (source_field, node_id)
            if key in seen:
                continue
            seen.add(key)
            result.append(key)
        return result


class WikiPathwaysPathwayContextAdapter(InputAdapter):
    def __init__(
        self,
        data_source=None,
        rdf_zip_file: Optional[str] = None,
        max_files: Optional[int] = None,
    ):
        if data_source is not None:
            rdf_zip_file = str(data_source.file())
            self.version_info = data_source.version_info()
        else:
            self.version_info = DatasourceVersionInfo(version=None, version_date=None, download_date=None)
        if rdf_zip_file is None:
            raise ValueError("WikiPathwaysPathwayContextAdapter requires data_source or rdf_zip_file")
        self.rdf_zip_file = Path(rdf_zip_file)
        self.max_files = max_files
        self._pathway_records_cache: Optional[List[Dict]] = None

    def get_datasource_name(self) -> DataSourceName:
        return DataSourceName.WikiPathways

    def get_version(self) -> DatasourceVersionInfo:
        return self.version_info

    def get_all(self) -> Generator[
        List[Union[
            GeneIdentifier,
            MetaboliteIdentifier,
            PathwayIdentifier,
            ProteinIdentifier,
            GenePathwayEdge,
            MetabolitePathwayEdge,
            ProteinPathwayEdge,
        ]],
        None,
        None,
    ]:
        yield from self._iter_pathway_node_batches()
        yield from self._iter_metabolite_node_batches()
        yield from self._iter_gene_node_batches()
        yield from self._iter_protein_node_batches()
        yield from self._iter_metabolite_pathway_edge_batches()
        yield from self._iter_gene_pathway_edge_batches()
        yield from self._iter_protein_pathway_edge_batches()

    def _iter_pathway_node_batches(self) -> Generator[List[PathwayIdentifier], None, None]:
        batch: List[PathwayIdentifier] = []
        emitted: Set[str] = set()
        for record in self._iter_pathway_records():
            pathway_id = record["pathway_id"]
            if pathway_id in emitted:
                continue
            emitted.add(pathway_id)
            names = [
                PathwayName(value=record["title"], source="WikiPathways", source_field="dc:title")
            ] if record.get("title") else []
            batch.append(
                PathwayIdentifier(
                    id=pathway_id,
                    stable_id=record["wp_id"],
                    url=f"https://www.wikipathways.org/pathways/{record['wp_id']}.html",
                    names=names,
                )
            )
            if len(batch) >= self.batch_size:
                yield batch
                batch = []
        if batch:
            yield batch

    def _iter_metabolite_node_batches(self) -> Generator[List[MetaboliteIdentifier], None, None]:
        batch: List[MetaboliteIdentifier] = []
        emitted: Set[str] = set()
        for record in self._iter_pathway_records():
            for metabolite_id in record["metabolite_ids"]:
                if metabolite_id in emitted:
                    continue
                emitted.add(metabolite_id)
                batch.append(MetaboliteIdentifier(id=metabolite_id))
                if len(batch) >= self.batch_size:
                    yield batch
                    batch = []
        if batch:
            yield batch

    def _iter_gene_node_batches(self) -> Generator[List[GeneIdentifier], None, None]:
        batch: List[GeneIdentifier] = []
        emitted: Set[str] = set()
        for record in self._iter_pathway_records():
            for gene_id in record["gene_ids"]:
                if gene_id in emitted:
                    continue
                emitted.add(gene_id)
                batch.append(GeneIdentifier(id=gene_id))
                if len(batch) >= self.batch_size:
                    yield batch
                    batch = []
        if batch:
            yield batch

    def _iter_protein_node_batches(self) -> Generator[List[ProteinIdentifier], None, None]:
        batch: List[ProteinIdentifier] = []
        emitted: Set[str] = set()
        for record in self._iter_pathway_records():
            for protein_id in record["protein_ids"]:
                if protein_id in emitted:
                    continue
                emitted.add(protein_id)
                batch.append(ProteinIdentifier(id=protein_id))
                if len(batch) >= self.batch_size:
                    yield batch
                    batch = []
        if batch:
            yield batch

    def _iter_metabolite_pathway_edge_batches(self) -> Generator[List[MetabolitePathwayEdge], None, None]:
        batch: List[MetabolitePathwayEdge] = []
        emitted: Set[Tuple[str, str]] = set()
        for record in self._iter_pathway_records():
            for metabolite_id in record["metabolite_ids"]:
                edge_key = (metabolite_id, record["pathway_id"])
                if edge_key in emitted:
                    continue
                emitted.add(edge_key)
                batch.append(
                    MetabolitePathwayEdge(
                        start_node=MetaboliteIdentifier(id=metabolite_id),
                        end_node=PathwayIdentifier(id=record["pathway_id"]),
                        details=[
                            MetabolitePathwayDetail(
                                source="WikiPathways",
                                source_field="wp:Metabolite",
                                source_id=metabolite_id,
                                pathway_id=record["pathway_id"],
                                pathway_name=record.get("title"),
                                url=f"https://www.wikipathways.org/pathways/{record['wp_id']}.html",
                                species="Homo sapiens",
                            )
                        ],
                    )
                )
                if len(batch) >= self.batch_size:
                    yield batch
                    batch = []
        if batch:
            yield batch

    def _iter_gene_pathway_edge_batches(self) -> Generator[List[GenePathwayEdge], None, None]:
        batch: List[GenePathwayEdge] = []
        emitted: Set[Tuple[str, str]] = set()
        for record in self._iter_pathway_records():
            for gene_id in record["gene_ids"]:
                edge_key = (gene_id, record["pathway_id"])
                if edge_key in emitted:
                    continue
                emitted.add(edge_key)
                batch.append(
                    GenePathwayEdge(
                        start_node=GeneIdentifier(id=gene_id),
                        end_node=PathwayIdentifier(id=record["pathway_id"]),
                        details=[
                            GenePathwayDetail(
                                source="WikiPathways",
                                source_field="wp:GeneProduct",
                                gene_id=gene_id,
                                pathway_id=record["pathway_id"],
                                pathway_name=record.get("title"),
                                url=f"https://www.wikipathways.org/pathways/{record['wp_id']}.html",
                                species="Homo sapiens",
                            )
                        ],
                    )
                )
                if len(batch) >= self.batch_size:
                    yield batch
                    batch = []
        if batch:
            yield batch

    def _iter_protein_pathway_edge_batches(self) -> Generator[List[ProteinPathwayEdge], None, None]:
        batch: List[ProteinPathwayEdge] = []
        emitted: Set[Tuple[str, str]] = set()
        for record in self._iter_pathway_records():
            for protein_id in record["protein_ids"]:
                edge_key = (protein_id, record["pathway_id"])
                if edge_key in emitted:
                    continue
                emitted.add(edge_key)
                batch.append(
                    ProteinPathwayEdge(
                        start_node=ProteinIdentifier(id=protein_id),
                        end_node=PathwayIdentifier(id=record["pathway_id"]),
                        details=[
                            ProteinPathwayDetail(
                                source="WikiPathways",
                                source_field="wp:Protein",
                                protein_id=protein_id,
                                pathway_id=record["pathway_id"],
                                pathway_name=record.get("title"),
                                url=f"https://www.wikipathways.org/pathways/{record['wp_id']}.html",
                                species="Homo sapiens",
                            )
                        ],
                    )
                )
                if len(batch) >= self.batch_size:
                    yield batch
                    batch = []
        if batch:
            yield batch

    def _iter_pathway_records(self) -> Iterable[Dict]:
        if self._pathway_records_cache is not None:
            yield from self._pathway_records_cache
            return

        records: List[Dict] = []
        with zipfile.ZipFile(self.rdf_zip_file) as archive:
            ttl_members = [info.filename for info in archive.infolist() if info.filename.endswith(".ttl")]
            for file_idx, member in enumerate(ttl_members):
                if self.max_files is not None and file_idx >= self.max_files:
                    break
                graph = Graph()
                graph.parse(data=archive.read(member).decode("utf-8", "ignore"), format="turtle")
                if not self._is_human_pathway(graph):
                    continue
                wp_id = Path(member).stem
                records.append({
                    "wp_id": wp_id,
                    "pathway_id": f"WikiPathways:{wp_id}",
                    "title": self._pathway_title(graph),
                    "metabolite_ids": self._metabolite_ids(graph),
                    "gene_ids": self._gene_ids(graph),
                    "protein_ids": self._protein_ids(graph),
                })
        self._pathway_records_cache = records
        yield from records

    @staticmethod
    def _is_human_pathway(graph: Graph) -> bool:
        organism = URIRef(WP + "organism")
        return any(str(obj).endswith("NCBITaxon_9606") for obj in graph.objects(None, organism))

    @staticmethod
    def _pathway_title(graph: Graph) -> Optional[str]:
        for title in graph.objects(None, DC.title):
            value = str(title).strip()
            if value:
                return value
        return None

    @classmethod
    def _metabolite_ids(cls, graph: Graph) -> List[str]:
        ids: List[str] = []
        for subject in graph.subjects(RDF.type, URIRef(WP + "Metabolite")):
            ids.extend(cls._metabolite_ids_for_subject(graph, subject))
        return WikiPathwaysMetaboliteEquivalenceAdapter._unique_clean_values(ids)

    @classmethod
    def _metabolite_ids_for_subject(cls, graph: Graph, subject) -> List[str]:
        subject_id = WikiPathwaysMetaboliteEquivalenceAdapter._normalize_uri(str(subject))
        return [subject_id] if subject_id is not None else []

    @classmethod
    def _gene_ids(cls, graph: Graph) -> List[str]:
        ids = []
        for rdf_type in (URIRef(WP + "GeneProduct"), URIRef(WP + "Protein")):
            for subject in graph.subjects(RDF.type, rdf_type):
                ids.extend(cls._ids_for_subject(graph, subject, GENE_IDENTIFIER_PREFIXES))
        return WikiPathwaysMetaboliteEquivalenceAdapter._unique_clean_values(ids)

    @classmethod
    def _protein_ids(cls, graph: Graph) -> List[str]:
        ids = []
        for rdf_type in (URIRef(WP + "GeneProduct"), URIRef(WP + "Protein")):
            for subject in graph.subjects(RDF.type, rdf_type):
                ids.extend(cls._ids_for_subject(graph, subject, PROTEIN_IDENTIFIER_PREFIXES))
        return WikiPathwaysMetaboliteEquivalenceAdapter._unique_clean_values(ids)

    @classmethod
    def _ids_for_subject(cls, graph: Graph, subject, prefix_map: Dict[str, str]) -> List[str]:
        subject_id = cls._normalize_identifier_uri(str(subject), prefix_map)
        return [subject_id] if subject_id is not None else []

    @staticmethod
    def _normalize_identifier_uri(uri: str, prefix_map: Dict[str, str]) -> Optional[str]:
        uri = uri.strip()
        if uri.startswith(IDENTIFIERS_ORG_PREFIX):
            suffix = uri[len(IDENTIFIERS_ORG_PREFIX):]
            if "/" not in suffix:
                return None
            family, value = suffix.split("/", 1)
            prefix = prefix_map.get(family)
            if prefix is None:
                return None
            return WikiPathwaysMetaboliteEquivalenceAdapter._prefixed_id(prefix, value)
        if uri.startswith(WIKIDATA_ENTITY_PREFIX):
            prefix = prefix_map.get("wikidata")
            if prefix is None:
                return None
            return WikiPathwaysMetaboliteEquivalenceAdapter._prefixed_id(prefix, uri[len(WIKIDATA_ENTITY_PREFIX):])
        return None
