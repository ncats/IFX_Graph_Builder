import copy
from abc import ABC, abstractmethod
from collections import Counter
from typing import List, Union, Dict, Generator, Any

from src.constants import DataSourceName
from src.interfaces.id_resolver import IdResolver
from src.models.datasource_version_info import DatasourceVersionInfo
from src.models.node import Node, Relationship
from src.shared.record_merger import FieldConflictBehavior


class InputAdapter(ABC):
    batch_size = 25000
    single_source: bool = False
    field_conflict_behavior: FieldConflictBehavior = FieldConflictBehavior.KeepFirst

    def set_unmatched_reporter(self, reporter):
        self._unmatched_reporter = reporter

    def _report_unmatched(self, node_type: str, role: str, counts: Counter, unmatched: dict):
        reporter = getattr(self, "_unmatched_reporter", None)
        if reporter and unmatched:
            reporter([
                {
                    "node_type": node_type,
                    "source_id": source_id,
                    "endpoint_role": role,
                    "occurrences": counts[source_id],
                }
                for source_id in unmatched
                if counts[source_id]
            ])

    @staticmethod
    def _canonicalize_relationship_class(
        rel: Relationship, start_node: Node, end_node: Node, *, already_copied: bool = False
    ) -> Relationship:
        from src.models.disease import GeneDiseaseEdge, ProteinDiseaseEdge
        from src.models.expression import GeneTissueExpressionEdge, ProteinTissueExpressionEdge
        from src.models.mouse_phenotype import GeneMousePhenotypeEdge, ProteinMousePhenotypeEdge
        from src.models.ortholog import GeneOrthologGeneEdge, ProteinOrthologGeneEdge
        from src.models.pathway import GenePathwayEdge, ProteinPathwayEdge
        from src.models.gwas_trait import GeneGwasTraitEdge, GwasTrait, ProteinGwasTraitEdge
        from src.models.protein import Protein
        from src.models.tissue import Tissue
        from src.models.disease import Disease
        from src.models.pathway import Pathway
        from src.models.mouse_phenotype import MousePhenotype
        from src.models.ortholog import OrthologGene

        mapping = {
            (GeneTissueExpressionEdge, Protein, Tissue): ProteinTissueExpressionEdge,
            (GeneDiseaseEdge, Protein, Disease): ProteinDiseaseEdge,
            (GenePathwayEdge, Protein, Pathway): ProteinPathwayEdge,
            (GeneOrthologGeneEdge, Protein, OrthologGene): ProteinOrthologGeneEdge,
            (GeneMousePhenotypeEdge, Protein, MousePhenotype): ProteinMousePhenotypeEdge,
            (GeneGwasTraitEdge, Protein, GwasTrait): ProteinGwasTraitEdge,
        }

        target_cls = mapping.get((type(rel), type(start_node), type(end_node)))
        if target_cls is None:
            return rel

        remapped = target_cls(
            start_node=start_node if already_copied else copy.deepcopy(start_node),
            end_node=end_node if already_copied else copy.deepcopy(end_node),
        )
        for key, value in rel.__dict__.items():
            if key in {"start_node", "end_node"}:
                continue
            setattr(remapped, key, value if already_copied else copy.deepcopy(value))
        if isinstance(rel, GeneDiseaseEdge) and target_cls is ProteinDiseaseEdge:
            # A scalar gene-pair assessment cannot represent multiple genes
            # after their edges merge on one protein-disease pair.
            remapped.disease_harmonizer_support = None
        return remapped

    def get_name(self) -> str:
        return f"{self.__class__.__name__} ({self.get_datasource_name().value})"

    def is_single_source(self) -> bool:
        return self.single_source

    def set_single_source(self, single_source: bool) -> 'InputAdapter':
        self.single_source = single_source
        return self

    def get_field_conflict_behavior(self) -> FieldConflictBehavior:
        return self.field_conflict_behavior

    @abstractmethod
    def get_all(self) -> Generator[List[Union[Node, Relationship]], None, None]:
        raise NotImplementedError("derived classes must implement get_all")

    def get_validators(self) -> list:
        return []

    def get_validation_data(self):
        return None

    @abstractmethod
    def get_datasource_name(self) -> DataSourceName:
        raise NotImplementedError("derived classes must implement get_datasource_name")

    @abstractmethod
    def get_version(self) -> DatasourceVersionInfo:
        raise NotImplementedError("derived classes must implement get_version")

    def get_resolved_and_provenanced_list(self, resolver_map: Dict[str, IdResolver]) -> Generator[list[Any], Any, None]:
        def get_and_delete_old_id(node):
            source_id = node.id
            if hasattr(node, 'old_id'):
                source_id = node.old_id
                delattr(node, 'old_id')
            return source_id

        for entries in self.get_all():

            for entry in entries:
                version_info = self.get_version()
                version_data = [self.get_datasource_name(), version_info.version, version_info.version_date, version_info.download_date]
                version_string = '\t'.join([str(e) for e in version_data])
                if not getattr(entry, 'provenance', None):
                    entry.provenance = version_string
                if self.get_datasource_name() != DataSourceName.PostProcessing and not getattr(entry, 'sources', None):
                    entry.sources = [version_string]

            nodes = [e for e in entries if isinstance(e, Node)]
            relationships = [e for e in entries if isinstance(e, Relationship)]

            if len(nodes) > 0:
                type_map = {}
                for node in nodes:
                    type = node.__class__.__name__
                    if type not in type_map:
                        type_map[type] = []
                    type_map[type].append(node)

                resolved_nodes = []
                for type, node_list in type_map.items():
                    if type in resolver_map:
                        resolver = resolver_map[type]
                        allow_retype = resolver.canonical_class is not None
                        entity_map = resolver.resolve_nodes(node_list, allow_retype=allow_retype)
                        self._report_unmatched(
                            type, "node", Counter(node.id for node in node_list),
                            entity_map[IdResolver.MatchKeys.unmatched],
                        )
                        resolved_nodes.extend(resolver.parse_flat_node_list_from_map(entity_map))
                    else:
                        resolved_nodes.extend(node_list)
                nodes = resolved_nodes

            for node in nodes:
                source_id = get_and_delete_old_id(node)
                node.entity_resolution = f"{self.get_datasource_name()}\t{self.__class__.__name__}\t{ source_id }"

            for i in range(0, len(nodes), self.batch_size):
                yield nodes[i:i + self.batch_size]

            for rel in relationships:
                start_source_id = get_and_delete_old_id(rel.start_node)
                end_source_id = get_and_delete_old_id(rel.end_node)
                rel.entity_resolution = f"{self.get_datasource_name()}\t{self.__class__.__name__}\t{ start_source_id }\t{ end_source_id }"

            if len(relationships) > 0:
                temp_nodes = [entry.start_node for entry in relationships] + [entry.end_node for entry in relationships]
                start_counts = Counter((entry.start_node.__class__.__name__, entry.start_node.id) for entry in relationships)
                end_counts = Counter((entry.end_node.__class__.__name__, entry.end_node.id) for entry in relationships)
                type_map = {}
                node_map = {}

                for node in temp_nodes:
                    type = node.__class__.__name__
                    if type not in type_map:
                        type_map[type] = []
                    type_map[type].append(node)

                for type in type_map:
                    if type in resolver_map:
                        resolver = resolver_map[type]
                        temp_node_map = resolver.resolve_nodes(type_map[type], allow_retype=True)
                        unmatched = temp_node_map[IdResolver.MatchKeys.unmatched]
                        self._report_unmatched(
                            type, "start", Counter({source_id: start_counts[(type, source_id)] for source_id in unmatched}),
                            unmatched,
                        )
                        self._report_unmatched(
                            type, "end", Counter({source_id: end_counts[(type, source_id)] for source_id in unmatched}),
                            unmatched,
                        )
                        parsed_node_map = resolver.parse_entity_map(temp_node_map)
                        for original_id, resolved_nodes in parsed_node_map.items():
                            node_map[(type, original_id)] = resolved_nodes

                return_relationships = []
                has_returned_batches = False
                for entry in relationships:
                    start_lookup = (entry.start_node.__class__.__name__, entry.start_node.id)
                    end_lookup = (entry.end_node.__class__.__name__, entry.end_node.id)

                    if entry.start_node.__class__.__name__ in resolver_map and start_lookup not in node_map:
                        continue
                    if entry.end_node.__class__.__name__ in resolver_map and end_lookup not in node_map:
                        continue

                    start_id = entry.start_node.id
                    start_nodes = [entry.start_node]

                    end_id = entry.end_node.id
                    end_nodes = [entry.end_node]

                    if start_lookup in node_map:
                        start_nodes = node_map[start_lookup]
                    if end_lookup in node_map:
                        end_nodes = node_map[end_lookup]

                    for start_node in start_nodes:
                        for end_node in end_nodes:
                            rel_copy = copy.deepcopy(entry)
                            if start_node.__class__ is not rel_copy.start_node.__class__:
                                rel_copy.start_node = copy.deepcopy(start_node)
                            else:
                                rel_copy.start_node.id = start_node.id
                            if end_node.__class__ is not rel_copy.end_node.__class__:
                                rel_copy.end_node = copy.deepcopy(end_node)
                            else:
                                rel_copy.end_node.id = end_node.id
                            rel_copy = self._canonicalize_relationship_class(
                                rel_copy, rel_copy.start_node, rel_copy.end_node, already_copied=True
                            )
                            return_relationships.append(rel_copy)

                    if len(return_relationships) >= self.batch_size:
                        has_returned_batches = True
                        print(f"prepared a batch of {len(return_relationships)} relationship records")
                        yield return_relationships
                        return_relationships = []

                if has_returned_batches:
                    print(f"final batch: {len(return_relationships)} relationship records")
                yield return_relationships
