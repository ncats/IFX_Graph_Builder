import copy
from abc import ABC, abstractmethod
from dataclasses import fields
from enum import Enum
from typing import List, Dict, Optional, Type
from dataclasses import dataclass, asdict, field

from src.interfaces.simple_enum import SimpleEnum
from src.models.node import Node, EquivalentId


class NoMatchBehavior(SimpleEnum):
    Skip = "Skip"
    Allow = "Allow"
    Error = "Error"


class MultiMatchBehavior(SimpleEnum):
    First = "First"
    All = "All"
    Error = "Error"


@dataclass
class IdMatch:
    input: str
    match: str
    equivalent_ids: List[str] = field(default_factory=list)
    context: List[str] = None

    def to_dict(self) -> Dict[str, str]:
        return asdict(self)

class IdResolver(ABC):
    class MatchKeys(Enum):
        matched = "matched"
        newborns = "newborns"
        unmatched = "unmatched"

    name: str
    no_match_behavior: NoMatchBehavior
    multi_match_behavior: MultiMatchBehavior
    resolve_cache: Dict[str, List[IdMatch]]
    types: List[str]

    def __init__(self,
                 types: List[str],
                 no_match_behavior = NoMatchBehavior.Allow,
                 multi_match_behavior = MultiMatchBehavior.All,
                 canonical_class: Optional[Type[Node]] = None):
        print(f'creating ID resolver: {self.__class__.__name__}')
        self.types = types
        self.no_match_behavior = NoMatchBehavior.parse(no_match_behavior)
        self.multi_match_behavior = MultiMatchBehavior.parse(multi_match_behavior)
        self.resolve_cache = {}
        self.canonical_class = canonical_class

    @staticmethod
    def curie_prefix(value: str) -> Optional[str]:
        if not isinstance(value, str) or ":" not in value:
            return None
        prefix = value.split(":", 1)[0].strip()
        return prefix or None

    def get_prefix_counts(self) -> List[Dict[str, object]]:
        return []

    def get_example_ids(self, limit: int = 5) -> List[str]:
        return []

    def _resolve_internal(self, input_nodes: List[Node]) -> (Dict[str, List[IdMatch]], set):
        output_dict = {}
        un_resolved_nodes = []
        unique_nodes_dict = {node.id: node for node in input_nodes}
        unique_nodes = list(unique_nodes_dict.values())
        for node in unique_nodes:
            if node.id in self.resolve_cache:
                output_dict[node.id] = self.resolve_cache[node.id]
            else:
                un_resolved_nodes.append(node)
        if len(un_resolved_nodes) == 0:
            return output_dict

        new_nodes = self.resolve_internal(un_resolved_nodes)
        for node_id, match_info in new_nodes.items():
            output_dict[node_id] = match_info
            self.resolve_cache[node_id] = match_info

        return output_dict

    def resolve_nodes(self, entries: List, allow_retype: bool = False):
        id_map = self._resolve_internal(entries)
        return self.get_merged_map(entries, id_map, allow_retype=allow_retype)

    @staticmethod
    def _value_is_meaningful(value):
        if value is None:
            return False
        if isinstance(value, str):
            return value.strip() != ''
        if isinstance(value, (list, dict, set, tuple)):
            return len(value) > 0
        return True

    def _can_retype_entry(self, entry: Node) -> bool:
        if self.canonical_class is None:
            return False
        canonical_fields = {f.name for f in fields(self.canonical_class)}
        ignored_fields = {'id', 'xref', 'xref_match_types', 'provenance', 'sources'}
        for f in fields(entry.__class__):
            if f.name in ignored_fields:
                continue
            if not self._value_is_meaningful(getattr(entry, f.name, None)):
                continue
            if f.name not in canonical_fields:
                return False
        return True

    def _coerce_entry_to_canonical(self, entry: Node, new_id: str, full_xref_list: List[EquivalentId]):
        coerced = self.canonical_class(id=new_id)
        coerced.xref = full_xref_list
        canonical_fields = {f.name for f in fields(self.canonical_class)}
        for f in fields(entry.__class__):
            if f.name in {'id', 'xref', 'xref_match_types'}:
                continue
            if f.name not in canonical_fields:
                continue
            value = copy.deepcopy(getattr(entry, f.name, None))
            if self._value_is_meaningful(value):
                setattr(coerced, f.name, value)
        return coerced

    def get_merged_map(self, entries, id_map, allow_retype: bool = False):
        updated_count, validated_count = 0, 0
        entity_map = {
            IdResolver.MatchKeys.matched: {},
            IdResolver.MatchKeys.newborns: {},
            IdResolver.MatchKeys.unmatched: {}
        }
        for entry in entries:
            if entry.id in id_map and len(id_map.get(entry.id)) > 0:
                matches = id_map[entry.id]

                old_id = entry.id
                if matches:
                    if len(matches) > 0:
                        first_match = matches[0]
                        new_id = first_match.match

                        is_cross_type = self.canonical_class is not None and not isinstance(entry, self.canonical_class)
                        can_retype = self._can_retype_entry(entry) if is_cross_type else False
                        if is_cross_type and not allow_retype:
                            print(f"WARNING: {type(entry).__name__} '{old_id}' resolves cross-type to "
                                  f"{self.canonical_class.__name__} '{new_id}'. Treating as unmatched — "
                                  f"use a {self.canonical_class.__name__} adapter for node data, or "
                                  f"canonical_type is only applied for edge endpoint resolution.")
                            if old_id not in entity_map[IdResolver.MatchKeys.unmatched]:
                                new_entry = copy.deepcopy(entry)
                                new_entry.old_id = old_id
                                entity_map[IdResolver.MatchKeys.unmatched][old_id] = new_entry
                            continue
                        if is_cross_type and allow_retype and not can_retype:
                            print(f"WARNING: {type(entry).__name__} '{old_id}' resolves cross-type to "
                                  f"{self.canonical_class.__name__} '{new_id}', but populated fields are not "
                                  f"compatible with {self.canonical_class.__name__}. Treating as unmatched.")
                            if old_id not in entity_map[IdResolver.MatchKeys.unmatched]:
                                new_entry = copy.deepcopy(entry)
                                new_entry.old_id = old_id
                                entity_map[IdResolver.MatchKeys.unmatched][old_id] = new_entry
                            continue

                        needs_match = old_id not in entity_map[IdResolver.MatchKeys.matched]
                        needs_newborns = len(matches) > 1 and old_id not in entity_map[IdResolver.MatchKeys.newborns]
                        if not needs_match and not needs_newborns:
                            continue

                        full_xref_list = list(set(
                            EquivalentId.parse(equiv_id) for m in matches for equiv_id in m.equivalent_ids
                        ))

                        if needs_match:

                            if is_cross_type and allow_retype:
                                first_entry = self._coerce_entry_to_canonical(
                                    entry, new_id, self._xrefs_for_match(first_match, full_xref_list)
                                )
                            else:
                                first_entry = copy.deepcopy(entry)
                                first_entry.id = new_id
                                first_entry.xref = self._xrefs_for_match(first_match, full_xref_list)

                            if new_id != old_id:
                                first_entry.old_id = old_id
                                updated_count += 1
                            else:
                                validated_count += 1

                            entity_map[IdResolver.MatchKeys.matched][old_id] = first_entry

                        if len(matches) > 1:
                            if needs_newborns:

                                entity_map[IdResolver.MatchKeys.newborns][old_id] = []
                                for subsequent_match in matches[1:]:
                                    new_id = subsequent_match.match
                                    if is_cross_type and allow_retype:
                                        new_entry = self._coerce_entry_to_canonical(
                                            entry, new_id, self._xrefs_for_match(subsequent_match, full_xref_list)
                                        )
                                    else:
                                        new_entry = copy.deepcopy(entry)
                                        new_entry.id = new_id
                                        new_entry.xref = self._xrefs_for_match(subsequent_match, full_xref_list)
                                    new_entry.old_id = old_id
                                    if new_entry.id not in [node.id for node in
                                                            entity_map[IdResolver.MatchKeys.newborns][old_id]]:
                                        entity_map[IdResolver.MatchKeys.newborns][old_id].append(new_entry)
            else:
                if entry.id not in entity_map[IdResolver.MatchKeys.unmatched]:
                    new_entry = copy.deepcopy(entry)
                    new_entry.old_id = entry.id
                    entity_map[IdResolver.MatchKeys.unmatched][entry.id] = new_entry
        print(f"updated {updated_count} ids, validated {validated_count} ids")
        return entity_map

    def _xrefs_for_match(self, match: IdMatch, merged_xrefs: List[EquivalentId]) -> List[EquivalentId]:
        """Choose xrefs for one fan-out result; existing resolvers use the union."""
        return merged_xrefs

    def parse_entity_map(self, entity_map):
        match_map = {k: [v] for k, v in entity_map[IdResolver.MatchKeys.matched].items()}
        newborn_map = entity_map[IdResolver.MatchKeys.newborns]
        if self.multi_match_behavior == MultiMatchBehavior.All:
            for key, value in newborn_map.items():
                if key in match_map:
                    match_map[key].extend(value)
                else:
                    match_map[key] = value
        elif self.multi_match_behavior == MultiMatchBehavior.Error and newborn_map:
            source_id, additional = next(iter(newborn_map.items()))
            targets = [node.id for node in match_map.get(source_id, []) + additional]
            raise ValueError(
                f"Identifier {source_id!r} resolves to multiple targets {targets}; "
                "multi_match_behavior = 'Error'"
            )

        if self.no_match_behavior == NoMatchBehavior.Allow:
            unmatched_map = {}
            for k, v in entity_map[IdResolver.MatchKeys.unmatched].items():
                v.resolver_miss = True
                unmatched_map[k] = [v]
            for key, value in unmatched_map.items():
                if key in match_map:
                    match_map[key].extend(value)
                else:
                    match_map[key] = value

        return match_map

    def parse_flat_node_list_from_map(self, entity_map):
        matched_entries = [v for k,v in entity_map[IdResolver.MatchKeys.matched].items()]
        newborn_entries = [node for key, node_list in entity_map[IdResolver.MatchKeys.newborns].items() for node in node_list]
        unmatched_entries = [v for k,v in entity_map[IdResolver.MatchKeys.unmatched].items()]
        if self.multi_match_behavior == MultiMatchBehavior.All:
            matched_entries.extend(newborn_entries)
            print(f"unmerged {len(newborn_entries)} ids - multi_match_behavior = 'All'")
        elif self.multi_match_behavior == MultiMatchBehavior.First:
            print(f"skipping {len(newborn_entries)} unmergable ids - multi_match_behavior = 'First'")
        elif self.multi_match_behavior == MultiMatchBehavior.Error:
            if len(newborn_entries) > 0:
                raise Exception(
                    "node list has degenerate identifiers - multi_match_behavior = 'Error'",
                    newborn_entries)

        if self.no_match_behavior == NoMatchBehavior.Skip:
            print(f"skipping {len(unmatched_entries)} unmatched ids - no_match_behavior = 'Skip'")
        elif self.no_match_behavior == NoMatchBehavior.Allow:
            for entry in unmatched_entries:
                entry.resolver_miss = True
            matched_entries.extend(unmatched_entries)
            print(f"passing along {len(unmatched_entries)} unmatched ids - no_match_behavior = 'Allow'")
        elif self.no_match_behavior == NoMatchBehavior.Error:
            if len(unmatched_entries) > 0:
                raise Exception(
                    "node list has unmatched identifiers - no_match_behavior = 'Error'",
                    unmatched_entries
                )

        return matched_entries

    @abstractmethod
    def resolve_internal(self, input_nodes: List[Node]) -> Dict[str, List[IdMatch]]:
        raise NotImplementedError("derived classes must implement resolve_internal")
