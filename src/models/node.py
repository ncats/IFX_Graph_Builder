from dataclasses import dataclass, field, fields, is_dataclass
from datetime import datetime, date
from typing import List, Union, get_origin, get_args, Optional

from src.constants import Prefix
from src.core.decorators import indexed, facets, search


def is_atypeof_field(type_hint, cls):
    if type_hint is cls:
        return True
    origin = get_origin(type_hint)
    args = get_args(type_hint)
    return origin in (Union, Optional) and cls in args

def is_datetime_field(type_hint):
    return is_atypeof_field(type_hint, datetime)

def is_date_field(type_hint):
    return is_atypeof_field(type_hint, date)

def unwrap_optional(type_hint):
    if get_origin(type_hint) is Union:
        args = get_args(type_hint)
        non_none = [arg for arg in args if arg is not type(None)]
        if len(non_none) == 1:
            return non_none[0]
    return type_hint

def generate_class_from_dict(cls, data: dict):
    result = {}
    for f in fields(cls):
        value = data.get(f.name)
        type_hint = unwrap_optional(f.type)

        if f.name == 'start_node':
            result[f.name] = generate_class_from_dict(type_hint, {'id': data['start_id']})
            continue

        if f.name == 'end_node':
            result[f.name] = generate_class_from_dict(type_hint, {'id': data['end_id']})
            continue

        if hasattr(type_hint, 'from_dict'):
            value = type_hint.from_dict(data)
            result[f.name] = None

        if value is None or (isinstance(value, str) and value.strip() == ''):
            result[f.name] = None
            continue

        if is_datetime_field(type_hint) and isinstance(value, str):
            try:
                value = datetime.fromisoformat(value)
            except ValueError:
                pass  # or raise a warning/log
        elif is_date_field(type_hint) and isinstance(value, str):
            try:
                value = date.fromisoformat(value)
            except ValueError:
                try:
                    value = datetime.fromisoformat(value).date()
                except ValueError:
                    pass

        elif get_origin(type_hint) is list and isinstance(value, list) and len(value) > 0:
            item_type = get_args(type_hint)[0]
            if isinstance(value[0], str) and hasattr(item_type, 'parse'):
                value = [item_type.parse(item) if isinstance(item, str) else item for item in value]
            elif hasattr(item_type, 'from_dict'):
                value = [item_type.from_dict(item) for item in value]
            elif isinstance(value[0], dict):
                value = [generate_class_from_dict(item_type, item) for item in value]

        elif isinstance(value, dict) and is_dataclass(type_hint):
            value = generate_class_from_dict(type_hint, value)
        elif isinstance(value, str) and hasattr(type_hint, 'parse'):
            value = type_hint.parse(value)

        result[f.name] = value

    return cls(**result)

@dataclass
class EquivalentId:
    id: str
    type: Prefix
    source: Optional[List[str]] = None
    status: Optional[str] = None

    def type_str(self):
        if hasattr(self.type, 'value'):
            return self.type.value
        return self.type
    def id_str(self):
        return f"{self.type_str()}:{self.id}"

    @staticmethod
    def parse(id: str):
        prefix_part, id_part = id.split(":", 1)
        if prefix_part == 'uniprot':
            prefix = Prefix.UniProtKB
        elif prefix_part == 'entrez':
            prefix = Prefix.NCBIGene
        elif prefix_part == 'gene_symbol':
            prefix = Prefix.Symbol
        elif prefix_part == 'ncbiprotein':
            prefix = Prefix.Unknown
        else:
            prefix = Prefix.parse(prefix_part)

        if prefix is None:
            print(f"unknown prefix: {prefix_part}")
            prefix = Prefix.Unknown
        return EquivalentId(id=id_part, type=prefix)

    def __hash__(self):
        # Convert the source list to a tuple for hashing
        return hash((self.id, self.type, tuple(self.source or []), self.status))

    def __eq__(self, other):
        if isinstance(other, EquivalentId):
            return (self.id == other.id and
                    self.type == other.type and
                    set(self.source or []) == set(other.source or []) and
                    self.status == other.status)
        return False


@dataclass
@indexed(fields=['id', 'xref'])
@facets(category_fields=['sources'])
@search(text_fields=['id'])
class Node:
    id: str
    xref: Optional[List[EquivalentId]] = field(default_factory=list)
    xref_match_types: dict[str, list[str]] = field(default_factory=dict, kw_only=True)
    provenance: Optional[str] = None
    sources: List[str] = field(default_factory=list)

    @classmethod
    def from_dict(cls, data: dict):
        return generate_class_from_dict(cls, data)


@dataclass
class Relationship:
    start_node: Node
    end_node: Node
    provenance: Optional[str] = None
    sources: List[str] = field(default_factory=list)

    @classmethod
    def from_dict(cls, data: dict):
        return generate_class_from_dict(cls, data)
