from enum import Enum
from functools import lru_cache
from typing import Union


@lru_cache(maxsize=None)
def _parse_members(enum_type):
    # Keep the first member for values that differ only in case.
    members = {}
    for index, member in enumerate(enum_type):
        members.setdefault(member.value.lower(), (index, member))
    return members


class SimpleEnum(str, Enum):
    def __str__(self):
        return self.value

    def __add__(self, other: str):
        return str.__add__(self.value, other)

    @classmethod
    def parse(cls, input_value: Union[str, Enum]):
        if input_value is None:
            return None
        if isinstance(input_value, Enum):
            return input_value
        value = input_value.strip().lower()
        members = _parse_members(cls)
        best = None
        for candidate_value in (value.replace('_', '.'), value.replace('-', '.'), value):
            candidate = members.get(candidate_value)
            if candidate is not None and (best is None or candidate[0] < best[0]):
                best = candidate
        if best is not None:
            return best[1]
        print(f"couldn't parse this value: {input_value}")


class LabeledIntEnum(int, Enum):
    """Enum whose .value is a stable integer ID (stored in DB) and .label is the display string.

    Use this when a controlled vocabulary needs a durable integer key in MySQL.
    The integer .value is serialized to ArangoDB and MySQL child tables.

    Usage::

        class MyType(LabeledIntEnum):
            FOO = (1, "Foo Display Name")
            BAR = (2, "Bar Display Name")
    """

    def __new__(cls, id: int, label: str):
        obj = int.__new__(cls, id)
        obj._value_ = id
        obj.label = label
        return obj

    def __str__(self):
        return self.label

    @classmethod
    def parse(cls, input_value: Union[str, int, 'LabeledIntEnum', None]):
        if input_value is None:
            return None
        if isinstance(input_value, cls):
            return input_value
        if isinstance(input_value, int):
            for member in cls:
                if member.value == input_value:
                    return member
        if isinstance(input_value, str):
            for member in cls:
                if member.label.lower() == input_value.strip().lower():
                    return member
        print(f"couldn't parse {cls.__name__} value: {input_value!r}")
        return None
