from src.constants import Prefix
from src.interfaces.simple_enum import SimpleEnum


def test_prefix_parse_accepts_existing_spellings():
    def previous_parse(value):
        for member in Prefix:
            if member.value.lower() == value.strip().lower().replace('_', '.'):
                return member
            if member.value.lower() == value.strip().lower().replace('-', '.'):
                return member
            if member.value.lower() == value.strip().lower():
                return member

    for member in Prefix:
        for value in (
            member.value,
            member.value.upper(),
            member.value.replace('.', '_'),
            member.value.replace('.', '-'),
        ):
            assert Prefix.parse(value) is previous_parse(value)


def test_parse_preserves_first_matching_member_order():
    class Overlapping(SimpleEnum):
        Underscore = 'a_b'
        Dotted = 'a.b'

    assert Overlapping.parse('a_b') is Overlapping.Underscore
    assert Overlapping.parse('a-b') is Overlapping.Dotted
