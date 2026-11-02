import re
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, List, Optional


@dataclass
class ValidationError:
    severity: str        # "error" | "warning"
    entity: str          # e.g. "project", "biosamples"
    field: str           # e.g. "name", "project_type"
    message: str
    sheet: Optional[str] = None        # e.g. "ProjectMeta" — tells user where to look
    column: Optional[str] = None       # e.g. "project_name"
    row: Optional[int] = None          # set for list entities
    source_file: Optional[str] = None  # path to the xlsx file that produced this issue


class Validator(ABC):
    def __init__(self, entity: str, field: str, message: str,
                 sheet: Optional[str] = None, column: Optional[str] = None,
                 severity: str = "error"):
        self.entity = entity
        self.field = field
        self.message = message
        self.sheet = sheet
        self.column = column
        self.severity = severity

    def _make_error(self, row: Optional[int] = None) -> ValidationError:
        return ValidationError(
            severity=self.severity,
            entity=self.entity,
            field=self.field,
            message=self.message,
            sheet=self.sheet,
            column=self.column,
            row=row,
        )

    @abstractmethod
    def validate(self, data: Any) -> List[ValidationError]:
        raise NotImplementedError


class RequiredValidator(Validator):
    def validate(self, data: Any) -> List[ValidationError]:
        entity_obj = getattr(data, self.entity, None)
        if entity_obj is None:
            return [self._make_error()]
        if isinstance(entity_obj, list):
            errors = []
            for i, item in enumerate(entity_obj):
                v = getattr(item, self.field, None)
                if v is None or v == "" or v == []:
                    errors.append(self._make_error(row=i))
            return errors
        value = getattr(entity_obj, self.field, None)
        if value is None or value == "" or value == []:
            return [self._make_error()]
        return []


class ConditionalRequiredValidator(Validator):
    """Check that a field is present only when another field has a specific value.

    Useful when a field is required for some experiment types but not others,
    e.g. ``metabolite_identification_description`` is required only when
    ``platform_type`` is ``metabolomics`` or ``lipidomics``.

    Comparison is case-insensitive.
    """

    def __init__(self, entity: str, field: str, when_field: str,
                 when_values: List[str], message: str,
                 sheet: Optional[str] = None, column: Optional[str] = None):
        super().__init__(entity, field, message, sheet, column)
        self.when_field = when_field
        self.when_values = [str(v).lower() for v in when_values]

    def _condition_met(self, obj: Any) -> bool:
        val = getattr(obj, self.when_field, None)
        if val is None:
            return False
        if isinstance(val, list):
            return any(str(v).lower() in self.when_values for v in val)
        return str(val).lower() in self.when_values

    def validate(self, data: Any) -> List[ValidationError]:
        entity_obj = getattr(data, self.entity, None)
        if entity_obj is None:
            return []
        if isinstance(entity_obj, list):
            errors = []
            for i, item in enumerate(entity_obj):
                if self._condition_met(item):
                    v = getattr(item, self.field, None)
                    if v is None or v == "" or v == []:
                        errors.append(self._make_error(row=i))
            return errors
        if self._condition_met(entity_obj):
            value = getattr(entity_obj, self.field, None)
            if value is None or value == "" or value == []:
                return [self._make_error()]
        return []


class RequiredMapKeyValidator(Validator):
    """Checks that a key is configured in data.param_maps[sheet].

    Validates that the submitter has mapped the given NCATSDPI key to a column
    name in their map sheet — i.e., the configuration is present, regardless of
    whether the data rows have values.
    """

    def validate(self, data: Any) -> List[ValidationError]:
        param_maps = getattr(data, "param_maps", None) or {}
        if self.sheet not in param_maps:
            return []
        param_map = param_maps[self.sheet]
        val = param_map.get(self.field)
        if not val or val in ("NA", "N/A", ""):
            return [self._make_error()]
        return []


class ConditionalRequiredMapKeyValidator(Validator):
    """Check that a map-sheet key is configured only when a condition on another entity is met.

    Combines the map-key check of ``RequiredMapKeyValidator`` with the conditional
    logic of ``ConditionalRequiredValidator``.  Used when a map sheet's required
    keys depend on the experiment type — e.g. ``EffectSize_Map`` must have
    ``metabolite_id`` configured only for metabolomics/lipidomics runs.

    *when_entity* is the entity name on ``ParsedPounceData`` (e.g. ``"experiments"``).
    *when_field* is the Python field name on that entity.
    Comparison is case-insensitive.
    """

    def __init__(self, sheet: str, field: str, when_entity: str, when_field: str,
                 when_values: List[str], message: str, column: Optional[str] = None):
        super().__init__(entity="param_maps", field=field, message=message,
                         sheet=sheet, column=column)
        self.when_entity = when_entity
        self.when_field = when_field
        self.when_values = [str(v).lower() for v in when_values]

    def _condition_met(self, data: Any) -> bool:
        entity_obj = getattr(data, self.when_entity, None)
        if not entity_obj:
            return False
        if isinstance(entity_obj, list):
            return any(
                str(getattr(item, self.when_field, None) or "").lower() in self.when_values
                for item in entity_obj
            )
        val = getattr(entity_obj, self.when_field, None)
        return val is not None and str(val).lower() in self.when_values

    def validate(self, data: Any) -> List[ValidationError]:
        if not self._condition_met(data):
            return []
        param_maps = getattr(data, "param_maps", None) or {}
        if self.sheet not in param_maps:
            return [self._make_error()]
        val = param_maps[self.sheet].get(self.field)
        if not val or val in ("NA", "N/A", ""):
            return [self._make_error()]
        return []


class ParallelListsValidator(Validator):
    """Check that a group of list-valued fields on the same entity all have the same length.

    Used when two or more fields are pipe-delimited lists that get zipped together —
    e.g. ``owner_name`` and ``owner_email``, where each name must pair with one email.
    An empty or missing field is treated as length 0.
    """

    def __init__(self, entity: str, fields: List[str], columns: List[str],
                 sheet: Optional[str] = None,
                 optional_if_empty: Optional[List[str]] = None):
        col_str = ", ".join(columns)
        super().__init__(
            entity=entity,
            field=", ".join(fields),
            message=f"Fields [{col_str}] must have the same number of pipe-delimited entries",
            sheet=sheet,
            column=col_str,
        )
        self.fields = fields
        self.columns = columns
        self.optional_if_empty = set(optional_if_empty or [])

    def validate(self, data: Any) -> List[ValidationError]:
        entity_obj = getattr(data, self.entity, None)
        if entity_obj is None:
            return []
        if isinstance(entity_obj, list):
            errors = []
            for i, item in enumerate(entity_obj):
                err = self._check_lengths(item, row=i)
                if err:
                    errors.append(err)
            return errors
        err = self._check_lengths(entity_obj)
        return [err] if err else []

    def _check_lengths(self, obj: Any, row: Optional[int] = None) -> Optional[ValidationError]:
        lengths = {}
        for f in self.fields:
            val = getattr(obj, f, None)
            if val is None:
                lengths[f] = 0
            elif isinstance(val, list):
                lengths[f] = len(val)
            else:
                lengths[f] = 1

        compared_lengths = {
            f: length for f, length in lengths.items()
            if length or f not in self.optional_if_empty
        }
        if len(set(compared_lengths.values())) <= 1:
            return None

        detail = ", ".join(
            f"{col}={lengths[f]}"
            for f, col in zip(self.fields, self.columns)
        )
        return ValidationError(
            severity=self.severity,
            entity=self.entity,
            field=self.field,
            message=f"Mismatched list lengths: {detail}",
            sheet=self.sheet,
            column=self.column,
            row=row,
        )


class AllowedValuesValidator(Validator):
    def __init__(self, entity: str, field: str, message: str, values: List[str],
                 sheet: Optional[str] = None, column: Optional[str] = None):
        super().__init__(entity, field, message, sheet, column)
        self.values = values

    def _value_allowed(self, value: Any) -> bool:
        if isinstance(value, list):
            return all(v in self.values for v in value if v is not None)
        return value in self.values

    def _make_value_error(self, value: Any, row: Optional[int] = None) -> ValidationError:
        bad = [v for v in value if v not in self.values] if isinstance(value, list) else [value]
        bad_str = ", ".join(f"'{v}'" for v in bad)
        return ValidationError(
            severity=self.severity,
            entity=self.entity,
            field=self.field,
            message=f"{self.message} (got: {bad_str})",
            sheet=self.sheet,
            column=self.column,
            row=row,
        )

    def validate(self, data: Any) -> List[ValidationError]:
        entity_obj = getattr(data, self.entity, None)
        if entity_obj is None:
            return []
        if isinstance(entity_obj, list):
            errors = []
            for i, item in enumerate(entity_obj):
                value = getattr(item, self.field, None)
                if value is not None and not self._value_allowed(value):
                    errors.append(self._make_value_error(value, row=i))
            return errors
        value = getattr(entity_obj, self.field, None)
        if value is None:
            return []
        if not self._value_allowed(value):
            return [self._make_value_error(value)]
        return []


_NA_VALUES = {"NA", "N/A", ""}


class IndexedGroupValidator(Validator):
    """Check that indexed field groups in a map sheet are complete.

    If any key for a given index is configured in ``param_maps[sheet]``
    (e.g. ``exposure1_concentration``), all ``required_templates`` for that
    index must also be configured (e.g. ``exposure1_names``, ``exposure1_type``,
    ``exposure1_category``).

    ``required_templates`` is a list of key patterns containing ``{}`` as the
    index placeholder, e.g. ``["exposure{}_names", "exposure{}_type"]``.
    All templates must share the same prefix before ``{}``.
    """

    def __init__(self, sheet: str, required_templates: List[str]):
        prefix = required_templates[0].split("{}")[0]
        self._prefix = prefix
        self._required_templates = required_templates
        self._index_pattern = re.compile(r"^" + re.escape(prefix) + r"(\d+)")
        super().__init__(
            entity="param_maps",
            field=", ".join(required_templates),
            message=f"Incomplete indexed group '{prefix}{{}}*' in sheet '{sheet}'",
            sheet=sheet,
            column=", ".join(required_templates),
        )

    def validate(self, data: Any) -> List[ValidationError]:
        param_maps = getattr(data, "param_maps", None) or {}
        param_map = param_maps.get(self.sheet, {})

        indices = self._find_active_indices(param_map)
        errors = []
        for idx in sorted(indices):
            for template in self._required_templates:
                key = template.replace("{}", str(idx))
                val = param_map.get(key)
                if not val or val in _NA_VALUES:
                    errors.append(ValidationError(
                        severity=self.severity,
                        entity=self.entity,
                        field=template,
                        message=(
                            f"'{self._prefix}{idx}' group is configured but "
                            f"required key '{key}' is missing"
                        ),
                        sheet=self.sheet,
                        column=key,
                    ))
        return errors

    def _find_active_indices(self, param_map: dict) -> set:
        """Return the set of integer indices for which any matching key is configured."""
        indices = set()
        for key, val in param_map.items():
            if val and val not in _NA_VALUES:
                m = self._index_pattern.match(key)
                if m:
                    indices.add(int(m.group(1)))
        return indices
