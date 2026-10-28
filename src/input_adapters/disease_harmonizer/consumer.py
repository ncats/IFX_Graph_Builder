"""Complete node payload from a pinned DISEASES release."""

from __future__ import annotations

from src.constants import DataSourceName
from src.interfaces.input_adapter import InputAdapter
from src.models.disease import Disease
from src.shared.disease_consumer_contract import disease_ids_file, iter_disease_rows


class DiseaseConsumerPayloadAdapter(InputAdapter):
    def __init__(self, data_source):
        self.path, _checksum = disease_ids_file(data_source)
        self.version_info = data_source.version_info()

    def get_datasource_name(self):
        return DataSourceName.DiseaseHarmonizer

    def get_version(self):
        return self.version_info

    def get_all(self):
        batch = []
        seen_ids = set()
        for row_number, row, _exact_aliases in iter_disease_rows(self.path):
            identifier = row["ncats_disease_id"].strip()
            if identifier in seen_ids:
                raise ValueError(f"{self.path}:{row_number} duplicates {identifier}")
            seen_ids.add(identifier)
            rare_value = row["is_rare"].strip().lower()
            if rare_value not in {"true", "false"}:
                raise ValueError(f"{self.path}:{row_number} has invalid is_rare {row['is_rare']!r}")
            batch.append(Disease(
                id=identifier,
                name=row["harmonized_name"].strip() or None,
                rare_disease=rare_value == "true",
                disease_harmonizer_payload=dict(row),
            ))
            if len(batch) >= self.batch_size:
                yield batch
                batch = []
        if not seen_ids:
            raise ValueError(f"{self.path} has no DISEASES rows")
        if batch:
            yield batch
