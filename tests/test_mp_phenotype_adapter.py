from datetime import date

from src.input_adapters.mp.phenotype_terms import MPPhenotypeAdapter
from tests.registry_fakes import registry_dataset


def test_mp_adapter_emits_mouse_phenotype_nodes(tmp_path):
    mp_path = tmp_path / "mp.obo"
    mp_path.write_text(
        """format-version: 1.2
data-version: releases/2026-05-01

[Term]
id: MP:0000001
name: mammalian phenotype

[Term]
id: MP:0000600
name: increased spleen weight

[Term]
id: GO:0008150
name: biological_process
""",
        encoding="utf-8",
    )

    adapter = MPPhenotypeAdapter(
        registry_dataset(
            tmp_path,
            mp_path.name,
            version="2026-05-01",
            version_date=date(2026, 5, 1),
        )
    )
    entries = [entry for batch in adapter.get_all() for entry in batch]

    assert len(entries) == 2
    assert {entry.id for entry in entries} == {"MP:0000001", "MP:0000600"}
    assert {entry.name for entry in entries} == {"mammalian phenotype", "increased spleen weight"}

    version = adapter.get_version()
    assert version.version == "2026-05-01"
    assert version.version_date.isoformat() == "2026-05-01"
