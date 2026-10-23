from src.qa_browser.app import _get_build_status, _summarize_checkpoint


def _checkpoint():
    return {
        "run_id": "src/use_cases/pharos/pharos.yaml",
        "last_updated": "2026-10-23T14:00:00+00:00",
        "adapters": {
            "post:SetPreferredSymbolAdapter": {
                "status": "completed", "adapter_position": 1, "adapter_total": 2,
            },
            "input:ProteinAdapter": {
                "status": "completed", "adapter_position": 1, "adapter_total": 2,
            },
            "input:GeneAdapter": {
                "status": "completed", "adapter_position": 2, "adapter_total": 2,
            },
            "post:TDLInputAdapter": {
                "status": "running", "adapter_position": 2, "adapter_total": 2,
            },
        },
    }


def test_build_status_groups_adapter_phases_and_sums_phase_totals():
    summary = _summarize_checkpoint(_checkpoint())
    assert summary["known_total"] == 4
    assert [phase["label"] for phase in summary["phases"]] == ["Input adapters", "Post adapters"]
    assert [row["name"] for row in summary["adapters"]] == [
        "input:ProteinAdapter", "input:GeneAdapter",
        "post:SetPreferredSymbolAdapter", "post:TDLInputAdapter",
    ]
    assert [row["display_name"] for row in summary["adapters"]] == [
        "ProteinAdapter", "GeneAdapter", "SetPreferredSymbolAdapter", "TDLInputAdapter",
    ]


def test_build_status_waits_for_final_metadata_and_ignores_older_failure():
    current = _checkpoint()
    current["adapters"]["post:TDLInputAdapter"]["status"] = "completed"
    old = {
        "run_id": "src/use_cases/pharos/old_post.yaml",
        "last_updated": "2026-10-22T14:00:00+00:00",
        "adapters": {"post:OldAdapter": {"status": "failed", "adapter_total": 1}},
    }

    class AQL:
        metadata = None

        def execute(self, query):
            if "etl_checkpoint" in query:
                return [current, old]
            return [self.metadata]

    class DB:
        aql = AQL()

        def has_collection(self, name):
            return name == "metadata_store"

    db = DB()
    status = _get_build_status(db)
    assert status["overall_status"] == "post_processing"
    assert status["runs"][0]["status"] == "adapters_complete"
    assert status["failed"] == 0
    db.aql.metadata = {
        "run_date": "2026-10-23T14:01:00+00:00",
        "source_yaml": current["run_id"],
    }
    assert _get_build_status(db)["overall_status"] == "completed"
    assert _get_build_status(db)["runs"][0]["status"] == "completed"
