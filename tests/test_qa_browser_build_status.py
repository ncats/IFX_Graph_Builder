from datetime import datetime

from src.qa_browser.app import _get_build_status, _get_build_diagnostics, _summarize_attempts, _summarize_checkpoint


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


def test_resumed_build_timing_separates_active_pause_and_total():
    timing = _summarize_attempts([
        {"id": "first", "status": "failed", "started_at": "2026-10-24T10:00:00+00:00",
         "ended_at": "2026-10-24T11:00:00+00:00", "elapsed_seconds": 3600},
        {"id": "second", "status": "completed", "started_at": "2026-10-24T13:00:00+00:00",
         "ended_at": "2026-10-24T16:00:00+00:00", "elapsed_seconds": 10800},
    ])
    assert timing["active"] == "4h 0m 0s"
    assert timing["paused"] == "2h 0m 0s"
    assert timing["total"] == "6h 0m 0s"
    assert timing["count"] == 2


def test_interrupted_attempt_does_not_claim_exact_active_time():
    timing = _summarize_attempts([
        {"status": "interrupted", "started_at": "2026-10-24T10:00:00+00:00"},
        {"status": "running", "started_at": "2026-10-24T13:00:00+00:00"},
    ], now=datetime.fromisoformat("2026-10-24T14:00:00+00:00").replace(tzinfo=None))
    assert timing["incomplete"] is True
    assert timing["active"] == ""
    assert timing["total"] == ""


def test_diagnostics_page_query_is_bounded_and_scoped_to_run():
    calls = []
    class AQL:
        def execute(self, query, bind_vars):
            calls.append((query, bind_vars))
            if "build_diagnostics_status" in query:
                return [{"status": "complete"}]
            if "COLLECT kind" in query:
                return [{"kind": "unmatched_id", "status": None, "count": 3}]
            return [{"_key": f"key{i}"} for i in range(101)]

    class DB:
        aql = AQL()
        def has_collection(self, name):
            return name in {"build_diagnostics", "metadata_store"}

    result = _get_build_diagnostics(DB(), "run.yaml", kind="unmatched_id", adapter="input:X")
    assert len(result["rows"]) == 100
    assert result["next_query"]
    assert all(call[1]["run_id"] == "run.yaml" for call in calls)
    assert calls[-1][1]["adapter"] == "input:X"


def test_failed_cleanup_attempt_overrides_completed_adapter_status():
    checkpoint = _checkpoint()
    for adapter in checkpoint["adapters"].values():
        adapter["status"] = "completed"
    checkpoint["attempts"] = [{
        "status": "failed", "started_at": "2026-10-23T12:00:00+00:00",
        "ended_at": "2026-10-23T14:00:00+00:00", "elapsed_seconds": 7200,
    }]

    class AQL:
        def execute(self, query):
            return [checkpoint] if "etl_checkpoint" in query else [{
                "run_date": "2026-10-23T14:01:00+00:00", "source_yaml": checkpoint["run_id"],
            }]

    class DB:
        aql = AQL()
        def has_collection(self, name):
            return name == "metadata_store"

    status = _get_build_status(DB())
    assert status["overall_status"] == "failed"
    assert status["runs"][0]["status"] == "failed"


def test_running_cleanup_attempt_keeps_page_auto_refresh_active():
    checkpoint = _checkpoint()
    for adapter in checkpoint["adapters"].values():
        adapter["status"] = "completed"
    checkpoint["attempts"] = [{"status": "running", "started_at": "2026-10-23T12:00:00+00:00"}]

    class AQL:
        def execute(self, query):
            return [checkpoint] if "etl_checkpoint" in query else [None]
    class DB:
        aql = AQL()
        def has_collection(self, name):
            return name == "metadata_store"

    status = _get_build_status(DB())
    assert status["overall_status"] == "running"
    assert status["is_active"] is True
