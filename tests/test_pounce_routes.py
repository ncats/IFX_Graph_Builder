import io
import os
import tempfile

import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient
from fastapi.templating import Jinja2Templates
from starlette.datastructures import UploadFile

from src.qa_browser import pounce_routes
from src.qa_browser.pounce_routes import _append_new_pairs, _build_session, _summarize_edit_submission


def _upload(filename: str, content: bytes = b"test") -> UploadFile:
    return UploadFile(filename=filename, file=io.BytesIO(content))


def test_append_new_pairs_appends_complete_pairs():
    with tempfile.TemporaryDirectory() as tmp:
        exp_paths = []
        stats_paths = []

        _append_new_pairs(
            [_upload("exp.xlsx")],
            [_upload("stats.xlsx")],
            tmp,
            exp_paths,
            stats_paths,
        )

        assert len(exp_paths) == 1
        assert len(stats_paths) == 1
        assert os.path.basename(exp_paths[0]) == "exp_new0_exp.xlsx"
        assert os.path.basename(stats_paths[0]) == "stats_new0_stats.xlsx"


def test_append_new_pairs_allows_experiment_without_stats():
    with tempfile.TemporaryDirectory() as tmp:
        exp_paths = []
        stats_paths = []
        _append_new_pairs([_upload("exp.xlsx")], [], tmp, exp_paths, stats_paths)
        assert len(exp_paths) == 1
        assert stats_paths == [None]


def test_append_new_pairs_rejects_stats_without_experiment():
    with tempfile.TemporaryDirectory() as tmp:
        with pytest.raises(ValueError, match="no experiment workbook"):
            _append_new_pairs([_upload("")], [_upload("stats.xlsx")], tmp, [], [])


def test_build_session_tracks_all_filenames():
    with tempfile.TemporaryDirectory() as tmp:
        project_path = os.path.join(tmp, "project.xlsx")
        exp_path = os.path.join(tmp, "exp.xlsx")
        stats_path = os.path.join(tmp, "stats.xlsx")
        for path in [project_path, exp_path, stats_path]:
            with open(path, "wb") as handle:
                handle.write(b"x")

        session_id, session = _build_session(project_path, [exp_path], [stats_path])

        assert session_id
        assert session["file_basenames"] == ["project.xlsx", "exp.xlsx", "stats.xlsx"]
        assert session["file_paths"] == [project_path, exp_path, stats_path]


def test_build_session_preserves_empty_stats_slot_without_attachment():
    with tempfile.TemporaryDirectory() as tmp:
        project_path = os.path.join(tmp, "project.xlsx")
        exp_path = os.path.join(tmp, "exp.xlsx")
        for path in [project_path, exp_path]:
            with open(path, "wb") as handle:
                handle.write(b"x")
        _, session = _build_session(project_path, [exp_path], [None])
        assert session["stats_paths"] == [None]
        assert session["file_paths"] == [project_path, exp_path]
        assert session["file_basenames"] == ["project.xlsx", "exp.xlsx"]


def test_summarize_edit_submission_tracks_replaced_unchanged_and_added():
    prior_session = {
        "project_path": "/tmp/project.xlsx",
        "exp_paths": ["/tmp/exp1.xlsx"],
        "stats_paths": ["/tmp/stats1.xlsx"],
    }

    summary = _summarize_edit_submission(
        prior_session,
        _upload("project-new.xlsx"),
        [_upload("")],
        [_upload("stats1-new.xlsx")],
        [_upload("exp2.xlsx")],
        [_upload("stats2.xlsx")],
    )

    assert "project.xlsx -> project-new.xlsx" in summary["replaced"]
    assert "stats1.xlsx -> stats1-new.xlsx" in summary["replaced"]
    assert "exp1.xlsx" in summary["unchanged"]
    assert "exp2.xlsx" in summary["added"]
    assert "stats2.xlsx" in summary["added"]


def test_upload_and_edit_preserve_optional_stats_slots():
    from pathlib import Path
    import shutil

    root = Path(__file__).resolve().parents[1]
    templates = root / "input_files/manual/pounce/OLD/templates"
    project = (templates / "POUNCE_Project_Metabo_Template_v2.xlsx").read_bytes()
    experiment = (templates / "POUNCE_Experiment_Metabo_Template_v2.xlsx").read_bytes()
    stats = (templates / "POUNCE_StatsResults_Metabo_Template_v2.xlsx").read_bytes()
    app = FastAPI()
    app.include_router(pounce_routes.router)
    pounce_routes.set_templates(Jinja2Templates(directory=str(root / "src/qa_browser/templates")))
    client = TestClient(app)
    existing_sessions = set(pounce_routes._sessions)

    try:
        response = client.post("/pounce/validate", files=[
            ("project_file", ("project.xlsx", project)),
            ("experiment_files", ("first.xlsx", experiment)),
            ("stats_files", ("", b"")),
            ("experiment_files", ("second.xlsx", experiment)),
            ("stats_files", ("second-stats.xlsx", stats)),
        ])
        assert response.status_code == 200
        assert "Could Not Parse Files" not in response.text
        assert "No StatsResults workbook" in response.text
        session_id = next(iter(set(pounce_routes._sessions) - existing_sessions))
        session = pounce_routes._sessions[session_id]
        assert session["stats_paths"][0] is None
        assert session["stats_paths"][1] is not None
        assert len(session["file_paths"]) == 4

        response = client.post(f"/pounce/edit/{session_id}", data={"remove_stats_1": "on"}, files=[
            ("project_file", ("", b"")),
            ("experiment_files", ("", b"")),
            ("stats_files", ("", b"")),
            ("experiment_files", ("", b"")),
            ("stats_files", ("", b"")),
        ])
        assert response.status_code == 200
        assert "Could Not Parse Files" not in response.text
        new_session_id = next(iter(set(pounce_routes._sessions) - existing_sessions - {session_id}))
        edited = pounce_routes._sessions[new_session_id]
        assert edited["stats_paths"] == [None, None]
        assert len(edited["file_paths"]) == 3
        assert any(name.endswith("second-stats.xlsx") for name in edited["file_change_summary"]["removed"])
    finally:
        for sid in set(pounce_routes._sessions) - existing_sessions:
            shutil.rmtree(pounce_routes._sessions[sid]["dir"], ignore_errors=True)
            del pounce_routes._sessions[sid]


def test_existing_project_loader_accepts_experiment_without_stats(monkeypatch):
    import shutil

    class Result:
        def __init__(self, rows):
            self.rows = rows

        def mappings(self):
            return self

        def first(self):
            return self.rows[0]

        def all(self):
            return self.rows

    class Connection:
        calls = 0

        def execute(self, query, params):
            self.calls += 1
            if self.calls == 1:
                return Result([{
                    "id": "p1", "name": "Example", "original_filename": "project.xlsx",
                    "content_blob": b"project",
                }])
            return Result([{
                "experiment_id": "e1", "experiment_name": "Experiment",
                "experiment_filename": "experiment.xlsx", "experiment_blob": b"experiment",
                "stats_result_id": None, "stats_filename": None, "stats_blob": None,
            }])

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

    class Engine:
        def connect(self):
            return Connection()

    monkeypatch.setattr(pounce_routes, "_get_mysql_engine", lambda: Engine())
    session_id = pounce_routes._create_existing_project_session("p1")
    try:
        session = pounce_routes._sessions[session_id]
        assert session["stats_paths"] == [None]
        assert len(session["file_paths"]) == 2
    finally:
        shutil.rmtree(pounce_routes._sessions[session_id]["dir"], ignore_errors=True)
        del pounce_routes._sessions[session_id]
