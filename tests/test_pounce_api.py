"""The JSON API must run the same workbook checks as Excel uploads."""

from datetime import date, datetime
from pathlib import Path

from fastapi import FastAPI
from fastapi.testclient import TestClient

from src.input_adapters.excel_sheet_adapter import ExcelsheetParser
from src.input_adapters.pounce_sheets.pounce_parser import PounceParser
from src.qa_browser.pounce_api_contract import PounceValidationRequest
from src.qa_browser.pounce_api_limits import PounceJsonSizeLimit
from src.qa_browser.pounce_routes import _sheet_workbook_loader, router


_TEMPLATES = Path(__file__).resolve().parents[1] / "input_files/manual/pounce/OLD/templates"


def _workbook(path):
    sheets = []
    for name, frame in ExcelsheetParser(str(path)).sheet_dfs.items():
        rows = [
            [
                value.isoformat() if isinstance(value, (date, datetime))
                else value.item() if hasattr(value, "item") else value
                for value in row
            ]
            for row in frame.itertuples(index=False, name=None)
        ]
        sheets.append({"name": name, "columns": list(frame.columns), "rows": rows})
    return {"source_id": path.name, "sheets": sheets}


def _submission():
    return {
        "schema_version": 1,
        "project": _workbook(_TEMPLATES / "POUNCE_Project_Metabo_Template_v2.xlsx"),
        "experiments": [_workbook(_TEMPLATES / "POUNCE_Experiment_Metabo_Template_v2.xlsx")],
        "stats_results": [_workbook(_TEMPLATES / "POUNCE_StatsResults_Metabo_Template_v2.xlsx")],
    }


def _client():
    app = FastAPI()
    app.add_middleware(PounceJsonSizeLimit)
    app.include_router(router)
    return TestClient(app)


def test_json_sheet_tables_match_excel_structural_checks():
    request = PounceValidationRequest.model_validate(_submission())
    files = [
        _TEMPLATES / "POUNCE_Project_Metabo_Template_v2.xlsx",
        _TEMPLATES / "POUNCE_Experiment_Metabo_Template_v2.xlsx",
        _TEMPLATES / "POUNCE_StatsResults_Metabo_Template_v2.xlsx",
    ]
    excel_data, excel_issues = PounceParser().parse_all(
        str(files[0]), [str(files[1])], [str(files[2])]
    )
    json_data, json_issues = PounceParser().parse_all(
        request.project.source_id,
        [request.experiments[0].source_id],
        [request.stats_results[0].source_id],
        workbook_loader=_sheet_workbook_loader(
            [request.project, *request.experiments, *request.stats_results]
        ),
    )
    assert [issue.message for issue in json_issues] == [issue.message for issue in excel_issues]
    assert len(json_data.metabolites) == len(excel_data.metabolites)
    assert set(json_data.workbook_versions.values()) == {2}


def test_api_detects_mixed_versions_without_client_version_field():
    payload = _submission()
    project_map = next(sheet for sheet in payload["project"]["sheets"] if sheet["name"] == "BioSampleMap")
    project_map["columns"] = [
        "Submitter_Column_Name" if column == "Submitter_Variable_Name" else column
        for column in project_map["columns"]
    ]
    response = _client().post("/pounce/api/validate", json=payload)
    assert response.status_code == 200
    assert [file["template_version"] for file in response.json()["files"]] == [3, 2, 2]
    assert response.json()["error_count"] == len([
        issue for issue in response.json()["issues"] if issue["severity"] == "error"
    ])


def test_api_keeps_optional_stats_aligned_to_its_experiment():
    payload = _submission()
    first_experiment = payload["experiments"][0]
    second_experiment = _workbook(_TEMPLATES / "POUNCE_Experiment_Metabo_Template_v2.xlsx")
    second_experiment["source_id"] = "second-experiment.xlsx"
    payload["experiments"] = [first_experiment, second_experiment]
    payload["stats_results"] = [None, payload["stats_results"][0]]
    stats_meta = next(
        sheet for sheet in payload["stats_results"][1]["sheets"]
        if sheet["name"] == "StatsResultsMeta"
    )
    key_column = stats_meta["columns"].index("NCATSDPI_Variable_Name")
    value_column = stats_meta["columns"].index("Submitter_Value")
    for row in stats_meta["rows"]:
        if row[key_column] == "pre_processing_description":
            row[value_column] = ""
            break

    response = _client().post("/pounce/api/validate", json=payload)
    assert response.status_code == 200
    files = response.json()["files"]
    assert [item["source_id"] for item in files] == [
        payload["project"]["source_id"], first_experiment["source_id"],
        second_experiment["source_id"], payload["stats_results"][1]["source_id"],
    ]
    stats_issues = [
        issue for issue in response.json()["issues"] if issue["entity"] == "stats_results"
    ]
    assert stats_issues
    assert all(issue["source_id"] == payload["stats_results"][1]["source_id"] for issue in stats_issues)

    payload["stats_results"] = [None]
    assert _client().post("/pounce/api/validate", json=payload).status_code == 422


def test_api_accepts_no_stats_for_one_experiment():
    payload = _submission()
    payload["stats_results"] = [None]
    response = _client().post("/pounce/api/validate", json=payload)
    assert response.status_code == 200
    assert all(item["role"] != "stats_results" for item in response.json()["files"])


def test_api_rejects_missing_matrix_cells_and_duplicate_sources():
    payload = _submission()
    payload["experiments"][0]["source_id"] = payload["project"]["source_id"]
    response = _client().post("/pounce/api/validate", json=payload)
    assert response.status_code == 422
    assert "source_id must be unique" in str(response.json())

    payload = _submission()
    sheet = payload["experiments"][0]["sheets"][0]
    sheet["rows"][0].pop()
    response = _client().post("/pounce/api/validate", json=payload)
    assert response.status_code == 422
    assert "every row must have" in str(response.json())


def test_api_rejects_oversized_body_before_parsing(monkeypatch):
    import src.qa_browser.pounce_api_limits as limits

    monkeypatch.setattr(limits, "MAX_POUNCE_JSON_BYTES", 20)
    response = _client().post("/pounce/api/validate", json={"schema_version": 1, "project": {}})
    assert response.status_code == 413

    response = _client().post(
        "/pounce/api/validate",
        content=b"{" + b" " * 25 + b"}",
        headers={"content-length": "1", "content-type": "application/json"},
    )
    assert response.status_code == 413


def test_api_size_limit_honors_asgi_root_path(monkeypatch):
    import asyncio
    import src.qa_browser.pounce_api_limits as limits

    monkeypatch.setattr(limits, "MAX_POUNCE_JSON_BYTES", 20)
    statuses = []

    async def downstream(scope, receive, send):
        statuses.append("downstream")

    async def receive():
        return {"type": "http.request", "body": b"", "more_body": False}

    async def send(message):
        if message["type"] == "http.response.start":
            statuses.append(message["status"])

    asyncio.run(limits.PounceJsonSizeLimit(downstream)(
        {
            "type": "http", "path": "/qa/pounce/api/validate", "root_path": "/qa",
            "headers": [(b"content-length", b"21")],
        }, receive, send,
    ))
    assert statuses == [413]


def test_json_blank_cells_do_not_turn_numeric_ids_into_floats():
    import pandas as pd

    from src.input_adapters.pounce_sheets.pounce_parser import PounceParser

    frame = pd.DataFrame([[123], [""], [456]], columns=["id"])
    values = [PounceParser._get_row_value(row, "id") for _, row in frame.iterrows()]
    assert values == ["123", None, "456"]

    payload = _submission()
    sheet = next(
        item for item in payload["project"]["sheets"]
        if item["name"] == "BioSampleMeta"
    )
    sheet["columns"] = ["id"]
    sheet["rows"] = [[123], [None], [456]]
    request = PounceValidationRequest.model_validate(payload)
    loader = _sheet_workbook_loader([request.project, *request.experiments, *request.stats_results])
    workbook = loader(request.project.source_id, "project")
    values = [
        PounceParser._get_row_value(row, "id")
        for _, row in workbook.parser.sheet_dfs["BioSampleMeta"].iterrows()
    ]
    assert values == ["123", None, "456"]
