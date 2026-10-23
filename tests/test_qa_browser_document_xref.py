import asyncio

from starlette.requests import Request

import src.qa_browser.app as qa_app


def _request(path, query_string=b""):
    return Request({
        "type": "http",
        "method": "GET",
        "path": path,
        "root_path": "/odin-qa",
        "query_string": query_string,
        "headers": [],
    })


class _Collection:
    def __init__(self, document=None):
        self.document = document
        self.requested_keys = []

    def get(self, key):
        self.requested_keys.append(key)
        return self.document

    def properties(self):
        return {"type": "document"}


class _AQL:
    def __init__(self, matches):
        self.matches = matches
        self.calls = []

    def execute(self, query, *, bind_vars, max_runtime):
        self.calls.append((query, bind_vars, max_runtime))
        return self.matches


class _DB:
    def __init__(self, matches, document=None, collection_exists=True):
        self.coll = _Collection(document)
        self.aql = _AQL(matches)
        self.collection_exists = collection_exists

    def collection(self, name):
        assert name == "Gene"
        return self.coll

    def has_graph(self, name):
        return False

    def has_collection(self, name):
        return self.collection_exists


def test_symbol_xref_url_redirects_to_canonical_gene(monkeypatch):
    lookup = "Symbol:ABC:DEF"
    db = _DB([{"_key": "IFXGene:00015HD", "id": "IFXGene:00015HD"}])
    monkeypatch.setattr(qa_app, "get_db", lambda name: db)

    response = asyncio.run(qa_app.document_detail(
        _request(f"/db/test/collection/Gene/doc/{lookup}"), "test", "Gene", lookup,
    ))

    assert response.status_code == 303
    assert response.headers["location"] == (
        "/odin-qa/db/test/collection/Gene/doc/IFXGene%3A00015HD"
        "?resolved_from=Symbol%3AABC%3ADEF"
    )
    assert db.coll.requested_keys == [lookup]
    assert db.aql.calls[0][1] == {"@collection": "Gene", "lookup_id": lookup}


def test_ambiguous_xref_shows_canonical_choices(monkeypatch):
    lookup = "Symbol:DUPLICATE"
    db = _DB([
        {"_key": "IFXGene:A", "id": "IFXGene:A", "symbol": "DUPLICATE", "xref_match_types": {"symbol": [lookup]}},
        {"_key": "IFXGene:B", "id": "IFXGene:B", "symbol": "DUPLICATE", "xref_match_types": {"symbol": [lookup]}},
    ])
    monkeypatch.setattr(qa_app, "get_db", lambda name: db)

    response = asyncio.run(qa_app.document_detail(
        _request(f"/db/test/collection/Gene/doc/{lookup}"), "test", "Gene", lookup,
    ))

    body = response.body.decode()
    assert response.status_code == 200
    assert "IFXGene:A" in body and "IFXGene:B" in body
    assert "More than one Gene matches" in body


def test_direct_key_does_not_query_xrefs(monkeypatch):
    db = _DB([], {"_key": "IFXGene:A", "_id": "Gene/IFXGene:A"})
    monkeypatch.setattr(qa_app, "get_db", lambda name: db)

    response = asyncio.run(qa_app.document_detail(
        _request("/db/test/collection/Gene/doc/IFXGene:A"), "test", "Gene", "IFXGene:A",
    ))

    assert response.status_code == 200
    assert db.aql.calls == []


def test_missing_key_and_xref_keep_not_found_page(monkeypatch):
    db = _DB([])
    monkeypatch.setattr(qa_app, "get_db", lambda name: db)

    response = asyncio.run(qa_app.document_detail(
        _request("/db/test/collection/Gene/doc/Symbol:MISSING"),
        "test", "Gene", "Symbol:MISSING",
    ))

    assert response.status_code == 200
    assert "not found in" in response.body.decode()


def test_missing_collection_keeps_not_found_page(monkeypatch):
    db = _DB([], collection_exists=False)
    monkeypatch.setattr(qa_app, "get_db", lambda name: db)

    response = asyncio.run(qa_app.document_detail(
        _request("/db/test/collection/Gene/doc/Symbol:MISSING"),
        "test", "Gene", "Symbol:MISSING",
    ))

    assert "not found in" in response.body.decode()
    assert db.aql.calls == []


def test_resolved_from_banner_requires_actual_xref(monkeypatch):
    db = _DB([], {
        "_key": "IFXGene:A", "_id": "Gene/IFXGene:A", "id": "IFXGene:A",
        "xref": ["IFXGene:A", "Symbol:REAL"],
    })
    monkeypatch.setattr(qa_app, "get_db", lambda name: db)

    spoofed = asyncio.run(qa_app.document_detail(
        _request("/db/test/collection/Gene/doc/IFXGene:A", b"resolved_from=Symbol%3AFAKE"),
        "test", "Gene", "IFXGene:A",
    ))
    assert "Opened IFXGene:A from" not in spoofed.body.decode()

    real = asyncio.run(qa_app.document_detail(
        _request("/db/test/collection/Gene/doc/IFXGene:A", b"resolved_from=Symbol%3AREAL"),
        "test", "Gene", "IFXGene:A",
    ))
    assert "Opened IFXGene:A from Symbol:REAL" in real.body.decode()
