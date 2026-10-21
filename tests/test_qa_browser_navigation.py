from types import SimpleNamespace

from jinja2 import Environment, FileSystemLoader

from src.qa_browser.app import app


def _templates(root_path: str = "/odin-qa") -> Environment:
    environment = Environment(loader=FileSystemLoader("src/qa_browser/templates"))
    environment.globals["root_path"] = root_path
    return environment


def _request(path: str) -> SimpleNamespace:
    return SimpleNamespace(scope={"path": path})


def test_harmonization_home_lists_only_harmonizer_destinations():
    rendered = _templates().get_template("home.html").render(
        request=_request("/"),
        style_version="test",
    )

    assert rendered.count('class="db-card workflow-card') == 2
    assert "Metabolite Harmonization Studio" in rendered
    assert "Disease Harmonizer Explorer" not in rendered
    assert "Target Harmonizer Explorer" not in rendered
    assert "Variant Harmonizer Explorer" not in rendered
    assert "Drug Harmonizer Explorer" not in rendered
    assert "Pathway Harmonizer Explorer" not in rendered
    assert "CURE-ID Entity Resolver" in rendered
    assert 'href="/odin-qa/cure-entity-resolver"' in rendered
    assert "Pounce Submission" not in rendered
    assert "Database Browser" not in rendered
    assert 'href="/" aria-label="Return to the Informatics Development Sandbox"' in rendered


def test_standalone_home_adds_local_tools():
    rendered = _templates(root_path="").get_template("home.html").render(
        request=_request("/"),
        style_version="test",
    )

    assert rendered.count('class="db-card ') == 4
    assert "Local tools" in rendered
    assert 'href="/qa-browser"' in rendered
    assert 'href="/pounce/validate"' in rendered
    assert 'aria-label="Return to local tools">IFX ODIN</a>' in rendered
    assert "IFX Sandbox" not in rendered


def test_global_breadcrumb_separates_frontpage_tools_from_harmonization():
    template = _templates().get_template("base.html")

    pounce = template.render(request=_request("/pounce/validate"), style_version="test")
    harmonizer = template.render(request=_request("/ramp-id-qa"), style_version="test")

    assert '<a href="/">Home</a>' in pounce
    assert '<a href="/odin-qa/">Harmonization</a>' not in pounce
    assert '<a href="/odin-qa/">Harmonization</a>' in harmonizer


def test_standalone_frontpage_tools_link_back_to_local_tools():
    rendered = _templates(root_path="").get_template("base.html").render(
        request=_request("/qa-browser"),
        style_version="test",
    )

    assert '<a href="/">Local tools</a>' in rendered
    assert "IFX Sandbox" not in rendered


def test_migrated_harmonizer_routes_are_not_registered():
    paths = {route.path for route in app.routes}
    retired_prefixes = (
        "/target-id-qa",
        "/disease-id-qa",
        "/drug-id-qa",
        "/pathway-id-qa",
        "/variant-id-qa",
        "/api/v1/target",
        "/api/v1/disease",
        "/api/v1/drug",
    )

    assert not any(path.startswith(retired_prefixes) for path in paths)
    assert "/ramp-id-qa" in paths
    assert "/cure-entity-resolver" in paths
