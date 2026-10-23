import yaml


def test_target_graph_current_tdls_exports_read_nested_tdl_metadata():
    config_paths = [
        "./src/use_cases/pharos/target_graph.yaml",
        "./src/use_cases/pharos/impatient_target_graph_aql_post.yaml",
    ]

    for config_path in config_paths:
        with open(config_path, "r", encoding="utf-8") as handle:
            config = yaml.safe_load(handle)

        graph_views = {view["id"]: view for view in config["graph_views"]}
        query = graph_views["current_tdls"]["query"]

        assert "pro.tdl_meta.tdl_ligand_count" in query
        assert "pro.tdl_meta.tdl_drug_count" in query
        assert "pro.tdl_meta.tdl_go_term_count" in query
        assert "pro.tdl_meta.tdl_generif_count" in query
        assert "pro.tdl_meta.tdl_pm_score" in query
        assert "pro.tdl_meta.tdl_antibody_count" in query
        assert "pro.tdl_ligand_count" not in query


def test_active_pharos_graphs_keep_post_adapters_and_views_in_primary_yaml():
    expected = {
        "pharos.yaml": [
            "SetPreferredSymbolAdapter", "SetLigandActivityFlagAdapter",
            "SetGoTermLeafFlagAdapter", "TDLInputAdapter", "TDLOverrideAdapter",
        ],
        "impatient_pharos.yaml": [
            "SetPreferredSymbolAdapter", "SetLigandActivityFlagAdapter",
            "SetGoTermLeafFlagAdapter", "TDLInputAdapter", "TDLOverrideAdapter",
        ],
        "target_graph.yaml": [
            "SetPreferredSymbolAdapter", "ExpandIDGFamilies",
            "SetLigandActivityFlagAdapter", "SetGoTermLeafFlagAdapter",
            "TDLInputAdapter", "TDLOverrideAdapter",
        ],
    }
    for primary_name, post_classes in expected.items():
        with open(f"src/use_cases/pharos/{primary_name}", encoding="utf-8") as stream:
            primary = yaml.safe_load(stream)

        assert [item["class"] for item in primary["post_adapters"]] == post_classes
        assert [view["id"] for view in primary["graph_views"]] == ["current_tdls"]
        assert primary["output_adapters"][0]["kwargs"]["database_name"] == (
            primary["post_adapters"][0]["kwargs"]["database_name"]
        )
