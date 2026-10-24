import yaml

from src.interfaces.resolver_metadata import resolver_fingerprints_by_type


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


def test_active_target_graphs_use_252_consumer_contract_resolvers():
    expected = {
        "pharos.yaml": ["TargetConsumerPharosCanonicalResolver"],
        "impatient_pharos.yaml": ["TargetConsumerPharosCanonicalResolver"],
        "target_graph.yaml": [
            "TargetConsumerGeneResolver",
            "TargetConsumerTranscriptResolver",
            "TargetConsumerProteinResolver",
        ],
    }
    for name, target_classes in expected.items():
        with open(f"src/use_cases/pharos/{name}", encoding="utf-8") as stream:
            config = yaml.safe_load(stream)
        resolvers = [
            entry for entry in config["resolvers"]
            if entry["class"] in target_classes
        ]
        assert [entry["class"] for entry in resolvers] == target_classes
        assert not any(
            entry["class"] in {
                "TCRDTargetResolver", "TargetGraphGeneResolver",
                "TargetGraphTranscriptResolver", "TargetGraphProteinResolver",
            }
            for entry in config["resolvers"]
        )
        assert all(
            entry["kwargs"]["data_source"]["snapshot_id"] == "ifx_harmonizers:targets:2.5.2"
            and entry["kwargs"]["no_match_behavior"] == "Skip"
            for entry in resolvers
        )
        if name != "target_graph.yaml":
            assert resolvers[0]["kwargs"]["multi_match_behavior"] == "All"
        else:
            assert all(
                entry["kwargs"]["multi_match_behavior"] == "Error"
                for entry in resolvers
            )


def test_tcrd_resolvers_match_source_pharos_graph():
    with open("src/use_cases/pharos/pharos.yaml", encoding="utf-8") as stream:
        pharos = yaml.safe_load(stream)
    with open("src/use_cases/pharos/tcrd.yaml", encoding="utf-8") as stream:
        tcrd = yaml.safe_load(stream)

    graph_resolvers = resolver_fingerprints_by_type(pharos["resolvers"])
    mysql_resolvers = resolver_fingerprints_by_type(tcrd["resolvers"])
    assert mysql_resolvers
    assert all(
        graph_resolvers[node_type]["fingerprint"] == expected["fingerprint"]
        for node_type, expected in mysql_resolvers.items()
    )
