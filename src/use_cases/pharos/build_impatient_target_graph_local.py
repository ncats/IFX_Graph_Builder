from src.use_cases.build_cli import run_arango_build_cli


def main():
    run_arango_build_cli(
        build_name="impatient_target_graph_local",
        primary_yaml="./src/use_cases/pharos/impatient_target_graph_local.yaml",
        post_yaml="./src/use_cases/pharos/impatient_target_graph_aql_post_local.yaml",
    )


if __name__ == "__main__":
    main()
