from src.use_cases.build_cli import run_arango_build_cli


if __name__ == "__main__":
    run_arango_build_cli(
        build_name="impatient_pharos",
        primary_yaml="src/use_cases/pharos/impatient_pharos.yaml",
    )

    run_arango_build_cli(
        build_name="target_graph",
        primary_yaml="src/use_cases/pharos/target_graph.yaml",
    )

    run_arango_build_cli(
        build_name="pharos",
        primary_yaml="src/use_cases/pharos/pharos.yaml",
    )
