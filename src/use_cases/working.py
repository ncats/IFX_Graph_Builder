from src.use_cases.build_cli import run_arango_build_cli


if __name__ == "__main__":
    # run_arango_build_cli(
    #     build_name="TARGETS 2.5.2 Target Graph legacy ID files",
    #     primary_yaml="src/use_cases/working_targets_252_legacy.yaml",
    # )

    run_arango_build_cli(
        build_name="TARGETS 2.5.2 Target Graph contract",
        primary_yaml="src/use_cases/working_targets_252_contract.yaml",
    )

    run_arango_build_cli(
        build_name="TARGETS 2.5.2 Pharos legacy canonical",
        primary_yaml="src/use_cases/working_targets_252_pharos_legacy.yaml",
    )

    run_arango_build_cli(
        build_name="TARGETS 2.5.2 Pharos contract canonical",
        primary_yaml="src/use_cases/working_targets_252_pharos_contract.yaml",
    )
