from src.use_cases.build_cli import run_arango_build_cli


def main():
    run_arango_build_cli(
        build_name="pharos_current_tdls_reconstruction",
        primary_yaml="./src/use_cases/pharos/pharos_current_tdls_reconstruction.yaml",
        post_yaml="./src/use_cases/pharos/pharos_current_tdls_reconstruction_aql_post.yaml",
    )


if __name__ == "__main__":
    main()
