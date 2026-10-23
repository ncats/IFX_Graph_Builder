from src.use_cases.build_cli import run_arango_build_cli


def main():
    run_arango_build_cli(
        build_name="pounce",
        primary_yaml="./src/use_cases/pounce/pounce.yaml",
        post_yaml="./src/use_cases/pounce/pounce_aql_post.yaml",
    )


if __name__ == "__main__":
    main()
