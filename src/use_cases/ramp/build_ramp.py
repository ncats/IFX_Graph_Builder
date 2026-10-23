from src.use_cases.build_cli import run_arango_build_cli


def main():
    run_arango_build_cli(
        build_name="ramp",
        primary_yaml="./src/use_cases/ramp/ramp.yaml",
    )


if __name__ == "__main__":
    main()
