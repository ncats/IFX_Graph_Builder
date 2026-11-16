from src.use_cases.pharos.build_impatient_pharos import main as build_impatient_pharos
from src.use_cases.pharos.build_target_graph import main as build_target_graph
from src.use_cases.pharos.build_pharos import main as build_pharos


if __name__ == "__main__":
    build_impatient_pharos()
    build_target_graph()
    build_pharos()
