import argparse
from core.experiment_controller import ExperimentController

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ZGP DFG Optimization Runner")
    parser.add_argument("--config", type=str, required=True, help="Ścieżka do pliku JSON z konfiguracją")
    parser.add_argument("--max_iterations", type=int, default=100, help="Maksymalna liczba iteracji")
    args = parser.parse_args()

    controller = ExperimentController(args.config)
    controller.run(max_iterations=args.max_iterations)