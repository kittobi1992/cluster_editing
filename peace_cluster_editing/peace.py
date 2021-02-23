#!/usr/bin/env python3

import argparse
from subprocess import run
from tempfile import NamedTemporaryFile
from pathlib import Path
from os import environ
from typing import Tuple, Dict, List, Optional
import numpy as np

BINARY_ENV_NAME = "PEACE_WEIGHTED_CLUSTER_EDITING_BINARY"


def convert_pace_to_cm(input_path: Path) -> str:
    with input_path.open() as input_file:
        lines = input_file.readlines()

    p_str, problem, n, m = lines[0].split()
    assert p_str == "p"
    assert problem == "cep"
    n, m = int(n), int(m)
    edges = [tuple(line.split()) for line in lines[1:]]
    edges = [(int(u) - 1, int(v) - 1) for u, v in edges]

    costs = -np.ones((n, n), dtype=float)
    for u, v in edges:
        costs[u, v] = 1
        costs[v, u] = 1

    instance_str = f"{n}\n"
    instance_str += "\n".join([f"{u}" for u in range(n)]) + "\n"
    for u in range(n - 1):
        instance_str += "\t".join(map(str, costs[u, u + 1:])) + "\n"
    return instance_str


def peace(input_path: Path, binary_path: Path, timeout: Optional[int] = None) -> Tuple[List[Dict], List[List[str]], float]:
    pace_format_instance = convert_pace_to_cm(input_path)

    with NamedTemporaryFile() as output_file, NamedTemporaryFile() as input_file:
        input_file.write(pace_format_instance.encode("utf8"))
        input_file.flush()

        out = run([binary_path, "--mode", "2", input_file.name, output_file.name], capture_output=True, timeout=timeout)

        if out.returncode != 0:
            raise RuntimeError(
                f"command failed with error code {out.returncode} and stderr '{out.stderr.decode('utf-8')}'.")

        stdout_output_str = out.stdout.decode("utf8")
        output_file_str = output_file.read().decode("utf8")

    statistics = []
    while True:
        idx = stdout_output_str.find("STATISTICS\n")
        if idx == -1:
            break
        stdout_output_str = stdout_output_str[idx + 11:]
        statistic_lines = stdout_output_str.split("\n")[:4]
        statistics.append(dict([tuple(line.split(": ")) for line in statistic_lines]))

    solution_lines = output_file_str.split("\n")[:-1]
    component_lines = solution_lines[:-1]
    components = [line.split(": ")[1].split(" ")[:-1] for line in component_lines]
    sum_of_costs = float(solution_lines[-1].split(": ")[1])

    return statistics, components, sum_of_costs


def parse_args() -> Tuple[Path, Path]:
    parser = argparse.ArgumentParser(
        description="Execute PEACE cluster editing.")

    parser.add_argument("--binary", type=str, default=None,
                        help="Path for PEACE binary.")

    parser.add_argument("--instance", type=str, required=True,
                        help="Path to instance in PACE format.")

    options = parser.parse_args()

    binary_path_str = options.binary if options.binary is not None else environ.get(BINARY_ENV_NAME)
    if binary_path_str is None:
        raise RuntimeError(
            f"binary must be specified as argument or environment variable {BINARY_ENV_NAME} must be specified.")
    binary_path = Path(binary_path_str)
    if not binary_path.exists():
        raise RuntimeError(f"binary {binary_path} does not exist.")

    instance_path = Path(options.instance)
    if not instance_path.exists():
        raise RuntimeError(f"instance {instance_path} does not exist.")

    return instance_path, binary_path


def main():
    instance_path, binary_path = parse_args()

    statistics, components, sum_of_costs = peace(instance_path, binary_path)

    print(f"time = {sum(float(s['time']) for s in statistics)}")
    print("components = ")
    for component in components:
        print(f"\t{' '.join(component)}")
    print(f"cost = {sum_of_costs}")
    print(f"{sum(float(s['costs']) for s in statistics)}")


if __name__ == '__main__':
    main()
