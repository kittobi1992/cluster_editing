import argparse
import subprocess
from collections import namedtuple
from hashlib import sha256
from multiprocessing import Pool
from pathlib import Path
from typing import List, Union, Dict, Tuple, Optional, Callable, Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from skopt import forest_minimize
from skopt.callbacks import DeadlineStopper
from skopt.plots import plot_objective
from skopt.space import Integer, Real
from skopt.utils import use_named_args

Config = namedtuple("Config",
                    ["instances", "time_limit", "max_num_unchanged", "n_calls", "n_initial_points", "deadline", "loss"])


def execute_experiment(binary_path: Path, instance: int, params: Dict[str, Union[int, float]],
                       time_limit: int) -> float:
    cmd = [binary_path, "--instance", str(instance), "--time_limit", str(time_limit)]
    for k, v in params.items():
        cmd += [f"--{k}", str(v)]
    out = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if out.returncode != 0:
        raise RuntimeError(out.stderr.decode("utf-8"))
    solved, time = out.stdout.decode("utf-8").strip().split("\t")
    time = float(time)
    if solved != "true" or time >= time_limit:
        time = time_limit
    return time


def build_callback(dimensions: List[Union[Integer, Real]], config: Config, state_path: Path,
                   figure_path: Optional[Path]) -> Callable[[Any], None]:
    def callback(res):
        df = pd.DataFrame(res.x_iters, columns=[d.name for d in dimensions])
        df["objective"] = res.func_vals
        for d in dimensions:
            df[d.name] = df[d.name].astype(d.dtype)

        df["config"] = str(config)

        df.to_csv(state_path, index=False)

        if len(res.models) >= 1 and figure_path is not None:
            plot_objective(res, size=3)
            plt.tight_layout()
            plt.savefig(figure_path)

    return callback


def build_objective(dimensions: List[Union[Integer, Real]], config: Config, binary_path: Path) \
        -> Callable[[Any], float]:
    assert binary_path.exists()

    loss_functions = dict(rmse=lambda xs: np.sqrt(np.mean(xs ** 2)), softmax=lambda xs: np.log(np.sum(np.exp(xs))))
    loss_function = loss_functions.get(config.loss)

    @use_named_args(dimensions)
    def objective(**params) -> float:
        with Pool() as pool:
            times = pool.starmap(execute_experiment,
                                 [(binary_path, i, params, config.time_limit) for i in config.instances])
        times = np.array(list(times))
        value = loss_function(times)
        print(times)
        print(value)
        return value

    return objective


def load_initial_points(dimensions: List[Union[Integer, Real]], state_path: Path) \
        -> Tuple[Optional[List[List[Union[int, float]]]], Optional[List[float]]]:
    x0 = None
    y0 = None
    if state_path.exists():
        df = pd.read_csv(state_path, index_col=False)
        x0 = [list(row) for row in np.array([df[d.name].values for d in dimensions]).T]
        y0 = list(df.loc[:, "objective"].values)
    return x0, y0


def parse_args(previous_results_path: Path) -> Config:
    train_set_lb = 1 * 60
    train_set_ub = 40 * 60

    parser = argparse.ArgumentParser()

    instance_group = parser.add_mutually_exclusive_group()
    instance_group.add_argument("--instances", nargs="+", type=int)
    instance_group.add_argument("--instances_from_previous_time", nargs=2, type=float)

    parser.add_argument("--time_limit", type=int, default=60 * 60)
    parser.add_argument("--max_num_unchanged", nargs=2, type=int, default=(0, 20))
    parser.add_argument("--n_calls", type=int, default=100)
    parser.add_argument("--n_initial_points", type=int, default=100)
    parser.add_argument("--deadline", type=int, default=2 * 24 * 60 * 60)
    parser.add_argument("--loss", choices=("rmse", "softmax"), default="rmse")

    options = parser.parse_args()

    instances = options.instances
    if instances is None:
        if options.instances_from_previous_time is not None:
            train_set_lb = options.instances_from_previous_time[0]
            train_set_ub = options.instances_from_previous_time[1]
        assert (previous_results_path.exists())
        df = pd.read_csv(previous_results_path)
        df.sort_values(by="time", ascending=False, inplace=True)
        instances = list(df.loc[(df["time"] > train_set_lb) & (df["time"] < train_set_ub), "instance"].values)

    max_num_unchanged = tuple(options.max_num_unchanged)

    config = Config(instances, options.time_limit, max_num_unchanged, options.n_calls, options.n_initial_points,
                    options.deadline, options.loss)
    return config


def main():
    previous_results_path = Path("../../results/2021-03-23-exact_with_star_bound.csv")

    config = parse_args(previous_results_path)

    print(len(config.instances), config.instances)

    hex_code = sha256(str(config).encode("utf8")).hexdigest()[:12]
    state_path = Path(f"state-{hex_code}.csv")
    figure_path = Path(f"objective-{hex_code}.png")
    binary_path = Path("../../cmake-build-release/cluster_editing/application/exact_parameter_tuning")
    # checkpoint_path = Path(f"checkpoint-{hex_code}.pkl")

    dimensions = [
        Integer(name="max_num_unchanged", low=config.max_num_unchanged[0], high=config.max_num_unchanged[1]),
        Real(name='probability', low=0.0, high=1.0)
    ]

    objective = build_objective(dimensions, config, binary_path)

    x0, y0 = load_initial_points(dimensions, state_path)

    callbacks = [
        build_callback(dimensions, config, state_path, figure_path),
        # CheckpointSaver(checkpoint_path, store_objective=False),
        DeadlineStopper(config.deadline)
    ]

    search_result = forest_minimize(
        objective, dimensions,
        x0=x0, y0=y0,
        n_calls=config.n_calls, n_initial_points=config.n_initial_points, initial_point_generator="sobol",
        callback=callbacks, verbose=True)

    print(search_result)


if __name__ == '__main__':
    main()
