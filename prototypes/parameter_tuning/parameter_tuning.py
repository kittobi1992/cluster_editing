import subprocess
from hashlib import sha256
from multiprocessing import Pool
from pathlib import Path
from typing import List, Union, Dict, Tuple, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from skopt import forest_minimize
from skopt.callbacks import DeadlineStopper
from skopt.plots import plot_objective
from skopt.space import Integer, Real
from skopt.utils import use_named_args

BINARY_PATH = Path("../../cmake-build-release/cluster_editing/application/exact_parameter_tuning")
assert BINARY_PATH.exists()


def execute_experiment(instance: int, params: Dict[str, Union[int, float]], time_limit: int) -> float:
    cmd = [BINARY_PATH, "--instance", str(instance), "--time_limit", str(time_limit)]
    for k, v in params.items():
        cmd += [f"--{k}", str(v)]
    out = subprocess.run(cmd, capture_output=True)
    assert out.returncode == 0
    solved, time = out.stdout.decode("utf-8").strip().split("\t")
    return float(time)


def build_callback(dimensions: List[Union[Integer, Real]], state_path: Path, figure_path: Path):
    def callback(res):
        df = pd.DataFrame(res.x_iters, columns=[d.name for d in dimensions])
        df["objective"] = res.func_vals
        for d in dimensions:
            df[d.name] = df[d.name].astype(d.dtype)

        df.to_csv(state_path, index=False)

        if len(res.models) >= 1:
            plot_objective(res, size=3)
            plt.tight_layout()
            plt.savefig(figure_path)

    return callback


def load_initial_points(dimensions: List[Union[Integer, Real]], state_path: Path) \
        -> Tuple[Optional[List[List[Union[int, float]]]], Optional[List[float]]]:
    x0 = None
    y0 = None
    if state_path.exists():
        df = pd.read_csv(state_path, index_col=False)
        x0 = [list(row) for row in np.array([df[d.name].values for d in dimensions]).T]
        y0 = list(df.loc[:, "objective"].values)
    return x0, y0


def main():
    df = pd.read_csv("../../results/2021-03-23-exact_with_star_bound.csv")
    df.sort_values(by="time", ascending=False, inplace=True)
    train_set = df.loc[(df["time"] > 5) & (df["time"] < 30), "instance"].values[:24]
    print(len(train_set), train_set)

    hex = sha256(",".join(map(str, train_set)).encode("utf8")).hexdigest()
    state_path = Path(f"state-{hex[:8]}.csv")
    figure_path = Path(f"objective-{hex[:8]}.png")
    checkpoint_path = Path(f"checkpoint-{hex[:8]}.pkl")

    dimensions = [
        Integer(name="max_num_unchanged", low=0, high=10),
        Real(name='probability', low=0.0, high=1.0)
    ]

    @use_named_args(dimensions)
    def objective(**params):
        time_limit = int(1.5 * 30)
        with Pool() as pool:
            times = pool.starmap(execute_experiment, [(i, params, time_limit) for i in train_set])
        times = np.array(list(times))
        return np.sqrt(np.sum(times ** 2))

    x0, y0 = load_initial_points(dimensions, state_path)

    callbacks = [
        build_callback(dimensions, state_path, figure_path),
        # CheckpointSaver(checkpoint_path, store_objective=False),
        DeadlineStopper(60 * 60)
    ]

    search_result = forest_minimize(
        objective, dimensions,
        x0=x0, y0=y0,
        n_calls=12, n_initial_points=10, initial_point_generator="random",
        callback=callbacks, verbose=True)

    print(search_result)


if __name__ == '__main__':
    main()
