#!/usr/bin/python3
# use python > 3.6
#
# Basic script to save solution sizes
#
import argparse
import glob
from os import path

import wrapper


def calc_scores(directory, timeout, verbose=True):
    scores = []
    for graph in glob.glob(path.join(directory, "*.gr")):
        k = wrapper.run(graph, timeout=timeout)
        entry = {
            "graph": path.basename(graph),
            "solution": k,
            "timeout": timeout
        }
        scores.append(entry)
        if verbose:
            print(entry)
    if verbose:
        print()
    return scores


def write_to_csv(scores, path):
    lines = []
    lines.append("graph,solution,timeout")
    for entry in scores:
        lines.append(f"{entry['graph']},{entry['solution']},{entry['timeout']}")

    lines = [l + "\n" for l in lines]
    with open(path, "w") as handle:
        handle.writelines(lines)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("graph_dir", help="directory with graphs in PACE format (*.gr)")
    parser.add_argument("timeout", help="timeout", default="10m")
    parser.add_argument("out_path", help="path to write stuff to")
    args = parser.parse_args()

    scores = calc_scores(args.graph_dir, args.timeout)
    write_to_csv(scores, args.out_path)

