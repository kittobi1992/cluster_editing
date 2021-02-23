#!/usr/bin/env python3
from pathlib import Path
from subprocess import run, TimeoutExpired
from peace import peace


def main():
    timeout = 60

    instances = list(Path("../instances/exact").glob("*.gr"))

    print("Graph,k,Total Time [s],Solved")
    for instance in sorted(instances, key=lambda x: x.name):
        try:
            statistics, components, sum_of_costs = peace(instance, Path("build/weighted_cluster_editing"), timeout)
            time = sum(float(s["time"]) for s in statistics)
            cost = int(sum_of_costs)
            solved = True
        except TimeoutExpired:
            time = timeout
            cost = ""
            solved = False
        print(f"{instance.name},{cost},{time},{solved}")


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        pass
