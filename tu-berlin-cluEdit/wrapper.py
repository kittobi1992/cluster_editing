#!/usr/bin/python3
# use python > 3.6
#
# Wrapper for java executable.
# Run function takes a pace graph and return the optimal value computed by the
# jar file.
#
from argparse import ArgumentParser
import subprocess
from os import path
import tempfile

import pace2cm

def call_tu_berlin_exec(graph, exe_dir='.', timeout=0):
    exe = path.join(exe_dir, "ClustEditOpt.jar")
    command = ["timeout", str(timeout), "java", "-jar", exe, graph, "NULL", "0", "0"]
    try:
        proc = subprocess.run(command, capture_output=True, text=True)
        lines = proc.stdout.splitlines()
        k = lines[-2].split(":")[-1]
        return float(k)
    except subprocess.CalledProcessError as e:
        if e.returncode == 124:
            print("Timeout")
        else:
            print(e)
        return -1
    except Exception as e:
        print(e)
        return -1

def run(pace_graph, **call_args):
    temp = tempfile.NamedTemporaryFile(suffix=".cm")
    pace2cm.convert(pace_graph, temp.name)
    return call_tu_berlin_exec(temp.name, **call_args)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("graph", help="Path to input graph in PACE format")
    parser.add_argument("--timeout", "-t", help="Timeout parameter (e.g. 10s, 5m, 24h)", default=None)
    parser.add_argument("--jar", "-j", help="directory of TU-Berlin jar executable", default=None)
    args = parser.parse_args()

    graph = args.graph
    timeout = args.timeout
    jar = args.jar
    kwargs = {}
    if timeout is not None:
        kwargs['timeout'] = timeout
    if jar is not None:
        kwargs['exe_dir'] = jar
    print(run(graph, **kwargs))
