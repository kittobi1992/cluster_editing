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

def call_tu_berlin_exec(graph, exe=None, timeout=0):
    if exe is None:
        exe = "./new_jar.jar"
    command = ["timeout", str(timeout), "java", "-jar", exe, graph, "NULL", "0", "0"]
    try:
        proc = subprocess.run(command, capture_output=True, text=True)
        return proc.stdout
    except subprocess.CalledProcessError as e:
        if e.returncode == 124:
            print("Timeout")
        else:
            print(e)
        return ""
    except Exception as e:
        print(e)
        return ""

def extract_key(stdout, key):
    try:
        for line in reversed(stdout.splitlines()):
            if key in line:
                val = line.strip().split(":")[-1]
                return float(val.strip())
        return -1
    except:
        return -1

def lowerbound(stdout):
    return extract_key(stdout, "lowerbound:")

def upperbound(stdout):
    return extract_key(stdout, "upperbound:")

def opt(stdout):
    return extract_key(stdout, "opt k:")

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
    out = run(graph, **kwargs)
    print(f"opt: {opt(out)}")
    print(f"lower: {lowerbound(out)}")
    print(f"upper: {upperbound(out)}")
