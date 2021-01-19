#!/usr/bin/python3
from argparse import ArgumentParser

def read_pace(path):
    with open(path) as handle:
        lines = handle.readlines()
    n = m = -1
    for line in lines:
        if line.startswith("c"):
            continue
        line = line.strip()
        splt = line.split(" ")
        a = int(splt[-2])
        b = int(splt[-1])
        if line.startswith("p"):
            n = a
            m = b
            adj = [list() for _ in range(n)]
        else:
            a -= 1
            b -= 1
            adj[a].append(b)
            adj[b].append(a)
    return (n, m, adj)


def write_lines(lines, path):
    lines = [line + "\n" for line in lines]
    with open(path, 'w') as handle:
        handle.writelines(lines)


def write_cm(graph, path):
    n, m, adj = graph
    lines = []

    lines.append(str(n))
    for i in range(n):
        lines.append(str(i))

    for node in range(n-1):
        line = ''
        for other in range(node+1, n):
            if other in adj[node]:
                line += "1\t"
            else:
                line += "-1\t"
        lines.append(line)

    write_lines(lines, path)

def convert(in_path, cm_path):
    graph = read_pace(in_path)
    write_cm(graph, cm_path)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("in_path", help="PACE format input file")
    parser.add_argument("cm_path", help="Jena / PEACE format (.cm) output file")

    args = parser.parse_args()
    in_path = args.in_path
    cm_path = args.cm_path

    print(f"Reading graph from '{in_path}'")
    convert(in_path, cm_path)
    print(f"Done. Wrote to '{cm_path}'")

