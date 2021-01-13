#!/bin/python3
import igraph
import numpy as np

import argparse
from glob import glob
from os import path
import os


def parse_graph(path):
    with open(path) as handle:
        for line in handle:
            if line.startswith("p cep"):
                split = line.strip().split(" ")
                n = int(split[2])
                m = int(split[3])
                break
        edges = []
        for line in handle:
            if line.startswith("c"):
                continue
            split = line.strip().split(" ")
            a = int(split[0]) - 1
            b = int(split[1]) - 1
            edges += [(a, b)]
        assert len(edges) == m
    return igraph.Graph(n, edges)

def graph_specs(graph: igraph.Graph):
    clustering = graph.transitivity_undirected()
#    diameter = graph.diameter()
    deg_dist_sd = graph.degree_distribution().sd
    deg_dist_var = graph.degree_distribution().var
    avg_deg = sum(graph.degree()) / len(graph.degree())
    max_deg = max(graph.degree())
#    max_cliques = graph.maximal_cliques()
#    num_max_cliques = len(max_cliques)
#    ω = max([len(c) for c in max_cliques])
    return clustering, deg_dist_var, deg_dist_sd, avg_deg, max_deg#, num_max_cliques, ω


def deg_dist(graph: igraph.Graph):
    degs = graph.degree()
    frequency = [0]*(max(degs)+1)
    for deg in degs:
        frequency[deg] += 1
    degree = list(range(len(frequency)))
    return frequency, degree


def write_csv(path, lines):
    with open(path, "w") as handle:
        csv = [line+"\n" for line in lines]
        handle.writelines(csv)


def build_fancy_csv_files(graph_dir, csv_dir):
    main_csv = []
    main_csv.append("graph,n,m,clustering,deg_dist_var,deg_dist_sd,avg_deg,max_deg")#,num_max_cliques,largest_clique")

    deg_csv = []
    deg_csv.append(f"graph,n,m,degree,frequency")

    if not path.exists(csv_dir):
        os.mkdir(csv_dir)

    for filename in glob(path.join(graph_dir, "*.gr")):
        graph = parse_graph(filename)
        specs = graph_specs(graph)
        specs = [str(item) for item in specs]

        name = path.basename(filename)
        print(f"Read graph {name}")
        n = graph.vcount()
        m = graph.ecount()
        line = f"{name},{n},{m},"
        line += ",".join(specs)
        main_csv.append(line)

        freqs, degs = deg_dist(graph)
        for freq, deg in zip(freqs, degs):
            if freq > 0:
                deg_csv.append(f"{name},{n},{m},{deg},{freq}")

    print("Writing Summary")
    write_csv(path.join(csv_dir, "summary.csv"), main_csv)
    write_csv(path.join(csv_dir, "degrees.csv"), deg_csv)
    print("Done.")


parser = argparse.ArgumentParser()
parser.add_argument("graph_dir")
parser.add_argument("csv_dir")
args = parser.parse_args()

build_fancy_csv_files(args.graph_dir, args.csv_dir)
