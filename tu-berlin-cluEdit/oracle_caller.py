import tempfile
from copy import deepcopy
from collections import deque

import pace2cm
import wrapper


def get_k(adj):
    temp = tempfile.NamedTemporaryFile(suffix=".cm")
    graph = (len(adj), -1, adj)
    pace2cm.write_cm(graph, temp.name)
    out = wrapper.call_tu_berlin_exec(temp.name)
    return wrapper.opt(out)


def edges(adj):
    edges = []
    for node, neighs in enumerate(adj):
        for neigh in neighs:
            if node < neigh:
                edges.append((node, neigh))
    return edges

def non_permanent_edge(adj, perm_edges):
    es = edges(adj)
    for edge in es:
        if edge not in perm_edges:
            return edge


def del_edge(adj, edge):
    adj2 = deepcopy(adj)
    a, b = edge
    adj2[a].remove(b)
    adj2[b].remove(a)
    return adj2


def reconstruct_solution(adj):
    opt = get_k(adj)
    print(f"First k: {opt}")

    permanent_edges = set()
    while len(edges(adj)) > len(permanent_edges):
        edge = non_permanent_edge(adj, permanent_edges)
        print(f"Trying without edge {edge}", end='')
        without_edge = del_edge(adj, edge)
        now_k = get_k(without_edge)
        if now_k < opt:
            adj = without_edge
            opt = now_k
            print(f" new k: {now_k} -> deleting edge")
        else:
            permanent_edges.add(edge)
            print(f" new k: {now_k} -> keeping edge")
    return adj


def partitions(adj):
    partitions = []
    component = [-1]*len(adj)

    comp = 0
    for node in range(len(adj)):
        if component[node] != -1:
            continue
        q = deque()
        q.append(node)
        while len(q) > 0:
            current = q.popleft()
            component[current] = comp
            for neigh in adj[current]:
                if component[neigh] == -1:
                    q.append(neigh)
        comp += 1
    return component


def clusters(components, plus=0):
    res = [list() for _ in range(max(components) + 1)]
    for node in range(len(components)):
        res[components[node]].append(node + plus)
    return res


def cost(adj, parti, clustering):
    adjacent = lambda a, b: b in adj[a]
    cost = 0

    for (a,b) in edges(adj):
        if parti[a] != parti[b]:
            cost += 1
    for cl in clustering:
        for a in cl:
            for b in cl:
                if a >= b:
                    continue
                if not adjacent(a, b):
                    cost += 1
    return cost




def dostuff(path):
    _,_,adj = pace2cm.read_pace(path)
    new_adj = reconstruct_solution(adj)
    print("Printing clusters:")
    parti = partitions(new_adj)
    clstrs = clusters(parti, 1)
    print(clstrs)
    print(cost(adj, parti, clstrs))

if __name__ == "__main__":
    import sys
    dostuff(sys.argv[1])
