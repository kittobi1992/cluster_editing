#!/usr/bin/env python3
# coding: utf-8

import itertools
import random
import collections
import copy

def readGraph(path):
    problem_read = False
    n = 0
    m = 0
    result = None
    with open(path, 'r') as fin:
        for line in fin:
            if line[0] == 'c':
                continue
            if line[0] == 'p':
                assert not problem_read
                problem_read = True
                assert line[0:5] == 'p cep'
                parts = line.split(' ')
                assert len(parts) == 4
                n = int(parts[2])
                m = int(parts[3])
                result = {u: set() for u in range(n)}
            else:
                assert problem_read
                parts = line.split(' ')
                assert(len(parts) == 2)
                u, v = map(int, parts)
                u -= 1
                v -= 1
                result[u].add(v)
                result[v].add(u)
    assert sum(map(len, result.values())) == 2 * m
    return result


def subGraph(g, nodes):
    return {u: {v for v in g[u] if v in nodes} for u in nodes}


def degeneracyOrdering(g):
    result = []
    degree = {u: len(g[u]) for u in g}
    node_by_degree = []
    for u, d in degree.items():
        while len(node_by_degree) <= d:
            node_by_degree.append(set())
        node_by_degree[d].add(u)
    for d in range(len(node_by_degree)):
        # Convert to list as we can append additional nodes here
        current_nodes = list(node_by_degree[d])
        for u in current_nodes:
            for v in g[u]:
                if v in degree:
                    dv = degree[v]
                    if dv > d:
                        node_by_degree[dv].remove(v)
                        if dv > d + 1:
                            node_by_degree[dv - 1].add(v)
                        else:
                            current_nodes.append(v)
                    degree[v] -= 1
            del degree[u]
            result.append(u)
    return result


def coloring(g):
    color = dict()
    colored_nodes = collections.defaultdict(list)
    order = degeneracyOrdering(g)
    order.reverse()
    for u in order:
        neighbor_colors = {color[v] for v in g[u] if v in color}
        c = 0
        while c in neighbor_colors:
            c += 1
        color[u] = c
        colored_nodes[c].append(u)

    for u in g:
        for v in g[u]:
            assert color[u] != color[v]
    return colored_nodes


def independentSet(g):
    result = set()
    # Greedy degree heuristic
    for u in sorted(g.keys(), key=lambda x: len(g[x])):
        if not any(v in result for v in g[u]):
            result.add(u)
    # 10 rounds two-improvements and prefer lower-degree replacements when possible
    for i in range(10):
        result_list = list(result)
        random.shuffle(result_list)
        for u in result_list:
            result.remove(u)
            candidate_found = False
            single_candidates = []
            for v in g[u]:
                if not any((x in result for x in g[v])):
                    result.add(v)
                    single_candidates.append(v)
                    for w in g[u]:
                        if w not in result and not any((x in result for x in g[w])):
                            result.add(w)
                            candidate_found = True
                    if candidate_found:
                        break
                    else:
                        result.remove(v)
            if not candidate_found:
                if single_candidates:
                    result.add(min(single_candidates, key=lambda x: len(g[x])))
                else:
                    result.add(u)
    # Ensure that we actually calculated an independent set
    for u in result:
        assert(not any((v in result for v in g[u])))
    return result


def lowerBound(g):
    # All node pairs that are already covered by the bound
    bound_pairs = set()
    bound = 0
    stars = []
    # Start with highest-degree node
    for i, u in enumerate(sorted(g.keys(), reverse=False, key=lambda x: len(g[x]))):
        # Exclude nodes where (u, v) is in the bound
        candidates = [v for v in g[u] if (u, v) not in bound_pairs]
        neighborgraph = subGraph(g, candidates)
        # Add edges for all node pairs that are in the bound
        for x, y in itertools.combinations(candidates, 2):
            if (x, y) in bound_pairs:
                neighborgraph[x].add(y)
                neighborgraph[y].add(x)
        while True:
            # A star where all node pairs are not covered can be added to the bound
            star = independentSet(neighborgraph)
            if len(star) < 2:
                break
            bound += len(star) - 1
            for x in star:
                for y in neighborgraph[x]:
                    neighborgraph[y].remove(x)
                del neighborgraph[x]
            star = [u] + list(star)
            stars.append(star)
            for x, y in itertools.combinations(star, 2):
                bound_pairs.add((x, y))
                bound_pairs.add((y, x))
        #print(i, u, bound)
    #print(stars)
    return bound


class LowerBound:
    def __init__(self, g):
        self.g = g
        # All node pairs that are already covered by the bound
        self.used_pairs = set()
        self.bound = 0
        self.stars = [set() for u in g]
        self.free_edges_g = copy.deepcopy(g)

    def free_neighbors(self, u):
        return self.free_edges_g[u]

    def pair_used(self, u, v):
        return (u, v) in self.used_pairs

    def has_star(self, star):
        star = tuple(star)
        return star in self.stars[star[0]]

    def stars_in_random_order(self):
        all_stars = [star for u_stars in self.stars for star in u_stars]
        random.shuffle(all_stars)
        return all_stars

    def can_add(self, star):
        return not any(p in self.used_pairs
                       for p in itertools.combinations(star, 2))

    def add_star(self, star):
        assert len(star) > 2
        star = tuple(star)
        for i, u in enumerate(star[1:]):
            assert star[0] in self.g[u]
            for v in star[i+2:]:
                assert v not in self.g[u]
        self.stars[star[0]].add(star)
        self.bound += len(star) - 2
        for x, y in itertools.combinations(star, 2):
            for u, v in [(x, y), (y, x)]:
                assert (u, v) not in self.used_pairs
                self.used_pairs.add((u, v))
                if v in self.g[u]:
                    self.free_edges_g[u].remove(v)

    def remove_star(self, star):
        star = tuple(star)
        self.stars[star[0]].remove(star)
        self.bound -= len(star) - 2
        for x, y in itertools.combinations(star, 2):
            for u, v in [(x, y), (y, x)]:
                self.used_pairs.remove((u, v))
                if v in self.g[u]:
                    self.free_edges_g[u].add(v)

    def get_candidates(self, u, v):
        if v in self.g[u]:
            return [(v, u, y) for y in self.free_edges_g[v] if y != u and
                    y not in self.g[u] and (u, y) not in self.used_pairs] + \
                   [(u, v, y) for y in self.free_edges_g[u] if y != v and
                    y not in self.g[v] and (v, y) not in self.used_pairs]
        return [(y, v, u) for y in self.free_edges_g[u] & self.free_edges_g[v]]


def lowerBoundColoring(g):
    bound = LowerBound(g)
    # Start with highest-degree node
    for i, u in enumerate(sorted(g.keys(), reverse=True, key=lambda x: len(g[x]))):
        # Exclude nodes where (u, v) is in the bound
        candidates = bound.free_neighbors(u)
        neighborgraph = subGraph(g, candidates)
        # Add edges for all node pairs that are in the bound
        # TODO: make this faster if not many pairs are used
        for x, y in itertools.combinations(candidates, 2):
            if bound.pair_used(x, y):
                neighborgraph[x].add(y)
                neighborgraph[y].add(x)
        # A star where all node pairs are not covered can be added to the bound
        u_stars = coloring(neighborgraph)
        for star in u_stars.values():
            if len(star) < 2:
                continue
            star = [u] + list(star)
            bound.add_star(star)
        # print(i, u, bound)
    for i in range(5):
        print(bound.bound)
        for star in bound.stars_in_random_order():
            # Check if for some reason we have removed this star
            if not bound.has_star(star):
                continue

            candidate_star = list(star)
            # Remove nodes one by one instead of removing all
            for v in list(star[1:]):  # Force a copy
                # Try removing v from star and inserting at least two P3
                bound.remove_star(candidate_star)

                is_p3 = len(candidate_star) == 3

                if not is_p3:
                    assert len(candidate_star) > 3
                    candidate_star.remove(v)
                    bound.add_star(candidate_star)

                candidates_per_pair = dict()
                if is_p3:
                    for x, y in itertools.combinations(candidate_star, 2):
                        candidates_per_pair[(x, y)] = bound.get_candidates(x, y)
                else:
                    for x in candidate_star:
                        assert v not in g[x] or x == star[0]
                        candidates_per_pair[(v, x)] = bound.get_candidates(v, x)

                two_found = False
                for pair, candidates in candidates_per_pair.items():
                    assert not two_found
                    for candidate in candidates:
                        assert not two_found
                        bound.add_star(candidate)

                        for pair2, candidates2 in candidates_per_pair.items():
                            if pair == pair2:
                                continue
                            for candidate2 in candidates2:
                                if bound.can_add(candidate2):
                                    two_found = True
                                    bound.add_star(candidate2)
                        if two_found:
                            break
                        bound.remove_star(candidate)
                    if two_found:
                        break
                # Re-insert v
                if not two_found:
                    # Replace by random candidate
                    all_candidates = list(set((c for c in candidates for candidates in candidates_per_pair.values())))
                    if all_candidates:
                        replacement = random.choice(all_candidates)
                        bound.add_star(replacement)
                    else:
                        if not is_p3:
                            bound.remove_star(candidate_star)
                            candidate_star.append(v)
                        bound.add_star(candidate_star)
                if is_p3:
                    break
    return bound.bound


graph = readGraph("/home/michael/graphs/PACE2021/exact/exact191.gr")
lowerBoundColoring(graph)
