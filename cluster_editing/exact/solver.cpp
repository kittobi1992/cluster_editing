
#include "solver.h"

#include <iostream>
#include <cassert>

#include "reductions.h"

using namespace std;

// isClique
bool isClique(const Edges &graph) {
    int n = size(graph);
    for (int u = 0; u < n; ++u)
        for (int v = u + 1; v < n; ++v)
            if (graph[u][v] <= 0)
                return false;

    return true;
}



auto selectBranchingEdge(const Instance &graph) {
    // find conflict triple
    int n = size(graph.edges);
    for (int v = 0; v < n; ++v) {
        for (int u = 0; u < n; ++u) {
            if (graph.edges[u][v] <= 0) continue; // uv must be edge
            for (int w = 0; w < n; ++w) {
                if (w == u || w == v) continue;
                if (graph.edges[u][w] <= 0) continue; // uw must be edge
                if (graph.edges[v][w] > 0) continue; // vw must be non-edge
                return pair(u, v);
            }
        }
    }
    return pair(-1, -1);
}

struct RunStatistics {
    int branchingNodes = 0;
    int numReducingNodes = 0;
    int sumReductions = 0;
    int numDisconnects = 0;
    int numPrunes = 0; // unsolvable leaf

};

ostream &operator<<(ostream &os, const RunStatistics &rhs) {
    os << "branching nodes: " << rhs.branchingNodes << endl;
    os << "reductions:      " << rhs.numReducingNodes << endl;
    os << "disconnects:     " << rhs.numDisconnects << endl;
    os << "prunes:          " << rhs.numPrunes << endl;
    os << "iters/reduction: " << rhs.sumReductions * 1.0 / rhs.numReducingNodes << endl;
    return os;
}

RunStatistics stats;

Solution solve(Instance graph, int budget, bool highL = false) {

    // apply reductions
    int num_reduces = 0;
    bool changed = false;
    for (bool repeat = true; repeat; changed |= repeat) {

        if (budget < graph.spendCost) {
            if (changed) {
                stats.sumReductions += num_reduces;
                stats.numReducingNodes++;
            }
            stats.numPrunes++;
            return {};
        }

        if (isClique(graph.edges)) {
            if (changed) {
                stats.sumReductions += num_reduces;
                stats.numReducingNodes++;
            }
            Solution solution;
            solution.cost = graph.spendCost;
            solution.worked = true;
            vector<int> expandedNodeSet;
            for (int v = 0; v < size(graph.edges); ++v)
                expandedNodeSet.insert(end(expandedNodeSet), begin(graph.idmap[v]), end(graph.idmap[v]));
            solution.cliques.push_back(expandedNodeSet);
            //cout << "found clique of size " << size(graph.edges) << endl;
            return solution;
        }

        repeat = false;
        auto r1 = icxReductions(graph, budget);
        if (r1) repeat = true, graph = *r1;
        // more reductions
        num_reduces++;
    }

    if (changed) {
        stats.sumReductions += num_reduces;
        stats.numReducingNodes++;
        return solveMaybeUnconnected(graph, budget);
    }

    //if(auto r = icxReductions(graph, budget); r)
    //return solveMaybeUnconnected(*r, budget);


    stats.branchingNodes++;
    // here we choose the edge to branch on
    auto[u, v] = selectBranchingEdge(graph);
    auto[icf, icp] = computeICFandICP(graph.edges);
    int bestReduction = min(icp[u][v], icf[u][v]);
    int n = size(graph.edges);
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (graph.edges[i][j] > 0 && min(icf[i][j], icp[i][j]) > bestReduction) {
                bestReduction = min(icf[i][j], icp[i][j]);
                u = i;
                v = j;
            }
        }
    }
    assert(u != v && graph.edges[u][v] > 0);


    auto minstance = merge(graph, u, v); // merged instance
    auto finstance = graph; // forbidden instance
    finstance.spendCost += max(0, graph.edges[u][v]); // cost for deletion
    finstance.edges[u][v] = -INF;
    finstance.edges[v][u] = -INF;

    for (int k = graph.spendCost; k <= budget; ++k) {

        if (highL) {
            cout << stats << endl;
            stats = RunStatistics{};
            cout << "===== " << k << "===== " << endl; // debug stuff
        }

        // try merging
        auto mergedSolution = solveMaybeUnconnected(minstance, k);
        assert(!mergedSolution.worked || mergedSolution.cost == k);
        if (mergedSolution.worked)
            return mergedSolution;

        // try permanent deletion
        auto forbiddenSolution = solveMaybeUnconnected(finstance, k);
        assert(!forbiddenSolution.worked || forbiddenSolution.cost == k);
        if (forbiddenSolution.worked)
            return forbiddenSolution;
    }

    return {};
}

Solution solveMaybeUnconnected(Instance graph, int budget, bool highL) {
    if (budget < graph.spendCost) return {}; // that should probably never happen
    Solution solution;
    solution.worked = true;
    solution.cost = graph.spendCost;
    auto comps = constructConnectedComponents(graph);
    if (size(comps) > 1) stats.numDisconnects++;
    for (auto comp : comps) {
        auto subsolution = solve(comp, budget - solution.cost, highL);
        if (!subsolution.worked || solution.cost + subsolution.cost > budget) {
            solution.worked = false;
            break;
        };
        solution.cost += subsolution.cost;
        // add cliques of component to solution for complete graph
        for (auto clique : subsolution.cliques)
            solution.cliques.push_back(clique);
    }
    return solution;
}
