
#include "thomas.h"

#include <vector>
#include <cassert>
#include <algorithm>

#include <cluster_editing/exact/solver.h>

using namespace std;


Instance createSubinst(Instance &graph, vector<int> &cluster) {
    auto inst = Instance(size(graph.edges));
    for(auto& row : inst.edges) fill(begin(row), end(row), 0);
    for (auto node: cluster) {
        for (int i = 0; i < size(graph.edges); ++i) {
            inst.edges[node][i] = graph.edges[node][i];
            inst.edges[i][node] = graph.edges[i][node];
        }
    }
    return inst;
}

int calc_together_cost(vector<int> &cluster, Instance &graph) {
    auto in_cluster = vector<bool>(size(graph.edges));
    for (auto node: cluster)
        in_cluster[node] = true;
    long long cost = 0;

    // negative edges in cluster
    for (auto node: cluster)
        for (auto other: cluster)
            if (node < other && graph.edges[node][other] < 0)
                cost -= graph.edges[node][other];

    // positive edges to the outside
    for (auto cluster_node: cluster)
        for (int i = 0; i < size(graph.edges); ++i)
            if (!in_cluster[i] && graph.edges[cluster_node][i] > 0)
                cost += graph.edges[cluster_node][i];

    return max<int>(-INF, cost);
}

// makes cluster to a clique in instance
void put_together(Instance &graph, vector<int> &cluster) {

    auto in_cluster = vector<bool>(size(graph.edges));
    for (auto node: cluster)
        in_cluster[node] = true;

    for (auto node: cluster) {
        for (int i = 0; i < size(graph.edges); ++i) {
            if (in_cluster[i]) {
                graph.edges[node][i] = 1;
                graph.edges[i][node] = 1;
            } else {
                graph.edges[node][i] = -1;
                graph.edges[i][node] = -1;
            }
        }
    }
}

int getSize(const Edges &g, const vector<int> &cluster) {
    int n = size(g);
    vector ininst(n, false);
    for (auto v : cluster)
        for (int i = 0; i < n; ++i)
            ininst[i] = ininst[i] | g[v][i] > 0;
    int res = 0;
    for (auto bit : ininst) res += bit;
    return res;
}

optional<Instance> thomas(Instance graph) {
    auto solution = solve_heuristic(graph);

    bool changed = false;

    sort(begin(solution.cliques), end(solution.cliques), [](auto &a, auto &b) { return size(a) < size(b); });
    for (auto cluster: solution.cliques) {
        if (size(cluster) < 2)
            continue;
        Instance subinst = createSubinst(graph, cluster);
        int put_together_cost = calc_together_cost(cluster, graph);
        if (put_together_cost == 0) {
            continue;
        }

        auto subinstSize = getSize(subinst.edges, cluster);
        if (subinstSize == size(graph.edges) || subinstSize > 25) { // magic number
            //cout << "skipping cluster because it's too large" << endl;
            continue;
        }
        //cout << "Checking heuristic cluster of size " << size(cluster) << " and Subinstance size: " << subinstSize << endl;
        // put_together cost is an upper bound for subinst because it is a valid solution to make one cluster
        // put_together cost is smaller than heuristic cost because the heuristic solution has this cluster (thus paid put_together cost)
        auto subsolution = solve_exact(subinst, put_together_cost); // TODO only do this until put_together_cost -1
        auto opt_cost = subsolution.cost;
        assert(subsolution.worked);
        if (put_together_cost == opt_cost) {
            //cout << "Found safe cluster of size " << size(cluster) << endl;
            //cout << size(cluster) << ' ' << flush;
            put_together(graph, cluster);
            graph.spendCost += put_together_cost;
            changed = true;
        } else {
            //cout << "reduction was not applicable" << endl;
        }
    }

    if(!changed) return {};
    return graph;
}
