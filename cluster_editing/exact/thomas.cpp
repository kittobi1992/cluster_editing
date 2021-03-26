
#include "thomas.h"

#include <vector>
#include <cassert>
#include <algorithm>

#include <cluster_editing/exact/solver.h>
#include <cluster_editing/exact/reductions.h>

using namespace std;


Instance createSubinst(Instance &graph, vector<int> &cluster) {
    auto inst = Instance(size(graph.edges));
    for(auto& row : inst.edges) fill(begin(row), end(row), 0);
    for(int i=0; i<size(inst.edges); ++i) inst.edges[i][i] = -1;
    inst.orig = inst.edges;
    for (auto node: cluster) {
        for (int i = 0; i < size(graph.edges); ++i) {
            inst.edges[i][node] = inst.edges[node][i] = graph.edges[i][node];
            inst.orig[i][node] = inst.orig[node][i] = graph.orig[i][node];
        }
    }
    for(int i=0; i<size(inst.edges); ++i)
        for(int j=0; j<size(inst.edges); ++j)
            assert(abs(inst.edges[i][j])==INF || inst.edges[i][j]==inst.orig[i][j]);
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
            if (in_cluster[i] && i!=node) {
                graph.edges[node][i] = graph.edges[i][node] = 1;
                graph.orig[node][i] = graph.orig[i][node] = 1;
            } else {
                graph.edges[node][i] = graph.edges[i][node] = -1;
                graph.orig[node][i] = graph.orig[i][node] = -1;
            }
        }
    }
}

int getSize(const Edges &g, const vector<int> &cluster) {
    int n = size(g);
    vector ininst(n, false);
    for (auto v : cluster) {
        ininst[v] = true;
        for (int i = 0; i < n; ++i)
            ininst[i] = ininst[i] | (g[v][i] > 0);
    }
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
        if (subinstSize == size(graph.edges) || subinstSize > 50) { // magic number
            //cout << "skipping cluster because it's too large" << endl;
            continue;
        }
        //cout << "Checking heuristic cluster of size " << size(cluster) << " and Subinstance size: " << subinstSize << endl;
        // put_together cost is an upper bound for subinst because it is a valid solution to make one cluster
        // put_together cost is smaller than heuristic cost because the heuristic solution has this cluster (thus paid put_together cost)
        auto subsolution = solve_exact(subinst, put_together_cost, 500); // TODO only do this until put_together_cost -1
        auto opt_cost = subsolution.cost;
        if(!subsolution.worked) continue;
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

std::optional<Instance> thomas_pairs(Instance inst) {
    // question: our thomas reduction does not only merge C but also excludes all other w from the cluster
    // the proof only says that C will end up in same cluster not that nothing else will
    // but in solving subinst we do not require C to be split
    // thus our implementation tests for a strictly stronger condition
    // namely that C will be a cluster on its own with nothing else in it

    // the following is strictly better than rule 3
    int n = size(inst.edges);
    auto& g = inst.edges;
    for(int u=0; u<n; u++) {
        for(int v=u+1; v<n; v++) {
            // try thomas reductions for cliques of size 2
            // for now we only calc what happens if u,v get separated
            if(g[u][v]<=0) continue;
            auto split_cost = g[u][v]; // lower bound on cost of putting u and v in different clusters
            auto keep_together_cost = 0; // lower bound on cost of keeping u,v in a cluster with possibly more nodes
            auto stand_alone_cost = 0; // cost of making {u,v} a cluster

            for(int w=0; w<n; ++w) {
                if(w==u || w==v) continue;
                if(g[u][w]<=0 && g[v][w]<=0) continue;
                keep_together_cost += g[u][w]>0&&g[v][w]>0 ? 0 : min(abs(g[u][w]), abs(g[v][w]));
                stand_alone_cost += max(0,g[u][w]) +max(0,g[v][w]);
                split_cost += min(max(0,g[u][w]), max(0,g[v][w]));
            }

            if(stand_alone_cost <= min(split_cost, keep_together_cost)) {
                auto res = merge(inst,u,v);
                for(int i=0; i<n-1;++i)
                    if(i!=u) {
                        res.spendCost += max(0,res.edges[u][i]);
                        res.edges[u][i] = res.edges[i][u] = -INF;
                    }
            }

            if(stand_alone_cost <= split_cost) {
                return merge(inst,u,v);
            }
        }
    }

    return {};
}


std::optional<Instance> heavy_edge_single_end(Instance inst) {

    int n = size(inst.edges);
    auto& g = inst.edges;
    for(int u=0; u<n; u++) {
        for(int v=0; v<n; v++) {
            if(v==u) continue;
            long long switch_cost = 0; // try to merge v into cluster of u
            for(int w=0; w<n; ++w) {
                if(w==u || w==v) continue;
                if(g[u][w]==-INF) switch_cost += max(0, g[v][w]);
                else switch_cost += abs(g[v][w]);
            }

            if(switch_cost <= g[u][v]) return merge(inst, u,v);
        }
    }

    return {};
}

std::optional<Instance> heavy_non_edge_single_end(Instance inst) {

    auto res = inst;
    int n = size(inst.edges);
    auto& g = res.edges;
    for(int u=0; u<n; u++) {
        for(int v=0; v<n; v++) {
            if(u==v) continue;
            if(g[u][v]>=0) continue;
            if(g[u][v]==-INF) continue;

            long long pull_out_cost = 0; // cost of pulling v out of u's cluster
            for(int w=0; w<n; ++w) {
                if(w==u || w==v) continue;
                pull_out_cost += max(0, g[v][w]);
            }

            if(pull_out_cost <= abs(g[u][v]))
                res.edges[u][v] = res.edges[v][u] = -INF;
        }
    }

    if(res.edges != inst.edges) return res;

    return {};
}
