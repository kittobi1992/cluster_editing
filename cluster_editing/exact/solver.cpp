
#include "solver.h"

#include <iostream>
#include <cassert>
#include <algorithm>
#include <chrono>

#include <cluster_editing/exact/reductions.h>
#include <cluster_editing/exact/lower_bounds.h>
#include <cluster_editing/exact/thomas.h>

#include "cluster_editing/multilevel.h"
#include "cluster_editing/io/graph_io.h"
#include "cluster_editing/metrics.h"
#include "cluster_editing/datastructures/graph_factory.h"


using namespace std;

// check if graph is union of cliques
std::optional<vector<vector<int>>> cliques(const Edges& graph) {
    auto compId = connectedComponents(graph);
    int n = size(graph);
    for(int v=0; v<n; ++v)
        for(int u=v+1; u<n; ++u)
            if((compId[v]==compId[u]) != (graph[v][u]>0))
                return {};
    int num = 1+ *max_element(begin(compId), end(compId));
    vector<vector<int>> res(num);
    for(int i=0; i<n; ++i)
        res[compId[i]].push_back(i);
    return res;
}

auto selectBranchingEdge(const Instance &graph) {
    int n = size(graph.edges);
    vector tripePerEdge(n, vector(n,0));
    for(int u=0; u<n; ++u) {
        for(int v=0; v<n; ++v) {
            auto& uv = graph.edges[u][v];
            for(int w=v+1; w<n && uv>0; ++w) {
                auto& uw = graph.edges[u][w];
                auto& vw = graph.edges[v][w];
                if(uw<=0 || vw>=0) continue;
                for(auto [a,b] : {pair(u,v), pair(u,w), pair(v,w)}) {
                    tripePerEdge[a][b]++;
                    tripePerEdge[b][a]++;
                }
            }
        }
    }


    auto u=0,v=1;
    auto value = [&](int a,int b) { return pair(graph.edges[a][b], tripePerEdge[a][b]); };
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (value(i,j) > value(u,v))
                u = i, v = j;
    assert(u != v);
    //assert(graph.edges[u][v] > 0);
    return pair(u,v);
}

Solution ExactSolver::solve_internal(Instance graph, int budget) {

    if(chrono::steady_clock::now()>time_limit) return {};

    // apply reductions
    int num_reduces = 0;
    bool changed = false;
    for (bool repeat = true; repeat; changed |= repeat) {

        //auto lower_bound = packing_lower_bound(graph.edges, budget-graph.spendCost);
        auto lower_bound = packing_local_search_bound(graph, budget-graph.spendCost);
        if(lower_bound + graph.spendCost > budget) {
            if (changed) {
                sumReductions += num_reduces;
                numReducingNodes++;
            }
            numPrunes++;
            return {};
        }

        if (auto opt=cliques(graph.edges); opt) { // check if graph is cliques
            if (changed) {
                sumReductions += num_reduces;
                numReducingNodes++;
            }
            Solution solution;
            solution.cost = graph.spendCost;
            solution.worked = true;
            for(auto& cl : *opt) {
                vector<int> expandedNodeSet;
                for(auto v : cl)
                    expandedNodeSet.insert(end(expandedNodeSet), begin(graph.idmap[v]), end(graph.idmap[v]));
                solution.cliques.push_back(expandedNodeSet);
            }
            return solution;
        }

        repeat = false;
        if(auto opt = forcedChoices(graph, budget); opt) graph = *opt, repeat=true;
        //if(auto opt = icxReductions(graph, budget); opt) graph = *opt, repeat=true;
        // more reductions
        num_reduces++;
    }

    if (changed) {
        sumReductions += num_reduces;
        numReducingNodes++;
    }

    //auto comps = connectedComponents(graph.edges);
    //auto num_comps = 1 + *max_element(begin(comps), end(comps));
    //if(num_comps>1) return solve_unconnected(graph, budget);

    branchingNodes++;

    auto [u,v] = selectBranchingEdge(graph);
    auto minstance = merge(graph, u, v); // merged instance
    auto finstance = graph; // forbidden instance
    finstance.spendCost += max(0, graph.edges[u][v]); // cost for deletion
    finstance.edges[u][v] = -INF;
    finstance.edges[v][u] = -INF;

    // try merging
    auto mergedSolution = solve_internal(minstance, budget);
    if (mergedSolution.worked) return mergedSolution;

    // try permanent deletion
    return solve_internal(finstance, budget);
}

Solution ExactSolver::solve(Instance inst, int budget_limit) {

    if(verbose) cout << "solving instance with n=" << size(inst.edges) << endl;
    // try some reductions for unweighted instances only
    auto isUnweighted = true;
    for(auto& row : inst.edges)
        for(auto val : row)
            isUnweighted &= val==1 || val==-1;
    if(isUnweighted) {
        if(auto opt = thomas(inst); opt) inst = *opt; // TODO multiple thomas reductions
        if(auto opt = forcedChoices(inst, solve_heuristic(inst).cost); opt) inst = *opt;
        if(auto opt = distance4Reduction(inst); opt) inst = *opt;
    }

    Solution s_comb;
    s_comb.cost = inst.spendCost;

    if(verbose) cout << "cost of initial reductions " << s_comb.cost << endl;
    if(verbose) cout << "instance size after reductions n=" << size(inst.edges) << endl;
    auto comps = constructConnectedComponents(inst);
    sort(begin(comps), end(comps), [](auto& a, auto& b){ return size(a.edges)<size(b.edges); });
    if(verbose) cout << "split into " << size(comps) << " CCs" << endl;
    for (auto& comp : comps) {

        Solution s;
        auto lower = packing_local_search_bound(comp, INF);
        if(verbose) cout << "start solving CC of size " << size(comp.edges) << " first bound " << lower << endl;

        for (int budget = lower; !s.worked && s_comb.cost+budget<=budget_limit; ++budget) {
            if(verbose) cout << budget << "   \r" << flush;
            s = solve_internal(comp, budget);

            if(chrono::steady_clock::now()>time_limit) {
                if(verbose) cout << endl << "time limit exceeded" << endl;
                return {};
            }
        }
        if(verbose) cout << endl;

        if(!s.worked) return {};

        // add cliques of component to solution for complete graph
        s_comb.cost += s.cost;
        for (auto& clique : s.cliques)
            s_comb.cliques.push_back(clique);
    }
    s_comb.worked = true;

    return s_comb;
}

Solution solve_exact(Instance inst, int budget_limit, int time_limit) {
    ExactSolver solver;
    solver.time_limit = chrono::steady_clock::now() + chrono::milliseconds(time_limit);
    return solver.solve(inst, budget_limit);
}

Solution solve_heuristic(const Instance &inst) {

    // ugly transform to unweighted instance
    int n = inst.edges.size();
    vector<vector<unsigned int>> adj(n);
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            assert(abs(inst.edges[i][j]) == 1);
            if (inst.edges[i][j] != 1) continue;
            adj[i].push_back(j);
            adj[j].push_back(i);
        }
    }

    cluster_editing::Context context;
    context.coarsening.algorithm = cluster_editing::CoarseningAlgorithm::lp_coarsener;
    context.refinement.lp.maximum_lp_iterations = 5;
    cluster_editing::Graph graph = cluster_editing::ds::GraphFactory::construct(adj);
    context.general.verbose_output = false;

    cluster_editing::multilevel::solve(graph, context);

    // build the solution
    Solution solution;
    const size_t edge_insertions = cluster_editing::metrics::edge_insertions(graph);
    const size_t edge_deletions = cluster_editing::metrics::edge_deletions(graph);
    solution.cost = inst.spendCost + edge_deletions + edge_insertions;

    vector<int> clique(n);
    int maxclq = 0;
    for (auto node: graph.nodes()) {
        clique[node] = graph.clique(node);
        if (clique[node] > maxclq)
            maxclq = clique[node];
    }
    vector<vector<int>> clusters(maxclq + 1);
    for (int i = 0; i < n; ++i)
        clusters[clique[i]].push_back(i);
    solution.cliques = clusters;

    return solution;
}


ostream &operator<<(ostream &os, const ExactSolver &rhs) {
    os << "branching nodes: " << rhs.branchingNodes << endl;
    os << "reductions:      " << rhs.numReducingNodes << endl;
    os << "disconnects:     " << rhs.numDisconnects << endl;
    os << "prunes:          " << rhs.numPrunes << endl;
    os << "iters/reduction: " << rhs.sumReductions * 1.0 / rhs.numReducingNodes << endl;
    return os;
}