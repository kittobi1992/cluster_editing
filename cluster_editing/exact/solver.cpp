/*******************************************************************************
 * This file is part of KaPoCE.
 *
 * KaPoCE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaPoCE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaPoCE.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "solver.h"

#include <iostream>
#include <cassert>
#include <algorithm>
#include <chrono>

#include <cluster_editing/exact/reductions.h>
#include <cluster_editing/exact/lower_bounds.h>
#include <cluster_editing/exact/thomas.h>
#include <cluster_editing/exact/star_bound.h>
#include <cluster_editing/heuristic/clustering.h>

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

    // apply reductions
    int num_reduces = 0;
    bool changed = false;
    for (bool repeat = true; repeat; changed |= repeat) {

        if(chrono::steady_clock::now()>time_limit) return {};

        //auto lower_bound = packing_lower_bound(graph.edges, budget-graph.spendCost);
        auto lower_bound = meta_lower_bound(graph, budget-graph.spendCost);
        if(lower_bound + graph.spendCost > budget) {
            numReducingNodes += changed;
            numPrunes++;
            return {};
        }

        if (auto opt=cliques(graph.edges); opt) { // check if graph is cliques
            numReducingNodes += changed;

            Solution solution;
            solution.cost = graph.spendCost;
            solution.worked = true;
            solution.cliques = graph.done_clusters;
            for(auto& cl : *opt) {
                vector<int> expandedNodeSet;
                for(auto v : cl)
                    expandedNodeSet.insert(end(expandedNodeSet), begin(graph.idmap[v]), end(graph.idmap[v]));
                solution.cliques.push_back(expandedNodeSet);
            }
            return solution;
        }

        repeat = false;
        if(auto opt = forcedChoicesStarBound(graph, budget, false); opt) graph = *opt, repeat=true, redForcedStar++;
        if(auto opt = forcedChoices(graph, budget); opt) graph = *opt, repeat=true, redForced++;
        if(auto opt = simpleTwin(graph); opt) graph = *opt, repeat=true, redTwin++;
        if(auto opt = complexTwin(graph,true); opt) graph = *opt, repeat=true, redTwin2++;
        if(auto opt = icxReductions(graph, budget); opt) graph = *opt, repeat=true, redICX++;
        if(auto opt = thomas_pairs(graph); opt) graph = *opt, repeat = true, redThomas2++;
        if(auto opt = heavy_edge_single_end(graph); opt) graph = *opt, repeat = true, redHeavyEdge++;
        if(auto opt = heavy_non_edge_single_end(graph); opt) graph = *opt, repeat = true, redHeavyNonEdge++;
        // more reductions
        num_reduces++;
    }

    if (changed) numReducingNodes++;
    branchingNodes++;

    auto comps = connectedComponents(graph.edges);
    auto num_comps = 1 + *max_element(begin(comps), end(comps));
    if(num_comps>1) numDisconnects++;

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
        auto t1 = chrono::steady_clock::now();
        auto upper = solve_heuristic(inst).cost;
        if(verbose) cout << "upper bound " << upper << endl;
        if(auto opt = thomas(inst); opt) inst = *opt; // TODO multiple thomas reductions
        if(auto opt = distance4Reduction(inst); opt) inst = *opt;
        while(true) {
            string applied = "";
            if(empty(applied)) if(auto opt = force_small_components(inst); opt) inst = *opt, applied = "small clean up";
            if(empty(applied)) if(auto opt = forcedChoices(inst, upper); opt) inst = *opt, applied = "force p3";
            if(empty(applied)) if(auto opt = forcedChoicesStarBound(inst, upper, false); opt) inst = *opt, applied = "force star";
            if(empty(applied)) if(auto opt = simpleTwin(inst); opt) inst = *opt, applied = "twin simple";
            if(empty(applied)) if(auto opt = complexTwin(inst,true); opt) inst = *opt, applied = "twin complex";
            if(empty(applied)) if(auto opt = icxReductions(inst, upper); opt) inst = *opt, applied = "icx";
            if(empty(applied)) if(auto opt = thomas_pairs(inst); opt) inst = *opt, applied = "heavy edge (b)";
            if(empty(applied)) if(auto opt = heavy_edge_single_end(inst); opt) inst = *opt, applied = "heavy edge (s)";
            if(empty(applied)) if(auto opt = heavy_non_edge_single_end(inst); opt) inst = *opt, applied = "heavy non-edge";
            if(empty(applied)) if(auto opt = forcedChoicesSingleMerge(inst, upper, false); opt) inst = *opt, applied = "forced single merge";
            if(empty(applied)) break;
            else if(verbose) cout << "reduced to n=" << size(inst.edges)
                << " -INFs=" << forbiddenEdges(inst)
                << " spent=" << inst.spendCost
                << " with " << applied << endl;
        }
        auto t2 = chrono::steady_clock::now();
        int lower = star_bound(inst,upper);
        if(verbose) {
            cout << "INITIAL REDUCTION FINISHED" << endl;
            cout << "time:  " << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() * 0.001 << " s\n";
            cout << "size:  " << size(inst.edges) << endl;
            cout << "lower: " << inst.spendCost + lower << endl;
            cout << "upper: " << upper << endl;
            cout << "gap:   " << upper-lower-inst.spendCost << endl;
            cout << "STARTING BRANCH AND BOUND" << endl;
        }
        root_size = size(inst.edges);
        root_time = chrono::duration_cast<chrono::milliseconds>(t2-t1).count();
        root_gap = upper-lower-inst.spendCost;
    }

    Solution s_comb;
    s_comb.cost = inst.spendCost;
    s_comb.cliques = inst.done_clusters;

    auto comps = constructConnectedComponents(inst);
    sort(begin(comps), end(comps), [](auto& a, auto& b){ return size(a.edges)<size(b.edges); });
    if(verbose) cout << "split into " << size(comps) << " CCs" << endl;
    for (auto& comp : comps) {

        Solution s;
        auto lower = meta_lower_bound(comp, INF);
        if(verbose) cout << "start solving CC of size " << size(comp.edges) << " first bound " << lower << endl;

        for (int budget = lower; !s.worked && s_comb.cost+budget<=budget_limit; ++budget) {
            if(verbose) cout << budget << "\t (total: " << s_comb.cost + budget << ")\r" << flush;
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

bool cluster_editing::utils::CommonOperations::dirty = false;
Solution solve_heuristic(const Instance &inst) {

    // ugly transform to unweighted instance
    int n = inst.edges.size();
    if(auto opt=cliques(inst.edges); opt) {
        Solution s;
        s.worked = true;
        s.cost = inst.spendCost;
        s.cliques = *opt; // TODO this (and the code below) does not yield expanded clusters but local indexed vertices instead
        return s;
    }

    vector<vector<unsigned int>> adj(n);
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            assert(abs(inst.edges[i][j]) == 1);
            if (inst.edges[i][j] != 1) continue;
            adj[i].push_back(j);
            adj[j].push_back(i);
        }
    }
    cluster_editing::Graph graph = cluster_editing::ds::GraphFactory::construct(adj);

    // TODO bring over evo
    // TODO stop if lower bound reached
    cluster_editing::utils::CommonOperations::dirty = true;

    cluster_editing::Context context;
    context.refinement.use_lp_refiner = true;
    context.refinement.lp.maximum_lp_iterations = 500;
    context.refinement.lp.random_shuffle_each_round = true;
    context.refinement.lp.node_order = cluster_editing::NodeOrdering::random_shuffle;
    context.refinement.lp.min_improvement = 5;
    context.refinement.lp.early_exit_window = 100;
    context.general.time_limit = 60.0 * 30.0;
    context.general.verbose_output = false;

    size_t num_reps = 200;
    size_t best_num_edits = graph.numEdges();
    vector<int> clique_assignment(n);
    auto update = [&] {
        size_t num_edits = cluster_editing::metrics::edits(graph);
        if (num_edits < best_num_edits) {
            best_num_edits = num_edits;
            for (auto u : graph.nodes()) {
                clique_assignment[u] = graph.clique(u);
            }
        }
        graph.reset();
    };

    for (size_t i = 0; i < num_reps; ++i) {
        cluster_editing::solve(graph, context);
        update();
    }

    // build the solution
    Solution solution;
    solution.cost = inst.spendCost + best_num_edits;

    int maxclq = *std::max_element(clique_assignment.begin(), clique_assignment.end());
    vector<vector<int>> clusters(maxclq + 1);
    for (int i = 0; i < n; ++i)
        clusters[clique_assignment[i]].push_back(i);
    solution.cliques = clusters;
    return solution;
}


ostream &operator<<(ostream &os, const ExactSolver &rhs) {
    os << "branching nodes: " << rhs.branchingNodes << endl;
    os << "reductions:      " << rhs.numReducingNodes << endl;
    os << "\tforced (star): " << rhs.redForcedStar << endl;
    os << "\tforced (p3):   " << rhs.redForced << endl;
    os << "\ttwin simple:   " << rhs.redTwin << endl;
    os << "\ttwin complex:  " << rhs.redTwin2 << endl;
    os << "\ticx:           " << rhs.redICX << endl;
    os << "\th edge (b):    " << rhs.redThomas2 << endl;
    os << "\th edge (s):    " << rhs.redHeavyEdge << endl;
    os << "\th non edge:    " << rhs.redHeavyNonEdge << endl;
    os << "disconnects:     " << rhs.numDisconnects << endl;
    os << "prunes:          " << rhs.numPrunes << endl;
    //os << "iters/reduction: " << rhs.sumReductions * 1.0 / rhs.numReducingNodes << endl;
    return os;
}
