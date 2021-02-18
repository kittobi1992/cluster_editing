
#include "solver.h"

#include <iostream>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <chrono>

#include <cluster_editing/exact/reductions.h>
#include <cluster_editing/exact/thomas.h>

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

auto packing_lower_bound(Edges edges, int budget) {
    int n = size(edges);

    long long best = 0;
    for(int u=0; u<n; ++u) {
        for(int v=u+1; v<n; ++v) {
            auto& uv = edges[u][v];
            if(uv<=0) continue;
            for(int w=v+1; w<n && uv>0; ++w) {
                auto& uw = edges[u][w];
                auto& vw = edges[v][w];
                if(uw<=0 || vw>=0) continue;
                auto mn = min(uv,min(uw,-vw));
                best += mn;
                if(best>budget) return best;
                uv -= mn;
                uw -= mn;
                vw += mn;
            }
        }
    }

    return best;
}

auto packing_local_search_bound(Edges edges, int limit) {
    using Triple = array<int,4>; // array (u,v,w,cost) represents triple v-u-w with cost being min(uv,uw,-vw)

    int n = size(edges);
    long long cost = 0; // maintains sum of packing[i].cost
    auto apply = [&](const Triple& t, bool undo=false) {
        auto [u,v,w,uvw_cost] = t;
        if(undo) uvw_cost = -uvw_cost;
        cost += uvw_cost;
        edges[v][u] = edges[u][v] = edges[u][v] - uvw_cost;
        edges[w][u] = edges[u][w] = edges[u][w] - uvw_cost;
        edges[w][v] = edges[v][w] = edges[v][w] + uvw_cost;

        if(!undo) {
            assert(edges[u][v]>=0);
            assert(edges[u][w]>=0);
            assert(edges[v][w]<=0);
        }
    };

    vector<Triple> packing;
    for(int u=0; u<n; ++u) {
        for(int v=0; v<n; ++v) {
            auto& uv = edges[u][v];
            if(uv<=0) continue;
            for(int w=v+1; w<n && uv>0; ++w) {
                auto& uw = edges[u][w];
                auto& vw = edges[v][w];
                if(uw<=0 || vw>=0) continue;
                packing.push_back({u,v,w, min(uv,min(uw,-vw))});
                assert(edges[u][v]>0 && edges[u][w]>0 && edges[v][w]<0);
                assert(packing.back()[3]>0);
                apply(packing.back());
                if(cost>limit) return cost;
            }
        }
    }

    // local search part
    int num_no_improvement = 0;
    for(int iter=0; iter<INF && cost<=limit && num_no_improvement<5; ++iter) {

        bool has_improved = false;

        for(int i=0; i<size(packing); ++i) {
            // try to replace this triple in packing
            auto [u,v,w,uvw_cost] = packing[i];
            // undo triple
            apply(packing[i], true);
            assert(edges[u][v]>0 && edges[u][w]>0 && edges[v][w]<0);

            optional<Triple> rep_uv;
            for(int x=0; x<n; ++x) { // find other triple using the edge uv
                if(x==w || x==v || x==u) continue;
                if(edges[v][x]==0 || edges[u][x]==0) continue;
                if((edges[v][x]>0)==(edges[u][x]>0)) continue;
                int c = min(edges[u][v], min(abs(edges[u][x]), abs(edges[v][x])));
                Triple t{u,v,x,c};
                if(edges[u][x]<0) swap(t[0],t[1]);
                assert(edges[t[0]][t[1]]>0 && edges[t[0]][t[2]]>0 && edges[t[1]][t[2]]<0);
                if(!rep_uv || rep_uv->at(3)<c) rep_uv = t;
            }

            optional<Triple> rep_uw;
            for(int x=0; x<n; ++x) { // find other triple using the edge uw
                if(x==w || x==v || x==u || (rep_uv && x==rep_uv->at(2))) continue;
                if(edges[w][x]==0 || edges[u][x]==0) continue;
                if((edges[w][x]>0)==(edges[u][x]>0)) continue;
                int c = min(edges[u][w], min(abs(edges[u][x]), abs(edges[w][x])));
                Triple t{u,w,x,c};
                if(edges[u][x]<0) swap(t[0],t[1]);
                assert(edges[t[0]][t[1]]>0 && edges[t[0]][t[2]]>0 && edges[t[1]][t[2]]<0);
                if(!rep_uw || rep_uw->at(3)<c) rep_uw = t;
            }

            optional<Triple> rep_vw;
            for(int x=0; x<n; ++x) { // find other triple using the edge vw
                if(x==w || x==v || x==u) continue;
                if(rep_uv && x==rep_uv->at(2)) continue;
                if(rep_uw && x==rep_uw->at(2)) continue;
                if(edges[w][x]<=0 || edges[v][x]<=0) continue; // both must be present
                int c = min(-edges[v][w], min(edges[v][x], edges[w][x]));
                Triple t{x,v,w,c};
                assert(edges[t[0]][t[1]]>0 && edges[t[0]][t[2]]>0 && edges[t[1]][t[2]]<0);
                if(!rep_vw || rep_vw->at(3)<c) rep_vw = t;
            }

            vector<Triple> reps;
            if(rep_uv) reps.push_back(*rep_uv);
            if(rep_uw) reps.push_back(*rep_uw);
            if(rep_vw) reps.push_back(*rep_vw);
            auto sum = accumulate(begin(reps), end(reps), 0, [](int acc, auto& t){ return acc+t[3]; });
            if(uvw_cost>=sum) { // uvw was good
                apply(packing[i]);
                continue;
            }

            // replace uvw
            has_improved = true;
            sort(begin(reps), end(reps), [](auto& t1, auto& t2){ return t1[3]>t2[3]; });
            packing[i] = reps[0];
            apply(reps[0]);
            for(int k=1; k<size(reps); ++k) {
                packing.push_back(reps[k]);
                apply(reps[k]);
            }
        }

        num_no_improvement = has_improved ? 0 : num_no_improvement+1;

    }

    return cost;
}

Solution ExactSolver::solve_internal(Instance graph, int budget) {

    if(chrono::steady_clock::now()>time_limit) return {};

    // apply reductions
    int num_reduces = 0;
    bool changed = false;
    for (bool repeat = true; repeat; changed |= repeat) {

        //auto lower_bound = packing_lower_bound(graph.edges, budget-graph.spendCost);
        auto lower_bound = packing_local_search_bound(graph.edges, budget-graph.spendCost);
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
        auto r1 = icxReductions(graph, budget);
        if (r1) repeat = true, graph = *r1;
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

    //if(auto r = icxReductions(graph, budget); r)
    //return solveMaybeUnconnected(*r, budget);


    branchingNodes++;
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

    // try merging
    auto mergedSolution = solve_internal(minstance, budget);
    assert(!mergedSolution.worked || mergedSolution.cost <= budget);
    if (mergedSolution.worked)
        return mergedSolution;

    // try permanent deletion
    auto forbiddenSolution = solve_internal(finstance, budget);
    assert(!forbiddenSolution.worked || forbiddenSolution.cost <= budget);
    return forbiddenSolution;
}

Solution ExactSolver::solve_unconnected(Instance graph, int budget) {
    if (budget < graph.spendCost) return {}; // that should probably never happen
    Solution solution;
    solution.worked = true;
    solution.cost = graph.spendCost;
    auto comps = constructConnectedComponents(graph);
    if (size(comps) > 1) numDisconnects++;
    for (auto comp : comps) {
        auto subsolution = solve_internal(comp, budget - solution.cost); // TODO this breaks if solve_internal is not guaranteed to find best solution
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

void ExactSolver::reset_stats() {
    branchingNodes = 0;
    numReducingNodes = 0;
    sumReductions = 0;
    numDisconnects = 0;
    numPrunes = 0;
}

Solution ExactSolver::solve(Instance inst, int budget_limit) {
    // try some reductions for unweighted instances only
    auto isUnweighted = true;
    for(auto& row : inst.edges)
        for(auto val : row)
            isUnweighted &= val==1 || val==-1;
    if(isUnweighted) {
        if(auto opt = thomas(inst); opt) inst = *opt;
        if(auto opt = distance4Reduction(inst); opt) inst = *opt;
    }

    Solution s_comb;
    s_comb.cost = inst.spendCost;

    if(verbose) cout << "cost of initial reductions " << s_comb.cost << endl;
    auto comps = constructConnectedComponents(inst);
    sort(begin(comps), end(comps), [](auto& a, auto& b){ return size(a.edges)<size(b.edges); });
    if(verbose) cout << "split into " << size(comps) << " CCs" << endl;
    for (auto& comp : comps) {

        Solution s;
        auto lower = packing_local_search_bound(comp.edges, INF);
        if(verbose) cout << "start solving CC of size " << size(comp.edges) << " first bound " << lower << endl;

        for (int budget = lower; !s.worked && s_comb.cost+budget<=budget_limit; ++budget) {
            //if(verbose) cout << "===== starting " << budget << " =====" << endl;
            if(verbose) cout << budget << "   \r" << flush;
            s = solve_internal(comp, budget);
            //if(verbose) cout << *this;
            reset_stats();
            //if(verbose) cout << "===== finished " << budget << " =====" << endl << endl;

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

ostream &operator<<(ostream &os, const ExactSolver &rhs) {
    os << "branching nodes: " << rhs.branchingNodes << endl;
    os << "reductions:      " << rhs.numReducingNodes << endl;
    os << "disconnects:     " << rhs.numDisconnects << endl;
    os << "prunes:          " << rhs.numPrunes << endl;
    os << "iters/reduction: " << rhs.sumReductions * 1.0 / rhs.numReducingNodes << endl;
    return os;
}