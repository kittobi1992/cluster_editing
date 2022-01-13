
#include <iostream>
#include <cassert>
#include <fstream>

#include <cluster_editing/exact/instance.h>
#include <cluster_editing/exact/solver.h>
#include <cluster_editing/io/output.h>
#include <cluster_editing/exact/reductions.h>
#include <cluster_editing/exact/star_bound.h>
#include <cluster_editing/exact/lower_bounds.h>
#include <cluster_editing/exact/thomas.h>

using namespace std;

auto get_edits(const Edges& edges, const Solution& sol) {
    assert(sol.worked);

    int n = size(edges);
    vector<int> cluster_id(n,-1);
    for(int i=0; i<int(size(sol.cliques)); ++i) {
        for(auto v : sol.cliques[i]) {
            assert(cluster_id[v]==-1);
            cluster_id[v] = i;
        }
    }
    for(auto id : cluster_id)
        assert(id != -1);

    vector<pair<int,int>> adds, rems;
    long long cost = 0;
    for(auto u=0; u<n; ++u)
        for(int v=u+1; v<n; ++v)
            if((cluster_id[u]==cluster_id[v]) == (edges[u][v]<0)) {
                cost += abs(edges[u][v]);
                (edges[u][v]<0 ? adds : rems).push_back({u,v});
            }
    assert(cost == sol.cost);
    return pair(adds, rems);
}

template<class F>
Instance reduceWithOne(Instance inst, F&& reduction, chrono::steady_clock::time_point limit) {
    if(auto opt = distance4Reduction(inst); opt) inst = *opt;
    while(true) {
        if(chrono::steady_clock::now()>limit) break;
        auto opt = reduction(inst);
        if(opt) inst = *opt;
        else break;
    }
    return inst;
}

Instance reduceWithAll(Instance inst, int upper, chrono::steady_clock::time_point limit) {
    auto t1 = chrono::steady_clock::now();
    while(true) {
        if(chrono::steady_clock::now()>limit) break;
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
    }
    auto t2 = chrono::steady_clock::now();
    return inst;
}

int main(int, char**) {

    vector reductions{"all reds", "force p3", "force star", "twin simple", "twin complex", "icx", "heavy edge (b)", "heavy edge (s)", "heavy non-edge", "forced single merge"};

    ofstream csv("pace_instances.csv");
    csv << "num,n,m,solved,opt,upper,dels,adds,low_star,low_p3,time,branches,dist4+,n-clean";
    for(auto rule : reductions)
        csv
            << ",after " << rule
            << ",forbs " << rule
            << ",time " << rule
            << ",spend " << rule
            << ",lower after " << rule;
    csv << endl;

    for(int i=1; i<=200; ++i) {
        cout << "\n\ninstance " << i << endl;

        // usual stuff
        auto inst = load_exact_instance(i);
        auto upper = solve_heuristic(inst).cost;
        long m = 0;
        for(auto& row : inst.edges) m += std::count(row.begin(), row.end(), 1);
        assert(m%2==0); m /= 2;
        ExactSolver solver;
        solver.verbose = true;
        solver.time_limit = chrono::steady_clock::now() + chrono::hours(1);
        auto t1 = chrono::steady_clock::now();
        auto sol = solver.solve(inst);
        auto t2 = chrono::steady_clock::now();
        auto [adds, dels] = get_edits(inst.edges, sol);
        cout << solver << endl;

        // dist 4 reduction applies to all
        int initial_n = size(inst.edges);
        if(auto opt = distance4Reduction(inst); opt) inst = *opt;
        if(auto opt = removeCliques(inst); opt) inst = *opt;

        long dist4num = 0;
        for(auto& row : inst.edges) dist4num += std::count(row.begin(), row.end(), -INF);
        dist4num /= 2;

        csv << i
            << ',' << initial_n
            << ',' << m
            << ',' << sol.worked
            << ',' << sol.cost
            << ',' << upper
            << ',' << size(dels)
            << ',' << size(adds)
            << ',' << star_bound(inst,upper)
            << ',' << packing_local_search_bound(inst,upper)
            << ',' << chrono::duration_cast<chrono::milliseconds>(t2-t1).count()
            << ',' << solver.branchingNodes
            << ',' << dist4num
            << ',' << size(inst.edges); // clean n after cliques are removed

        // reduction effectiveness
        for(auto rule : reductions) {
            cout << "reduce with " << rule << endl;

            Instance reduced;
            auto t3 = chrono::steady_clock::now();
            auto limit = t3 + chrono::hours(1);
            if(rule=="all reds"s) reduced = reduceWithAll(inst, upper, limit);
            if(rule=="force p3"s) reduced = reduceWithOne(inst, [upper](auto& a){ return forcedChoices(a, upper); }, limit);
            if(rule== "force star"s) reduced = reduceWithOne(inst, [upper](auto& a){ return forcedChoicesStarBound(a, upper, false); }, limit);
            if(rule== "twin simple"s) reduced = reduceWithOne(inst, [upper](auto& a){ return simpleTwin(a);}, limit);
            if(rule== "twin complex"s) reduced = reduceWithOne(inst, [upper](auto& a){ return complexTwin(a, true);}, limit);
            if(rule== "icx"s) reduced = reduceWithOne(inst, [upper](auto& a){ return icxReductions(a, upper);}, limit);
            if(rule== "heavy edge (b)"s) reduced = reduceWithOne(inst, [upper](auto& a){ return thomas_pairs(a);}, limit);
            if(rule== "heavy edge (s)"s) reduced = reduceWithOne(inst, [upper](auto& a){ return heavy_edge_single_end(a);}, limit);
            if(rule== "heavy non-edge"s) reduced = reduceWithOne(inst, [upper](auto& a){ return heavy_non_edge_single_end(a);}, limit);
            if(rule== "forced single merge"s) reduced = reduceWithOne(inst, [upper](auto& a){ return forcedChoicesSingleMerge(a, upper, false);}, limit);
            auto t4 = chrono::steady_clock::now();

            if(auto opt = removeCliques(reduced); opt) reduced = *opt;

            auto lower_after = reduced.spendCost + meta_lower_bound(reduced,upper);

            long forbs = 0;
            for(auto& row : reduced.edges) forbs += std::count(row.begin(), row.end(), -INF);
            assert(forbs%2==0); forbs /= 2;
            csv << ',' << size(reduced.edges)
               << ',' << forbs
               << ',' << chrono::duration_cast<chrono::milliseconds>(t4-t3).count()
               << ',' << reduced.spendCost
               << ',' << lower_after;
        }

        csv << endl;
    }


    return 0;
}
