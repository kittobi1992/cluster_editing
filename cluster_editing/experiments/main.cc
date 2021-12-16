
#include <iostream>
#include <cassert>
#include <map>
#include <fstream>

#include <cluster_editing/exact/instance.h>
#include <cluster_editing/exact/solver.h>
#include <cluster_editing/io/output.h>
#include <cluster_editing/exact/reductions.h>
#include <cluster_editing/exact/star_bound.h>
#include <cluster_editing/exact/lower_bounds.h>
#include <cluster_editing/exact/thomas.h>
#include <cluster_editing/utils/common_operations.h>

#include <girgs/Generator.h>

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

void make_pdf(
        std::vector<std::vector<double>> pos,
        std::vector<std::pair<int, int>> graph,
        std::string file,
        map<int,string> nodeColors = {},
        map<int,string> nodeLabels = {},
        map<pair<int,int>,string> edgeColors = {},
        bool torus=true) {

    std::ofstream f{file+".dot"};
    //f << "graph girg {\n\toverlap=scale;\n\n";
    f << "graph girg {\n\tscale=2000\n\n";

    if(torus) {
        f << "\tviewport=\"2000,2000\"\n\n";
        auto dist = [](auto a, auto b) { return max(abs(a[0]-b[0]), abs(a[1]-b[1])); };
        auto npos = [&](int v, int i) { return vector{pos[v][0] + (i&1), pos[v][1] + ((i>>1)&1)}; };

        // duplicate nodes 4 times
        int n = size(pos);
        for(int i=1;i<4;++i) {
            for(int v=0; v<n; ++v) {
                pos.push_back(npos(v,i));
                if(nodeColors.count(v)) nodeColors[v+i*n] = nodeColors[v];
                nodeLabels[v+i*n] = nodeLabels.count(v) ? nodeLabels[v] : to_string(v);
            }
        }

        // duplicate edges; only draw if dist <= 0.5
        auto old = graph;
        graph.clear();
        for(auto [u,v] : old) {
            auto it = edgeColors.find(minmax(u,v));
            for(int i=0;i<4;++i) {
                for(int j=0; j<4; ++j) {
                    if(dist(npos(u,i), npos(v,j))>0.5) continue;
                    graph.emplace_back(u + i*n, v + j*n);
                    if(it==end(edgeColors)) continue;
                    edgeColors[minmax(u+i*n,v+j*n)] = it->second;
                }
            }
        }
    }

    f << std::fixed;
    for (int i = 0; i < (int)pos.size(); ++i) {
        f << '\t' << i << " [pos=\"" << pos[i][0] << ',' << pos[i][1] << "\"";
        if(nodeColors.count(i) > 0)
            f << " style =\"filled\" fillcolor=\"" << nodeColors[i] << "\"";
        if(nodeLabels.count(i) > 0)
            f << " label=\"" << nodeLabels[i] << "\"";
        f << "];\n";
    }
    f << '\n';
    for (auto e : graph) {
        auto [a,b] = minmax(e.first,e.second);
        f << '\t' << a << "\t-- " << b;
        if(edgeColors.count({a,b}) > 0)
            f << "\t[color=\"" << edgeColors[{a,b}] << "\"]";
        f << ";\n";
    }
    f << "}\n";
    f.flush();

    // create pdf
    string command = "neato -n -Tpdf " + file + ".dot -o " + file + ".pdf";
    system(command.c_str());
}

Instance fromEdgeList(int n, const vector<pair<int,int>>& edges) {
    vector<vector<int>> res(n, vector<int>(n, -1));
    Instance inst(n);
    for(auto [u,v] : edges) inst.edges[u][v] = inst.edges[v][u] = 1;
    inst.orig = inst.edges;
    return inst;
}

template<class F>
Instance reduceWithOne(Instance inst, F&& reduction) {
    if(auto opt = distance4Reduction(inst); opt) inst = *opt;
    while(true) {
        auto opt = reduction(inst);
        if(opt) inst = *opt;
        else break;
    }
    return inst;
}

Instance reduceWithAll(Instance inst, int upper, bool verbose = false) {
    if(verbose) cout << "solving instance with n=" << size(inst.edges) << endl;
    auto t1 = chrono::steady_clock::now();
    if(verbose) cout << "upper bound " << upper << endl;
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
    if(verbose) {
        cout << "INITIAL REDUCTION FINISHED" << endl;
        cout << "time:  " << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() * 0.001 << " s\n";
        cout << "size:  " << size(inst.edges) << endl;
        int lower = star_bound(inst,upper);
        cout << "lower: " << inst.spendCost + lower << endl;
        cout << "upper: " << upper << endl;
        cout << "gap:   " << upper-lower-inst.spendCost << endl;
    }
    return inst;
}


void girgs_stuff() {
    // default girg config
    int n = 150;
    int deg = 10;
    double ple = 2.9;
    auto alpha = numeric_limits<double>::infinity();
    int d = 2;
    int seed1 = 12951243, seed2 = 613241, seed3 = 90192471;

    auto pos = girgs::generatePositions(n,d,seed1);
    auto wei = girgs::generateWeights(n,ple,seed2);
    girgs::scaleWeights(wei,deg,d,alpha);
    auto edges = girgs::generateEdges(wei, pos, alpha, seed3);

    cout << size(edges) << endl;

    Instance inst = fromEdgeList(n,edges);
    ExactSolver solver;
    solver.verbose = true;
    auto sol= solver.solve(inst);

    sort(begin(sol.cliques), end(sol.cliques), [](auto a, auto b) { return size(a)>size(b); });
    vector which_clique(n,0);
    for(auto i=0u; i<size(sol.cliques); ++i) {
        for(auto v : sol.cliques[i])
            which_clique[v] = i;
    }

    map<int,string> node_col;
    map<pair<int,int>,string> edge_col;

    // color nodes by clique
    for(int i=0;i<n;++i)
        if(which_clique[i]<12)
            node_col[i] = "/paired12/"s + to_string(which_clique[i]+1);

    auto [adds, rems] = get_edits(inst.edges, sol);
    for(auto [u,v] : adds) edges.emplace_back(u,v);

    for(auto [u,v] : adds) edge_col[{u, v}] = "green";
    for(auto [u,v] : rems) edge_col[{u,v}] = "red";

    make_pdf(pos,edges,"graphimage",node_col, {}, edge_col);
}


bool cluster_editing::utils::CommonOperations::dirty = false;

int main(int argc, char *argv[]) {

    vector reductions{"all reds", "force p3", "force star", "twin simple", "twin complex", "icx", "heavy edge (b)", "heavy edge (s)", "heavy non-edge", "forced single merge"};

    ofstream csv("pace_instances.csv");
    csv << "num,n,m,solved,opt,upper,dels,adds,low_star,low_p3,time,branches,dist4+";
    for(auto rule : reductions)
        csv
            << ",after " << rule
            << ",forbs " << rule
            << ",time " << rule
            << ",spend " << rule
            << ",lower after " << rule;
    csv << endl;

    for(int i=1; i<50; ++i) {
        cluster_editing::utils::CommonOperations::dirty = true;

        cout << "\n\ninstance " << i << endl;

        // usual stuff
        auto inst = load_exact_instance(i);
        auto upper = solve_heuristic(inst).cost;
        long m = 0;
        for(auto& row : inst.edges) m += std::count(row.begin(), row.end(), 1);
        assert(m%2==0); m /= 2;
        ExactSolver solver;
        solver.verbose = true;
        solver.time_limit = chrono::steady_clock::now() + chrono::seconds(30);
        auto t1 = chrono::steady_clock::now();
        auto sol = solver.solve(inst);
        auto t2 = chrono::steady_clock::now();
        auto [adds, dels] = get_edits(inst.edges, sol);
        cout << solver << endl;

        // dist 4 reduction applies to all
        if(auto opt = distance4Reduction(inst); opt) inst = *opt;
        long dist4num = 0;
        for(auto& row : inst.edges) dist4num += std::count(row.begin(), row.end(), -INF);
        dist4num /= 2;

        csv << i
            << ',' << size(inst.edges)
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
            << ',' << dist4num;

        // reduction effectiveness
        for(auto rule : reductions) {
            cout << "reduce with " << rule << endl;

            Instance reduced;
            auto t3 = chrono::steady_clock::now();
            if(rule=="all reds"s) reduced = reduceWithAll(inst, upper);
            if(rule=="force p3"s) reduced = reduceWithOne(inst, [upper](auto& a){ return forcedChoices(a, upper); });
            if(rule== "force star"s) reduced = reduceWithOne(inst, [upper](auto& a){ return forcedChoicesStarBound(a, upper, false); });
            if(rule== "twin simple"s) reduced = reduceWithOne(inst, [upper](auto& a){ return simpleTwin(a);});
            if(rule== "twin complex"s) reduced = reduceWithOne(inst, [upper](auto& a){ return complexTwin(a, true);});
            if(rule== "icx"s) reduced = reduceWithOne(inst, [upper](auto& a){ return icxReductions(a, upper);});
            if(rule== "heavy edge (b)"s) reduced = reduceWithOne(inst, [upper](auto& a){ return thomas_pairs(a);});
            if(rule== "heavy edge (s)"s) reduced = reduceWithOne(inst, [upper](auto& a){ return heavy_edge_single_end(a);});
            if(rule== "heavy non-edge"s) reduced = reduceWithOne(inst, [upper](auto& a){ return heavy_non_edge_single_end(a);});
            if(rule== "forced single merge"s) reduced = reduceWithOne(inst, [upper](auto& a){ return forcedChoicesSingleMerge(a, upper, false);});
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
