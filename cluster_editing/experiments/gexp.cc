
#include <iostream>
#include <cassert>
#include <map>
#include <fstream>

#include <cluster_editing/exact/instance.h>
#include <cluster_editing/exact/solver.h>
#include <cluster_editing/io/output.h>
#include <cluster_editing/exact/star_bound.h>

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
    Instance inst(n);
    for(auto [u,v] : edges) inst.edges[u][v] = inst.edges[v][u] = 1;
    inst.orig = inst.edges;
    return inst;
}

void runExperiment(int n, int deg, double ple, double T, int d, int s1, int s2, int s3, int rep, ofstream& csv, string type) {
    cout << "starting iteration" << endl;
    cout << type << ' ' << n << ' ' << deg << ' ' << T << ' ' << rep << endl;

    // sampling
    auto alpha = T==0 ? numeric_limits<double>::infinity() : 1/T;
    auto pos = girgs::generatePositions(n,d,s1);
    auto wei = girgs::generateWeights(n,ple,s2);
    girgs::scaleWeights(wei,deg,d,alpha);
    auto edges = girgs::generateEdges(wei, pos, alpha, s3);
    auto inst = fromEdgeList(n, edges);

    // solving
    ExactSolver solver;
    solver.verbose = true;
    solver.time_limit = chrono::steady_clock::now() + chrono::hours(1);
    auto t1 = chrono::steady_clock::now();
    auto sol = solver.solve(inst);
    auto t2 = chrono::steady_clock::now();
    cout << solver;
    cout << "solved (or terminated) after " << std::chrono::duration_cast<chrono::seconds>(t2-t1).count() << 's' << endl;
    cout << endl;

    // printing
    auto upper = solve_heuristic(inst).cost;
    auto lower = star_bound(inst,upper);
    auto [adds,dels] = get_edits(inst.edges, sol);
    csv << type
        << ',' << n
        << ',' << deg
        << ',' << ple
        << ',' << T
        << ',' << d
        << ',' << s1
        << ',' << s2
        << ',' << s3
        << ',' << rep
        << ',' << chrono::duration_cast<chrono::milliseconds>(t2-t1).count()
        << ',' << solver.branchingNodes
        << ',' << solver.root_size
        << ',' << solver.root_gap
        << ',' << solver.root_time
        << ',' << lower
        << ',' << upper
        << ',' << sol.worked
        << ',' << sol.cost
        << ',' << size(adds)
        << ',' << size(dels)
        << endl;
}

int main() {

    // default girg config
    int n = 150;
    int deg = 10;
    double ple = 2.9;
    double T = 0;
    int d = 2;
    int s1 = 12951243, s2 = 613241, s3 = 90192471;

    ofstream csv("girg.csv");
    csv << "type,n,deg,ple,T,d,s1,s2,s3,rep,time,branches,root_size,root_gap,root_time,lower,upper,solved,opt,adds,dels" << endl;

    // varying n
    for(int n=100; n<=200; n+=10) {
        for(int rep=0; rep<10; ++rep) {
            s1 += n ^ deg ^ rep;
            s2 += n ^ deg ^ rep;
            s3 += n ^ deg ^ rep;
            runExperiment(n,deg,ple,T,d,s1,s2,s3,rep,csv,"n");
        }
    }

    // varying deg
    for(int deg : {7,8,9,10,11,12,13}) {
        for(auto T : {0.0, 0.2, 0.4, 0.6, 0.8}) {
            for(int rep=0; rep<10; ++rep) {
                s1 += n ^ deg ^ rep;
                s2 += n ^ deg ^ rep;
                s3 += n ^ deg ^ rep;
                runExperiment(n,deg,ple,T,d,s1,s2,s3,rep,csv,"grid");
            }
        }
    }

    /*
    // print the default girg
    auto pos = girgs::generatePositions(n,d,s1);
    auto wei = girgs::generateWeights(n,ple,s2);
    girgs::scaleWeights(wei,deg,d,alpha);
    auto edges = girgs::generateEdges(wei, pos, alpha, s3);
    auto inst = fromEdgeList(n, edges);

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
    */

    return 0;
}
