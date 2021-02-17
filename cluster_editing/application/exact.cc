
#include <algorithm>
#include <iostream>
#include <numeric>
#include <sstream>
#include <cassert>

#include <cluster_editing/exact/instance.h>
#include <cluster_editing/exact/reductions.h>
#include <cluster_editing/exact/thomas.h>
#include <cluster_editing/exact/solver.h>

#include <cluster_editing/data_path.h>
#include <cluster_editing/utils/timer.h>

using namespace std;

// read graph into matrix
vector<vector<int>> readGraph() {

    std::istringstream sstream;
    auto getline = [&]() {
        std::string line;
        do {
            std::getline(cin, line);
        } while (line[0] == 'c');
        sstream = std::istringstream(line);
    };


    getline();
    std::string skip;
    int n, m;
    sstream >> skip >> skip >> n >> m;

    vector<vector<int>> res(n, vector<int>(n, -1));

    for (int i = 0; i < m; ++i) {
        getline();

        int u, v;
        sstream >> u >> v;
        --u;
        --v;
        res[u][v] = 1;
        res[v][u] = 1;
    }

    return res;
}

void verifySolution(const Edges& edges, const Solution& sol) {
    if(!sol.worked) return;

    int n = size(edges);
    vector<int> cluster_id(n,-1);
    for(int i=0; i<size(sol.cliques); ++i) {
        for(auto v : sol.cliques[i]) {
            assert(cluster_id[v]==-1);
            cluster_id[v] = i;
        }
    }
    for(auto id : cluster_id)
        assert(id!=-1);
    long long cost = 0;
    for(auto u=0; u<n; ++u)
        for(int v=u+1; v<n; ++v)
            if((cluster_id[u]==cluster_id[v]) == (edges[u][v]<0))
                cost += abs(edges[u][v]);
    assert(cost == sol.cost);
}


int main(int argc, char *argv[]) {

    /*
    we can solve
    - 001 (k=3)
    - 003 (k=42)
    - 005 (k=46)
    - 007 (k=46) 226s
    - 009 (k=86) 231s
    - 011 (k=81) 77s
    - 041 (k=184) 15s
    - 057 (k=122)
    - 065 (k=128) 10s
    - 075 (k=132) 61s
    - 077 (k=78)
    - 079 (k=48)
    - 097 (k=95 auch wenn TU Berlin solver 96 sagt)
    - 115 (k=306) 10s
    - 119 (k=184) 9s
    - 127 (k=131 aber TU sagt 130?)
    - 137 (k=16)
    - 153 (k=6)
    - 155 (k=63)
    - 161 (k=333) 88s
    - 173 (k=100)
    - 185 (k=203) 51s
    */
    string file = EXACT_DATA_DIR + "exact179.gr"s;
    if (argc > 1)
        file = argv[1];
    freopen(file.c_str(), "r", stdin);
    cout << "reading file " << file << endl;
    auto edges = readGraph();
    cout << "n = " << size(edges) << endl;

    Instance graph(edges.size());
    graph.edges = edges;

    {
        auto CCs = constructConnectedComponents(graph);
        cout << "CCs before reduction: " << size(CCs) << endl;
        sort(begin(CCs), end(CCs), [](auto& a, auto& b) { return size(a.edges) < size(b.edges); });
        for(auto& comp : CCs) cout << size(comp.edges) << " ";
        cout << endl;
    }

    bool do_again = true;
    while(do_again) {

        auto reduced = thomas(graph);
        if(reduced) graph = *reduced;
        cout << "Spendcost in reduction " << graph.spendCost << endl;

        {
            auto CCs = constructConnectedComponents(graph);
            cout << "CCs after reduction: " << size(CCs) << endl;
            sort(begin(CCs), end(CCs), [](auto& a, auto& b) { return size(a.edges) < size(b.edges); });
            for(auto& comp : CCs) cout << size(comp.edges) << " ";
            cout << endl;
        }

        do_again = reduced.has_value();
    }

    auto distReduced = distance4Reduction(graph);
    if (distReduced) graph = *distReduced;

    auto& timer = cluster_editing::utils::Timer::instance();
    timer.start_timer("solving", "solving instance");

    ExactSolver solver;
    solver.verbose = true;
    auto solution = solver.solve(graph);
    cout << solution.worked << endl;
    cout << "k=" << solution.cost << endl;
    //cout << stats << endl;

    verifySolution(edges, solution);

    timer.stop_timer("solving");
    cout << timer.get("solving") << endl;

    return 0;
}
