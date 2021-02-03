
#include <algorithm>
#include <iostream>
#include <numeric>
#include <sstream>
#include <cassert>

#include <cluster_editing/exact/instance.h>
#include <cluster_editing/exact/reductions.h>
#include <cluster_editing/exact/thomas.h>
#include <cluster_editing/exact/solver.h>

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



int main(int argc, char *argv[]) {

    /*
    we can solve
    - 001 (k=3)
    - 065 (k=128) longest with 1s
    - 077 (k=78)
    - 137 (k=16)
    - 153 (k=6)
    - 155 (k=63)
    - 173 (k=100)
    */
    string file = "../../../instances/exact/exact179.gr";
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


    auto solution = solveMaybeUnconnected(graph, INF, true);
    cout << solution.worked << endl;
    cout << "k=" << solution.cost << endl;
    //cout << stats << endl;
    return 0;
}
