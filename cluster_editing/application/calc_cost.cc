
#include <algorithm>
#include <iostream>
#include <cassert>
#include <fstream>

#include <cluster_editing/exact/instance.h>
#include <cluster_editing/exact/reductions.h>
#include <cluster_editing/exact/thomas.h>
#include <cluster_editing/exact/solver.h>

#include <cluster_editing/data_path.h>
#include <cluster_editing/utils/timer.h>

using namespace std;

// read graph into matrix

vector<vector<int>> read_graph(int num) {
    auto suf = to_string(num);
    while(size(suf)<3) suf = '0'+suf;
    auto file_name = EXACT_DATA_DIR + ("exact" + suf + ".gr");
    ifstream in(file_name);

    std::istringstream sstream;
    auto getline = [&]() {
        std::string line;
        do {
            std::getline(in, line);
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

    int seconds = 1;
    if(argc>1) seconds = atoi(argv[1]);

    int solved = 0;
    for(int i=1; i<200; i+=2) {
        auto edges = read_graph(i);
        Instance inst(size(edges));
        inst.edges = edges;
        auto t1 = chrono::steady_clock::now();
        auto s = solve_exact(inst, INF, 1'000*seconds);
        auto t2 = chrono::steady_clock::now();
        auto dur = chrono::duration_cast<chrono::milliseconds>(t2-t1).count() * 0.001;
        if(s.worked) {
            solved++;
            cout << "solved " << i << " in " << dur << " seconds" << endl;
        }
    }

    cout << "total: " << solved << endl;

    return 0;
}
