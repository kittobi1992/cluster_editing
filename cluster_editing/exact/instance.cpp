
#include <cluster_editing/exact/instance.h>

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>

#include <cluster_editing/data_path.h>

using namespace std;

Instance load_exact_instance(int num) {
    assert(1<=num && num<=199 && num%2);
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

    Instance inst(n);
    inst.edges = res;
    inst.orig = res;
    return inst;
}

int forbiddenEdges(const Instance &inst) {
    int res = 0;
    for(auto& row : inst.edges)
        for(auto val : row)
            res += (val==-INF);
    return res/2;
}

Instance remove_nodes(const Instance &inst, const vector<int>& nodes) {
    int old_n = size(inst.edges);

    vector<bool> keep(old_n,true);
    for(auto v : nodes) keep[v] = false;

    vector<int> rem_nodes;
    for(int v=0; v<old_n; ++v) if(keep[v]) rem_nodes.push_back(v);

    Instance res(size(rem_nodes));
    res.done_clusters = inst.done_clusters;
    res.spendCost = inst.spendCost;

    for(int i=0; i<size(rem_nodes); ++i) {
        res.idmap[i] = inst.idmap[rem_nodes[i]];
        for(int j=0; j<size(rem_nodes); ++j) {
            res.edges[i][j] = inst.edges[rem_nodes[i]][rem_nodes[j]];
            res.orig[i][j] = inst.orig[rem_nodes[i]][rem_nodes[j]];
        }
    }

    return res;
}
