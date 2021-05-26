
#include <iostream>
// #include <cassert>

#include <cluster_editing/exact/instance.h>
#include <cluster_editing/exact/solver.h>

using namespace std;

// we redefine assert for the final submission
// should never happen but we like to get TLE instead of WA
#define assert(cond) if(!(cond)) while(1);

void print_edits(const Edges& edges, const Solution& sol) {
    assert(sol.worked);

    int n = size(edges);
    vector<int> cluster_id(n,-1);
    for(int i=0; i<size(sol.cliques); ++i) {
        for(auto v : sol.cliques[i]) {
            assert(cluster_id[v]==-1);
            cluster_id[v] = i;
        }
    }
    for(auto id : cluster_id)
        assert(id != -1);
    long long cost = 0;
    for(auto u=0; u<n; ++u)
        for(int v=u+1; v<n; ++v)
            if((cluster_id[u]==cluster_id[v]) == (edges[u][v]<0)) {
                cost += abs(edges[u][v]);
                cout << u+1 << ' ' << v+1 << '\n';
            }
    assert(cost == sol.cost);
}


int main(int argc, char *argv[]) {

    auto inst = load_exact_instance(); // read from stdin

    ExactSolver solver;
    auto solution = solver.solve(inst);

    print_edits(inst.edges, solution);

    return 0;
}
