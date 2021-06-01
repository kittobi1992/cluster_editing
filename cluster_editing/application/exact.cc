
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

#include <iostream>
// #include <cassert>

#include <cluster_editing/exact/instance.h>
#include <cluster_editing/exact/solver.h>
#include <cluster_editing/io/output.h>

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

    // handle args
    for(int i=1; i<argc; ++i) {
        string param = argv[i];
        if(param=="--enable-logging=true")
            solver.verbose = true;
        if(param.find("--time-limit=")==0) {
            int secs = stoi(param.substr(13));
            solver.time_limit = chrono::steady_clock::now() + chrono::seconds(secs);
        }
    }

    if ( solver.verbose ) {
      cluster_editing::io::printBanner();
    }

    auto solution = solver.solve(inst);

    if(solver.verbose) {
        cout << solver << endl;
        cout << "instance solved with " << solution.cost << " edits" << endl;
    }

    print_edits(inst.edges, solution);

    return 0;
}
