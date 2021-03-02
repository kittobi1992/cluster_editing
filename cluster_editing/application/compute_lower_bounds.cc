
#include <algorithm>
#include <iostream>

#include <cluster_editing/exact/instance.h>
#include <cluster_editing/exact/lower_bounds.h>
#include <cluster_editing/exact/reductions.h>
#include <cluster_editing/exact/thomas.h>
#include <cluster_editing/exact/solver.h>

using namespace std;

int main(int argc, char *argv[]) {

    for(int i=1; i<200; i+=2) {
        auto inst = load_exact_instance(i);
        cout << "instance " << i << " of size " << size(inst.edges) << endl;
        cout << "Lower " << packing_local_search_bound(inst, INF) << endl;
        cout << "Upper " << solve_heuristic(inst).cost << endl;
        if(auto opt = thomas(inst); opt) inst = *opt; // TODO multiple thomas reductions
        if(auto opt = distance4Reduction(inst); opt) inst = *opt;
        if(auto opt = simpleNeighbor(inst); opt) inst = *opt;
        if(auto opt = forcedChoices(inst, solve_heuristic(inst).cost, true); opt) inst = *opt;
        cout << "After Reductions n=" << size(inst.edges) << endl;
        auto lower = packing_local_search_bound(inst, INF);
        cout << "Spent         " << inst.spendCost << endl;
        cout << "Lower         " << lower << endl;
        cout << "Spend + Lower " << lower+inst.spendCost << endl;
        cout << endl;
    }

    return 0;
}
