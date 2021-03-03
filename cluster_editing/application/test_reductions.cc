
#include <algorithm>
#include <iostream>
#include <set>

#include <cluster_editing/exact/instance.h>
#include <cluster_editing/exact/lower_bounds.h>
#include <cluster_editing/exact/reductions.h>
#include <cluster_editing/exact/thomas.h>
#include <cluster_editing/exact/solver.h>

using namespace std;

set not_solved {
    43,
    45,
    51,
    53,
    69,
    81,
    83,
    91,
    93,
    99,
    101,
    103,
    141,
    167,
    169,
    179,
    181,
    183,
    187,
    191,
    193,
    195,
    197,
};


int main(int argc, char *argv[]) {

    for(int i=1; i<200; i+=2) {
        auto inst = load_exact_instance(i);
        cout << "+-"[not_solved.count(i)];
        cout << "instance\t" << i << "\tn=\t" << size(inst.edges);
        auto upper = solve_heuristic(inst).cost;
        //if(auto opt = thomas(inst); opt) inst = *opt;
        //if(auto opt = distance4Reduction(inst); opt) inst = *opt;
        //if(auto opt = simpleTwin(inst); opt) inst = *opt;
        if(auto opt = complexTwin(inst,true); opt) inst = *opt;
        //if(auto opt = forcedChoices(inst, upper, false); opt) inst = *opt;
        if(inst.edges!=load_exact_instance(i).edges)
            cout << "\t\tn'=\t" << size(inst.edges) << "\tspent\t" << inst.spendCost << endl;
        else
            cout << endl;
    }

    return 0;
}
