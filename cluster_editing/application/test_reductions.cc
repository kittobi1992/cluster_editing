
#include <algorithm>
#include <iostream>
#include <set>
#include <cassert>

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

bool reduce(Instance& inst, int upper) {
    bool changed = false;
    if(auto opt = complexTwin(inst,true); opt) inst = *opt, changed = true;
    if(auto opt = forcedChoices(inst, upper, false); opt) inst = *opt, changed = true;
    if(auto opt = icxReductions(inst, upper); opt) inst = *opt, changed = true;
    return changed;
}

int real_inst_size(const Instance& inst) {
    int n = size(inst.edges);
    auto compId = connectedComponents(inst.edges);
    int num = 1+ *max_element(begin(compId), end(compId));
    vector isClique(num, true);
    vector compSize(num, 0);
    for(int c : compId) compSize[c]++;
    for(int v=0; v<n; ++v)
        for(int u=v+1; u<n; ++u)
            if(compId[v]==compId[u] && inst.edges[v][u]<0)
                isClique[compId[v]] = false;

    auto inst_size = 0;
    for(int v=0; v<n; ++v)
        inst_size += !isClique[compId[v]];
    return inst_size;
}


int main(int argc, char *argv[]) {

    int total_size = 0;
    for(int i=1; i<200; i+=2) {
        auto inst = load_exact_instance(i);
        cout << "+-"[not_solved.count(i)];
        cout << "instance\t" << i << "\tn=\t" << size(inst.edges);
        auto upper = solve_heuristic(inst).cost;
        if(auto opt = thomas(inst); opt) inst = *opt;
        if(auto opt = distance4Reduction(inst); opt) inst = *opt;
        if(auto opt = simpleTwin(inst); opt) inst = *opt;
        while(reduce(inst,upper));

        if(inst.edges!=load_exact_instance(i).edges) {
            auto left = upper - inst.spendCost - packing_local_search_bound(inst, INF);
            cout << "\t\tn'=\t" << real_inst_size(inst) << "\tleft\t" << left << "\tspent\t" << inst.spendCost << endl;
        } else
            cout << endl;
        total_size += real_inst_size(inst);
    }

    cout << "Total remaining size: " << total_size << endl;

    return 0;
}
