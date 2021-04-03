
#include <algorithm>
#include <iostream>
#include <set>
#include <fstream>
#include <iomanip>

#include <cluster_editing/exact/instance.h>
#include <cluster_editing/exact/reductions.h>
#include <cluster_editing/exact/thomas.h>
#include <cluster_editing/exact/solver.h>
#include <cluster_editing/exact/star_bound.h>

using namespace std;

set not_solved {
    43, // did it once; 984 is solution
    51,
    69,
    83,
    91,
    93,
    101,
    103,
    105, // was solved by threshold editing; needs faster branching :)
    167,
    169,
    179,
    181,
    183,
    191,
    195,
    197,
};

bool reduce(Instance& inst, int upper) {
    string applied = "";
    if(empty(applied)) if(auto opt = force_small_components(inst); opt) inst = *opt, applied = "small clean up";
    if(empty(applied)) if(auto opt = forcedChoicesStarBound(inst, upper, false); opt) inst = *opt, applied = "force star";
    if(empty(applied)) if(auto opt = forcedChoices(inst, upper); opt) inst = *opt, applied = "force p3";
    if(empty(applied)) if(auto opt = simpleTwin(inst); opt) inst = *opt, applied = "twin simple";
    if(empty(applied)) if(auto opt = complexTwin(inst,true); opt) inst = *opt, applied = "twin complex";
    if(empty(applied)) if(auto opt = icxReductions(inst, upper); opt) inst = *opt, applied = "icx";
    if(empty(applied)) if(auto opt = thomas_pairs(inst); opt) inst = *opt, applied = "heavy edge (b)";
    if(empty(applied)) if(auto opt = heavy_edge_single_end(inst); opt) inst = *opt, applied = "heavy edge (s)";
    if(empty(applied)) if(auto opt = heavy_non_edge_single_end(inst); opt) inst = *opt, applied = "heavy non-edge";
    //if(empty(applied)) if(auto opt = forcedChoicesStarBound(inst, upper, true, 60); opt) inst = *opt, applied = "force star (slow)";
    if(empty(applied)) if(auto opt = forcedChoicesSingleMerge(inst, upper, false); opt) inst = *opt, applied = "forced single merge";
    if(empty(applied)) if(auto opt = weightedKernel(inst); opt) inst = *opt, applied = "weighted kernel";
    //if(applied!="") cout << "reduced to n=" << size(inst.edges) << " -INFs=" << forbiddenEdges(inst) << " spent=" << inst.spendCost << " with " << applied << endl;
    return !empty(applied);
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

    bool write_log_csv = true;
    ofstream log_file("reduction_results.csv");
    if(write_log_csv) log_file << "can_solve,instance,n,size after reductions,upper-lower\n";

    int total_size = 0;
    int total_left = 0;
    double total_time = 0;
    for(int i=1; i<200; i+=2) {
        auto inst = load_exact_instance(i);
        int old_n = size(inst.edges);
        cout << "+-"[not_solved.count(i)];
        cout << "instance\t" << i << "\tn=\t" << size(inst.edges);
        auto t1 = chrono::steady_clock::now();
        auto upper = solve_heuristic(inst).cost;
        if(auto opt = thomas(inst); opt) inst = *opt;
        if(auto opt = distance4Reduction(inst); opt) inst = *opt;
        if(auto opt = simpleTwin(inst); opt) inst = *opt;
        while(reduce(inst,upper));
        auto t2 = chrono::steady_clock::now();

        auto left = upper - inst.spendCost - star_bound(inst, INF);
        auto real_size = real_inst_size(inst);
        auto time = chrono::duration_cast<chrono::duration<double>>(t2-t1).count();
        total_left += left;
        total_size += real_size;
        total_time += time;
        cout << "\t\tn'=\t" << real_size << "\tleft\t" << left << "\ttime\t" << time << endl;
        if(write_log_csv) log_file
                    << "+-"[not_solved.count(i)]
                    << ',' << i
                    << ',' << old_n
                    << ',' << real_size
                    << ',' << left
                    << endl;
    }

    cout << "Total remaining size:   " << total_size << endl;
    cout << "Total upper-lower diff: " << total_left << endl;
    cout << "Total time:             " << total_time << endl;

    /*
    Total remaining size:   7176
    Total upper-lower diff: 18199
    Total time:             319.864
     */

    return 0;
}
