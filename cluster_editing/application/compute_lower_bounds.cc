
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
        if(auto opt = thomas(inst); opt) inst = *opt; // TODO multiple thomas reductions
        if(auto opt = distance4Reduction(inst); opt) inst = *opt;
        auto t1 = chrono::steady_clock::now();
        auto lower_bound = packing_local_search_bound(inst, INF);
        auto t2 = chrono::steady_clock::now();
        auto dur = chrono::duration_cast<chrono::milliseconds>(t2-t1).count() * 0.001;
        cout << "lower bound for\t" << i << '\t' << inst.spendCost + lower_bound << "\tin " << dur << " seconds" << endl;
    }

    return 0;
}
