
#include <algorithm>
#include <iostream>

#include <cluster_editing/exact/instance.h>
#include <cluster_editing/exact/solver.h>

using namespace std;


int main(int argc, char *argv[]) {

    int seconds = 1;
    if(argc>1) seconds = atoi(argv[1]);

    int solved = 0;
    for(int i=1; i<200; i+=2) {
        auto inst = load_exact_instance(i);
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
