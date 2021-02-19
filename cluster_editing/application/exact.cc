
#include <iostream>

#include <cluster_editing/exact/instance.h>
#include <cluster_editing/exact/solver.h>

#include <cluster_editing/utils/timer.h>

using namespace std;


int main(int argc, char *argv[]) {

    auto inst = load_exact_instance(179);

    auto& timer = cluster_editing::utils::Timer::instance();
    timer.start_timer("solving", "solving instance");

    ExactSolver solver;
    solver.verbose = true;
    auto solution = solver.solve(inst);
    cout << solver << endl;
    cout << solution.worked << endl;
    cout << "k=" << solution.cost << endl;

    timer.stop_timer("solving");
    cout << timer.get("solving") << endl;

    return 0;
}
