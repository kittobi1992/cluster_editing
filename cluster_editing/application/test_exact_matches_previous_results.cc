#include <algorithm>
#include <iostream>
#include <map>
#include <iomanip>

#include <cluster_editing/exact/instance.h>
#include <cluster_editing/exact/reductions.h>
#include <cluster_editing/exact/solver.h>
#include <cluster_editing/utils/timer.h>


struct Result {
    int nr;
    int k;
    bool solved;
};
constexpr Result results[] = {
        {1,   3,     true},
        {3,   42,    true},
        {5,   46,    true},
        {7,   86,    true},
        {9,   90,    true},
        {11,  81,    true},
        {13,  181,   true},
        {15,  164,   true},
        {17,  236,   true},
        {19,  298,   true},
        {21,  322,   true},
        {23,  281,   true},
        {25,  439,   true},
        {27,  432,   true},
        {29,  509,   true},
        {31,  285,   true},
        {33,  672,   true},
        {35,  385,   true},
        {37,  703,   true},
        {39,  665,   true},
        {41,  184,   true},
        {43,  899,   false},
        {45,  1066,  false},
        {47,  749,   true},
        {49,  854,   true},
        {51,  1219,  false},
        {53,  1389,  false},
        {55,  1410,  true},
        {57,  122,   true},
        {59,  976,   true},
        {61,  116,   true},
        {63,  665,   true},
        {65,  128,   true},
        {67,  792,   true},
        {69,  1562,  false},
        {71,  2131,  true},
        {73,  1612,  true},
        {75,  132,   true},
        {77,  78,    true},
        {79,  48,    true},
        {81,  1882,  false},
        {83,  2243,  false},
        {85,  1047,  true},
        {87,  3243,  true},
        {89,  534,   true},
        {91,  2419,  false},
        {93,  3373,  false},
        {95,  358,   true},
        {97,  95,    true},
        {99,  644,   false},
        {101, 3490,  false},
        {103, 3658,  false},
        {105, 5320,  true},
        {107, 274,   true},
        {109, 2436,  true},
        {111, 2412,  true},
        {113, 999,   true},
        {115, 306,   true},
        {117, 2777,  true},
        {119, 184,   true},
        {121, 742,   true},
        {123, 416,   true},
        {125, 1413,  true},
        {127, 130,   true},
        {129, 2211,  true},
        {131, 386,   true},
        {133, 4836,  true},
        {135, 772,   true},
        {137, 16,    true},
        {139, 5514,  true},
        {141, 682,   false},
        {143, 1475,  true},
        {145, 1761,  true},
        {147, 3975,  true},
        {149, 5807,  true},
        {151, 1258,  true},
        {153, 6,     true},
        {155, 63,    true},
        {157, 1117,  true},
        {159, 1483,  true},
        {161, 333,   true},
        {163, 1538,  true},
        {165, 6094,  true},
        {167, 7235,  false},
        {169, 7864,  false},
        {171, 2901,  true},
        {173, 100,   true},
        {175, 557,   true},
        {177, 4089,  true},
        {179, 620,   false},
        {181, 1445,  false},
        {183, 2609,  false},
        {185, 203,   true},
        {187, 3654,  false},
        {189, 8563,  true},
        {191, 15757, false},
        {193, 20403, false},
        {195, 2068,  false},
        {197, 14841, false}
};

int main(int argc, char *argv[]) {
    int time_limit = 60;
    if (argc > 1) time_limit = atoi(argv[1]);

    std::map<int, Result> results_map;
    for (auto result : results) {
        results_map[result.nr] = result;
    }

    for (auto[nr, result] : results_map) {
        auto inst = load_exact_instance(nr);

        ExactSolver solver;
        solver.time_limit = std::chrono::steady_clock::now() + std::chrono::seconds(time_limit);
        solver.verbose = false;

        auto &timer = cluster_editing::utils::Timer::instance();
        timer.start_timer("solving", "solving instance");
        auto solution = solver.solve(inst);
        timer.stop_timer("solving");

        std::cout << std::fixed << std::setw(3) << nr << "\t\t";
        std::cout << "current: k=" << std::fixed << std::setw(10) << solution.cost << " solved=" << std::boolalpha
                  << solution.worked << "\t";
        std::cout << "previous: k=" << std::fixed << std::setw(10) << result.k << " solved=" << std::boolalpha
                  << result.solved << "\t\t";
        std::cout << "took=" << timer.get("solving") << std::endl;

        if (solution.worked && result.solved && solution.cost != result.k) {
            throw std::runtime_error("mismatch");
        }
        if (solution.worked && !result.solved && solution.cost < result.k) {
            throw std::runtime_error("mismatch");
        }
    }

    return 0;
}
