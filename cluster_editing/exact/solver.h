
#pragma once

#include <vector>
#include <iostream>
#include <chrono>

#include "instance.h"

struct Solution {
    bool worked = false;
    int cost = INF;
    std::vector<std::vector<int>> cliques;
};

class ExactSolver {
public:
    Solution solve(Instance inst, int budget_limit = INF);

    // stats
    bool verbose = false;
    int branchingNodes = 0;
    int numReducingNodes = 0;
    int numDisconnects = 0;
    int redForcedStar = 0;
    int redForced = 0;
    int redTwin = 0;
    int redTwin2 = 0;
    int redICX = 0;
    int redThomas2 = 0;
    int redHeavyEdge = 0;
    int redHeavyNonEdge = 0;
    int numPrunes = 0; // unsolvable leaf

    std::chrono::steady_clock::time_point time_limit = std::chrono::steady_clock::time_point::max();
    Solution solve_internal(Instance graph, int budget);
};

Solution solve_exact(Instance inst, int budget_limit = INF, int time_limit = INF);

// TODO works only for unweighted instances for now
Solution solve_heuristic(const Instance& inst);

std::ostream &operator<<(std::ostream& os, const ExactSolver& rhs);
