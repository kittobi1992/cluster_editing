
#pragma once

#include <vector>
#include <iostream>

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
    int sumReductions = 0;
    int numDisconnects = 0;
    int numPrunes = 0; // unsolvable leaf

    void reset_stats();
private:
    Solution solve_unconnected(Instance graph, int budget);
    Solution solve_internal(Instance graph, int budget);
};

Solution solve_exact(Instance inst, int budget_limit = INF);

std::ostream &operator<<(std::ostream& os, const ExactSolver& rhs);
