/*******************************************************************************
 * This file is part of KaPoCE.
 *
 * KaPoCE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaPoCE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaPoCE.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

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
