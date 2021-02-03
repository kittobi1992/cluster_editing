
#include <vector>

#include "instance.h"

struct Solution {
    bool worked = false;
    int cost = INF;
    std::vector<std::vector<int>> cliques;
};

Solution solveMaybeUnconnected(Instance graph, int budget, bool highL = false);
