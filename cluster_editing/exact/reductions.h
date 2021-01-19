
#include "instance.h"

std::vector<int> connectedComponents(const Edges& graph);

std::vector<Instance> constructConnectedComponents(const Instance& graph);

// contracts the edge u,v
Instance merge(const Instance& inst, int u, int v);