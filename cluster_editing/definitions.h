#pragma once

#include <chrono>

#include "cluster_editing/datastructures/graph.h"

namespace cluster_editing {

using Graph = ds::Graph;
using HighResClockTimepoint = std::chrono::time_point<std::chrono::high_resolution_clock>;

} // namespace cluster_editing