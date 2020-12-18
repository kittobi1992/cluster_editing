#pragma once

#include <vector>

#include "cluster_editing/datastructures/graph.h"

namespace cluster_editing::ds {

class GraphFactory {

 using AdjacencyList = std::vector<std::vector<NodeID>>;

 public:
  static Graph construct(const AdjacencyList& adj_list);

 private:
  GraphFactory() { }
};

} // namespace cluster_editing::ds