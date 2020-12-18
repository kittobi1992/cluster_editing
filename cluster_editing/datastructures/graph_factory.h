#pragma once

#include <vector>

#include "cluster_editing/datastructures/graph.h"

namespace cluster_editing::ds {

class GraphFactory {

 public:
  static Graph construct(const std::vector<size_t>& index_array,
                         const std::vector<NodeID> edges);

 private:
  GraphFactory() { }
};

} // namespace cluster_editing::ds