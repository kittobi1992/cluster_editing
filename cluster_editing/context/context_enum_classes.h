#pragma once

#include <iostream>
#include <string>

namespace cluster_editing {

enum class CoarseningAlgorithm {
  do_nothing,
  lp_coarsener,
  UNDEFINED
};

enum class NodeOrdering : uint8_t {
  none = 0,
  random_shuffle = 1,
  degree_increasing = 2,
  degree_decreasing = 3
};

std::ostream & operator<< (std::ostream& os, const CoarseningAlgorithm& algo);

std::ostream & operator<< (std::ostream& os, const NodeOrdering& order);

CoarseningAlgorithm coarseningAlgorithmFromString(const std::string& algo);

NodeOrdering nodeOrderingFromString(const std::string& order);

} // namespace cluster_editing