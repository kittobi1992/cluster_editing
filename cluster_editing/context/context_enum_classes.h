#pragma once

#include <iostream>
#include <string>

namespace cluster_editing {

enum class CoarseningAlgorithm {
  do_nothing,
  lp_coarsener,
  UNDEFINED
};

enum class NodeOrdering {
  none,
  random_shuffle,
  degree_increasing,
  degree_decreasing
};

std::ostream & operator<< (std::ostream& os, const CoarseningAlgorithm& algo);

std::ostream & operator<< (std::ostream& os, const NodeOrdering& order);

CoarseningAlgorithm coarseningAlgorithmFromString(const std::string& algo);

NodeOrdering nodeOrderingFromString(const std::string& order);

} // namespace cluster_editing