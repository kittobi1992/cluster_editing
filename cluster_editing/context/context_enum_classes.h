#pragma once

#include <iostream>
#include <string>

namespace cluster_editing {

enum class NodeOrdering {
  none,
  random_shuffle,
  degree_increasing,
  degree_decreasing
};

std::ostream & operator<< (std::ostream& os, const NodeOrdering& order);

NodeOrdering nodeOrderingFromString(const std::string& order);

} // namespace cluster_editing