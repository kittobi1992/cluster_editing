#pragma once

#include <iostream>
#include <string>

namespace cluster_editing {

enum class FMType {
  boundary,
  localized
};

enum class NodeOrdering : uint8_t {
  none = 0,
  random_shuffle = 1,
  degree_increasing = 2,
  degree_decreasing = 3
};

enum class EvoAction : uint8_t {
  INTESIVATE = 0,
  MUTATE = 1,
  UNDEFINED
};

enum class Mutation : uint8_t {
  LARGE_CLIQUE_ISOLATOR = 0,
  LARGE_CLIQUE_WITH_NEIGHBOR_ISOLATOR = 1,
  RANDOM_NODE_ISOLATOR = 2,
  RANDOM_NODE_MOVER = 3,
  CLIQUE_SPLITTER = 4,
  NUM_MUTATIONS = 5
};

std::ostream & operator<< (std::ostream& os, const NodeOrdering& order);

std::ostream & operator<< (std::ostream& os, const EvoAction& order);

std::ostream & operator<< (std::ostream& os, const Mutation& mutation);

NodeOrdering nodeOrderingFromString(const std::string& order);

} // namespace cluster_editing