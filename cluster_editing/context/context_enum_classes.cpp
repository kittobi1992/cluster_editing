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

#include "context_enum_classes.h"

#include "cluster_editing/macros.h"

namespace cluster_editing {

  std::ostream & operator<< (std::ostream& os, const NodeOrdering& order) {
    switch (order) {
      case NodeOrdering::none: return os << "none";
      case NodeOrdering::random_shuffle: return os << "random_shuffle";
      case NodeOrdering::degree_increasing: return os << "degree_increasing";
      case NodeOrdering::degree_decreasing: return os << "degree_decreasing";
        // omit default case to trigger compiler warning for missing cases
    }
    return os << static_cast<uint8_t>(order);
  }

  std::ostream & operator<< (std::ostream& os, const EvoAction& action) {
    switch (action) {
      case EvoAction::INTESIVATE: return os << "INTESIVATE";
      case EvoAction::MUTATE: return os << "MUTATE";
      case EvoAction::UNDEFINED: return os << "UNDEFINED";
        // omit default case to trigger compiler warning for missing cases
    }
    return os << static_cast<uint8_t>(action);
  }

  std::ostream & operator<< (std::ostream& os, const Mutation& mutation) {
    switch (mutation) {
      case Mutation::LARGE_CLIQUE_ISOLATOR: return os << "LARGE_CLIQUE_ISOLATOR";
      case Mutation::LARGE_CLIQUE_WITH_NEIGHBOR_ISOLATOR: return os << "LARGE_CLIQUE_WITH_NEIGHBOR_ISOLATOR";
      case Mutation::RANDOM_NODE_ISOLATOR: return os << "RANDOM_NODE_ISOLATOR";
      case Mutation::RANDOM_NODE_MOVER: return os << "RANDOM_NODE_MOVER";
      case Mutation::CLIQUE_SPLITTER: return os << "CLIQUE_SPLITTER";
      case Mutation::NUM_MUTATIONS: return os << "UNDEFINED";
        // omit default case to trigger compiler warning for missing cases
    }
    return os << static_cast<uint8_t>(mutation);
  }

  NodeOrdering nodeOrderingFromString(const std::string& order) {
    if (order == "none") {
      return NodeOrdering::none;
    } else if(order == "random_shuffle") {
      return NodeOrdering::random_shuffle;
    } else if(order == "degree_increasing") {
      return NodeOrdering::degree_increasing;
    } else if(order == "degree_decreasing") {
      return NodeOrdering::degree_decreasing;
    }
    ERROR("Illegal option: " + order);
    return NodeOrdering::none;
  }

} // namespace cluster_editing