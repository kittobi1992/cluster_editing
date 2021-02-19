#include "context_enum_classes.h"

#include "cluster_editing/macros.h"

namespace cluster_editing {

  std::ostream & operator<< (std::ostream& os, const CoarseningAlgorithm& algo) {
    switch (algo) {
      case CoarseningAlgorithm::do_nothing: return os << "do_nothing";
      case CoarseningAlgorithm::lp_coarsener: return os << "lp_coarsener";
      case CoarseningAlgorithm::UNDEFINED: return os << "UNDEFINED";
        // omit default case to trigger compiler warning for missing cases
    }
    return os << static_cast<uint8_t>(algo);
  }

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

  CoarseningAlgorithm coarseningAlgorithmFromString(const std::string& algo) {
    if (algo == "do_nothing") {
      return CoarseningAlgorithm::do_nothing;
    } else if(algo == "lp_coarsener") {
      return CoarseningAlgorithm::lp_coarsener;
    }
    ERROR("Illegal option: " + algo);
    return CoarseningAlgorithm::UNDEFINED;
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