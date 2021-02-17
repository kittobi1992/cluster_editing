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

  CoarseningAlgorithm coarseningAlgorithmFromString(const std::string& algo) {
    if (algo == "do_nothing") {
      return CoarseningAlgorithm::do_nothing;
    } else if(algo == "lp_coarsener") {
      return CoarseningAlgorithm::lp_coarsener;
    }
    ERROR("Illegal option: " + algo);
    return CoarseningAlgorithm::UNDEFINED;
  }

} // namespace cluster_editing