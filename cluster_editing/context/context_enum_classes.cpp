#include "context_enum_classes.h"

#include "cluster_editing/macros.h"

namespace cluster_editing {

  std::ostream & operator<< (std::ostream& os, const CoarseningAlgorithm& algo) {
    switch (algo) {
      case CoarseningAlgorithm::do_nothing: return os << "do_nothing";
      case CoarseningAlgorithm::UNDEFINED: return os << "UNDEFINED";
        // omit default case to trigger compiler warning for missing cases
    }
    return os << static_cast<uint8_t>(algo);
  }

  std::ostream & operator<< (std::ostream& os, const RefinementAlgorithm& algo) {
    switch (algo) {
      case RefinementAlgorithm::do_nothing: return os << "do_nothing";
      case RefinementAlgorithm::UNDEFINED: return os << "UNDEFINED";
        // omit default case to trigger compiler warning for missing cases
    }
    return os << static_cast<uint8_t>(algo);
  }

  CoarseningAlgorithm coarseningAlgorithmFromString(const std::string& algo) {
    if (algo == "do_nothing") {
      return CoarseningAlgorithm::do_nothing;
    }
    ERROR("Illegal option: " + algo);
    return CoarseningAlgorithm::UNDEFINED;
  }

  RefinementAlgorithm refinementAlgorithmFromString(const std::string& algo) {
    if (algo == "do_nothing") {
      return RefinementAlgorithm::do_nothing;
    }
    ERROR("Illegal option: " + algo);
    return RefinementAlgorithm::UNDEFINED;
  }

} // namespace cluster_editing