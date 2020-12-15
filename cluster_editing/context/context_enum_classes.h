#pragma once

#include <iostream>
#include <string>

namespace cluster_editing {

enum class CoarseningAlgorithm {
  do_nothing,
  UNDEFINED
};

enum class RefinementAlgorithm {
  do_nothing,
  UNDEFINED
};

std::ostream & operator<< (std::ostream& os, const CoarseningAlgorithm& algo);

std::ostream & operator<< (std::ostream& os, const RefinementAlgorithm& algo);

CoarseningAlgorithm coarseningAlgorithmFromString(const std::string& algo);

RefinementAlgorithm refinementAlgorithmFromString(const std::string& algo);

} // namespace cluster_editing