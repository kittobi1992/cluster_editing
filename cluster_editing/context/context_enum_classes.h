#pragma once

#include <iostream>
#include <string>

namespace cluster_editing {

enum class CoarseningAlgorithm {
  do_nothing,
  lp_coarsener,
  UNDEFINED
};

std::ostream & operator<< (std::ostream& os, const CoarseningAlgorithm& algo);

CoarseningAlgorithm coarseningAlgorithmFromString(const std::string& algo);

} // namespace cluster_editing