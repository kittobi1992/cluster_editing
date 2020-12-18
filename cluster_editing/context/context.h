#pragma once

#include <string>

#include "context_enum_classes.h"

namespace cluster_editing {

struct GeneralParameters {
  std::string graph_filename {};
  std::string config_file {};
  std::string output_file {};
  bool verbose_output = true;
};

std::ostream & operator<< (std::ostream& str, const GeneralParameters& params);

struct CoarseningParameters {
  CoarseningAlgorithm algorithm = CoarseningAlgorithm::UNDEFINED;
};

std::ostream & operator<< (std::ostream& str, const CoarseningParameters& params);

struct RefinementParameters {
  RefinementAlgorithm algorithm = RefinementAlgorithm::UNDEFINED;
};

std::ostream & operator<< (std::ostream& str, const RefinementParameters& params);

class Context {
 public:
  GeneralParameters general { };
  CoarseningParameters coarsening { };
  RefinementParameters refinement { };

  Context() { }
};

std::ostream & operator<< (std::ostream& str, const Context& params);

} // namespace cluster_editing