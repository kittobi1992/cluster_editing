#pragma once

#include <string>
#include <limits>

#include "context_enum_classes.h"

namespace cluster_editing {

struct GeneralParameters {
  std::string graph_filename {};
  std::string config_file {};
  std::string output_file {};
  bool verbose_output = true;
  bool print_result_line = false;
  int seed = 0;
  bool use_multilevel = false;
  int num_repititions = 0;
};

std::ostream & operator<< (std::ostream& str, const GeneralParameters& params);

struct CoarseningParameters {
  CoarseningAlgorithm algorithm = CoarseningAlgorithm::UNDEFINED;
};

std::ostream & operator<< (std::ostream& str, const CoarseningParameters& params);

struct LabelPropagationRefinerParameters {
  int maximum_lp_iterations = std::numeric_limits<int>::max();
  int activate_all_cliques_after_rounds = std::numeric_limits<int>::max();
};

std::ostream & operator<< (std::ostream& str, const LabelPropagationRefinerParameters& params);

struct RefinementParameters {
  bool use_lp_refiner = false;
  LabelPropagationRefinerParameters lp;
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
