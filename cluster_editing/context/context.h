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
  bool print_csv = false;
  int seed = 0;
  int num_repititions = 0;
  int num_fruitless_repititions = 0;
};

std::ostream & operator<< (std::ostream& str, const GeneralParameters& params);

struct LabelPropagationRefinerParameters {
  int maximum_lp_iterations = std::numeric_limits<int>::max();
  int activate_all_cliques_after_rounds = std::numeric_limits<int>::max();
  bool random_shuffle_each_round = false;
  NodeOrdering node_order = NodeOrdering::none;
  int min_improvement = 0;
};

std::ostream & operator<< (std::ostream& str, const LabelPropagationRefinerParameters& params);

struct FMParameters {
  int maximum_fm_iterations = std::numeric_limits<int>::max();
  double fraction_of_fruitless_moves = 1.0;
  size_t max_fruitless_moves = std::numeric_limits<size_t>::max();
  size_t num_seed_nodes = 0;
};

std::ostream & operator<< (std::ostream& str, const FMParameters& params);

struct RefinementParameters {
  bool use_lp_refiner = false;
  bool use_boundary_fm_refiner = false;
  bool use_localized_fm_refiner = false;
  LabelPropagationRefinerParameters lp;
  FMParameters fm;
};

std::ostream & operator<< (std::ostream& str, const RefinementParameters& params);

class Context {
 public:
  GeneralParameters general { };
  RefinementParameters refinement { };

  Context() { }

  void computeParameters(const int num_nodes);
};

std::ostream & operator<< (std::ostream& str, const Context& params);

} // namespace cluster_editing
