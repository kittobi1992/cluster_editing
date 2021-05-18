#pragma once

#include <string>
#include <limits>

#include "cluster_editing/definitions.h"
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
  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  double time_limit = 0.0;
};

std::ostream & operator<< (std::ostream& str, const GeneralParameters& params);

struct EvolutionaryParameters {
  bool enable_detailed_output = false;
  int solution_pool_size = 0;
  int evolutionary_steps = 0;
  int initial_lp_iterations = 0;
  int intensivate_lp_iterations = 0;
  int lp_iterations_after_mutate = 0;
  bool use_random_node_ordering = true;
  std::string enabled_mutations = "000";
  float min_node_isolation_prob = 0.0;
  float max_node_isolation_prob = 0.0;
  float min_node_move_prob = 0.0;
  float max_node_move_prob = 0.0;
  float min_test_mutation_prob = 0.0;
  float max_test_mutation_prob = 0.0;
  float random_prob_selection_prob = 0.0;
};

std::ostream & operator<< (std::ostream& str, const EvolutionaryParameters& params);

struct LabelPropagationRefinerParameters {
  int maximum_lp_iterations = std::numeric_limits<int>::max();
  bool random_shuffle_each_round = false;
  NodeOrdering node_order = NodeOrdering::none;
  NodeID rating_map_degree_threshold = 0;
  int min_target_edit_distance = 0;
  int min_improvement = 0;
  size_t early_exit_window = 0;
};

std::ostream & operator<< (std::ostream& str, const LabelPropagationRefinerParameters& params);

struct FMParameters {
  int maximum_fm_iterations = std::numeric_limits<int>::max();
  double fraction_of_fruitless_moves = 1.0;
  size_t max_fruitless_moves = std::numeric_limits<size_t>::max();
  size_t num_seed_nodes = 0;
  int min_improvement = 0;
  size_t early_exit_window = 0;
};

std::ostream & operator<< (std::ostream& str, const FMParameters& params);

struct RefinementParameters {
  bool use_evo = false;
  bool use_lp_refiner = false;
  bool use_boundary_fm_refiner = false;
  bool use_localized_fm_refiner = false;
  EvolutionaryParameters evo;
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

  bool isTimeLimitReached() const {
    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds(end - general.start);
    return elapsed_seconds.count() > general.time_limit;
  }
};

std::ostream & operator<< (std::ostream& str, const Context& params);

} // namespace cluster_editing
