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
  bool enable_logging = false;
  bool verbose_output = true;
  bool print_result_line = false;
  bool print_csv = false;
  bool print_edits = false;
  bool write_to_file = false;
  bool read_from_file = false;
  int seed = 0;
  int num_repititions = 0;
  int num_fruitless_repititions = 0;
  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  double time_limit = 0.0;
};

std::ostream & operator<< (std::ostream& str, const GeneralParameters& params);

struct EvolutionaryParameters {
  bool enable_detailed_output = false;
  double time_limit = 0;
  int solution_pool_size = 0;
  int evolutionary_steps = 0;
  int initial_lp_iterations = 0;
  int initial_node_swapper_iterations = 0;
  int intensivate_lp_iterations = 0;
  int lp_iterations_after_mutate = 0;
  int clique_remover_iterations = 0;
  int clique_splitter_iterations = 0;
  int node_swapper_iterations = 0;
  int node_swapper_max_cluster_size = 0;
  bool use_random_node_ordering = true;
  int enable_all_mutations_after_steps = 0;
  std::string enabled_mutations = "000000";
  size_t large_clique_threshold = 0;
  float min_clique_isolate_prob = 0.0;
  float max_clique_isolate_prob = 0.0;
  float min_neighbor_clique_isolate_prob = 0.0;
  float max_neighbor_clique_isolate_prob = 0.0;
  float min_node_isolation_prob = 0.0;
  float max_node_isolation_prob = 0.0;
  float min_node_move_prob = 0.0;
  float max_node_move_prob = 0.0;
  float min_clique_split_mutation_prob = 0.0;
  float max_clique_split_mutation_prob = 0.0;
  float min_test_mutation_prob = 0.0;
  float max_test_mutation_prob = 0.0;
};

std::ostream & operator<< (std::ostream& str, const EvolutionaryParameters& params);

struct LocalizedEvolutionaryParameters {
  size_t steps = 0;
  bool run_until_time_limit = false;
  int max_lp_iterations = 0;
  int min_mutations_nodes = 0;
  int max_mutations_nodes = 0;
  float choose_adjacent_mutation_node_prob = 0.0f;
  int max_distance_to_mutation_node = 0;
  int degree_sampling_threshold = 0;
};

std::ostream & operator<< (std::ostream& str, const LocalizedEvolutionaryParameters& params);

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
  bool use_localized_evo = false;
  bool use_lp_refiner = false;
  bool use_boundary_fm_refiner = false;
  bool use_localized_fm_refiner = false;
  EvolutionaryParameters evo;
  LocalizedEvolutionaryParameters localized_evo;
  LabelPropagationRefinerParameters lp;
  FMParameters fm;
};

std::ostream & operator<< (std::ostream& str, const RefinementParameters& params);

class Context {
 public:
  GeneralParameters general { };
  RefinementParameters refinement { };

  Context() { }

  void configureLogging() {
    if ( !general.enable_logging ) {
      general.print_edits = true;
      general.verbose_output = false;
      general.print_result_line = false;
      general.read_from_file = false;
      general.print_csv = false;
      general.write_to_file = false;
      refinement.evo.enable_detailed_output = false;
    }
  }

  void configureAlgorithm(const Graph& graph) {
    if ( refinement.use_boundary_fm_refiner ) {
      refinement.fm.max_fruitless_moves = std::max(
        refinement.fm.fraction_of_fruitless_moves * graph.numNodes(), 10000.0);
    }

    if ( refinement.use_localized_evo ) {
      if ( graph.numEdges() / 2 > 1000000 ) {
        refinement.evo.time_limit *= 1.5;
      }

      if ( graph.numNodes() < 10000 ) {
        refinement.localized_evo.steps /= 100;
      } else if ( graph.numNodes() > 50000 ) {
        refinement.localized_evo.run_until_time_limit = true;
      }
    }
  }

  bool isTimeLimitReached() const {
    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds(end - general.start);
    return elapsed_seconds.count() > general.time_limit;
  }
};

std::ostream & operator<< (std::ostream& str, const Context& params);

} // namespace cluster_editing
