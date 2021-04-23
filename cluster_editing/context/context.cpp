#include <iostream>

#include "context.h"

namespace cluster_editing {

std::ostream & operator<< (std::ostream& str, const GeneralParameters& params) {
  str << "General Parameters:" << std::endl;
  str << "  Graph:                         " << params.graph_filename << std::endl;
  str << "  Config File:                   " << params.config_file << std::endl;
  str << "  Output File:                   " << params.output_file << std::endl;
  str << "  Verbose Output:                " << std::boolalpha << params.verbose_output << std::endl;
  str << "  Seed:                          " << params.seed << std::endl;
  str << "  Number of Repititions:         " << params.num_repititions << std::endl;
  str << "  Number Fruitless Repititions:  " << params.num_fruitless_repititions << std::endl;
  return str;
}


std::ostream & operator<< (std::ostream& str, const EvolutionaryParameters& params) {
  str << "\n  Evolutionary Parameters:" << std::endl;
  str << "    Enable Detailed Output:      " << std::boolalpha << params.enable_detailed_output << std::endl;
  str << "    Solution Pool Size:          " << params.solution_pool_size << std::endl;
  str << "    Evolutionary Steps:          " << params.evolutionary_steps << std::endl;
  str << "    Initial LP Iterations:       " << params.initial_lp_iterations << std::endl;
  str << "    Intensivate LP Iterations:   " << params.intensivate_lp_iterations << std::endl;
  str << "    LP Iterations After Mutate:  " << params.lp_iterations_after_mutate << std::endl;
  str << "    Use Random Node Ordering:    " << std::boolalpha << params.use_random_node_ordering << std::endl;
  str << "    Enabled Mutations:           " << params.enabled_mutations << std::endl;
  str << "    Intensivate Probability:     " << params.intensivate_prob << std::endl;
  str << "    Muate Probability:           " << params.mutate_prob << std::endl;
  str << "    Combine Probability:         " << (1.0 - params.intensivate_prob - params.mutate_prob) << std::endl;
  str << "    Large Clique Threshold:      " << params.large_clique_threshold << std::endl;
  str << "    Clique Isolation Prob.:      " << params.clique_isolate_prob << std::endl;
  str << "    Neigh.Clique Isolation Prob: " << params.neighbor_clique_isolate_prob << std::endl;
  str << "    Node Isolation Prob.:        " << params.node_isolation_prob << std::endl;
  str << "    Node Move Prob.:             " << params.node_move_prob << std::endl;
  return str;
}

std::ostream & operator<< (std::ostream& str, const LabelPropagationRefinerParameters& params) {
  str << "\n  Label Propagation Refiner Parameters:" << std::endl;
  str << "    Maximum LP Iterations:       " << params.maximum_lp_iterations << std::endl;
  str << "    Random Shuffle each Round:   " << std::boolalpha << params.random_shuffle_each_round << std::endl;
  str << "    Node Ordering:               " << params.node_order << std::endl;
  str << "    Min Improvement:             " << params.min_improvement << std::endl;
  str << "    Early Exit Window:           " << params.early_exit_window << std::endl;
  return str;
}

std::ostream & operator<< (std::ostream& str, const FMParameters& params) {
  str << "\n  FM Refiner Parameters:" << std::endl;
  str << "    Maximum FM Iterations:       " << params.maximum_fm_iterations << std::endl;
  str << "    Max. Num Fruitless Moves:    " << params.max_fruitless_moves << std::endl;
  str << "    Num Seed Nodes:              " << params.num_seed_nodes << std::endl;
  str << "    Min Improvement:             " << params.min_improvement << std::endl;
  str << "    Early Exit Window:           " << params.early_exit_window << std::endl;
  return str;
}


std::ostream & operator<< (std::ostream& str, const RefinementParameters& params) {
  str << "Refinement Parameters:" << std::endl;
  if ( params.use_evo ) {
    str << params.evo;
  }
  if ( params.use_lp_refiner ) {
    str << params.lp;
  }
  if ( params.use_boundary_fm_refiner || params.use_localized_fm_refiner ) {
    str << params.fm;
  }
  return str;
}

std::ostream & operator<< (std::ostream& str, const Context& context) {
  str << "*******************************************************************************\n"
      << "*                         Cluster Editing Context                             *\n"
      << "*******************************************************************************\n"
      << context.general
      << "-------------------------------------------------------------------------------\n"
      << context.refinement
      << "-------------------------------------------------------------------------------\n";
  return str;
}

void Context::computeParameters(const int num_nodes) {
  if ( refinement.use_boundary_fm_refiner ) {
    refinement.fm.max_fruitless_moves = std::max(
      refinement.fm.fraction_of_fruitless_moves * num_nodes, 10000.0);
  }
}

} // namespace cluster_editing