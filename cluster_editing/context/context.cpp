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


std::ostream & operator<< (std::ostream& str, const LabelPropagationRefinerParameters& params) {
  str << "\n  Label Propagation Refiner Parameters:" << std::endl;
  str << "    Maximum LP Repititions:       " << params.maximum_lp_repititions << std::endl;
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