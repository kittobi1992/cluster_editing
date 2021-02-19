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
  str << "  Use Multilevel:                " << std::boolalpha << params.use_multilevel << std::endl;
  str << "  Number of Repititions:         " << params.num_repititions << std::endl;
  return str;
}

std::ostream & operator<< (std::ostream& str, const CoarseningParameters& params) {
  str << "Coarsening Parameters:" << std::endl;
  str << "  Algorithm:                     " << params.algorithm << std::endl;
  return str;
}

std::ostream & operator<< (std::ostream& str, const LabelPropagationRefinerParameters& params) {
  str << "\n  Label Propagation Refiner Parameters:" << std::endl;
  str << "    Maximum LP Iterations:       " << params.maximum_lp_iterations << std::endl;
  str << "    Act. Cliques After Round:    " << params.activate_all_cliques_after_rounds << std::endl;
  str << "    Random Shuffle each Round:   " << std::boolalpha << params.random_shuffle_each_round << std::endl;
  str << "    Node Ordering:               " << params.node_order << std::endl;
  return str;
}

std::ostream & operator<< (std::ostream& str, const RefinementParameters& params) {
  str << "Refinement Parameters:" << std::endl;
  if ( params.use_lp_refiner ) {
    str << params.lp;
  }
  return str;
}

std::ostream & operator<< (std::ostream& str, const Context& context) {
  str << "*******************************************************************************\n"
      << "*                         Cluster Editing Context                             *\n"
      << "*******************************************************************************\n"
      << context.general
      << "-------------------------------------------------------------------------------\n"
      << context.coarsening
      << "-------------------------------------------------------------------------------\n"
      << context.refinement
      << "-------------------------------------------------------------------------------\n";
  return str;
}

} // namespace cluster_editing