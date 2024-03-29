/*******************************************************************************
 * This file is part of KaPoCE.
 *
 * KaPoCE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaPoCE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaPoCE.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include <iostream>

#include "context.h"

namespace cluster_editing {

std::ostream & operator<< (std::ostream& str, const GeneralParameters& params) {
  str << "General Parameters:" << std::endl;
  str << "  Graph:                         " << params.graph_filename << std::endl;
  str << "  Config File:                   " << params.config_file << std::endl;
  str << "  Output File:                   " << params.output_file << std::endl;
  str << "  Write to File:                 " << std::boolalpha << params.write_to_file << std::endl;
  str << "  Read from File:                " << std::boolalpha << params.read_from_file << std::endl;
  str << "  Verbose Output:                " << std::boolalpha << params.verbose_output << std::endl;
  str << "  Seed:                          " << params.seed << std::endl;
  str << "  Number of Repititions:         " << params.num_repititions << std::endl;
  str << "  Number Fruitless Repititions:  " << params.num_fruitless_repititions << std::endl;
  str << "  Time Limit:                    " << params.time_limit << std::endl;
  return str;
}


std::ostream & operator<< (std::ostream& str, const EvolutionaryParameters& params) {
  str << "\n  Evolutionary Parameters:" << std::endl;
  str << "    Enable Detailed Output:      " << std::boolalpha << params.enable_detailed_output << std::endl;
  str << "    Evolutionary Time Limit:     " << params.time_limit << std::endl;
  str << "    Run until Time Limit:        " << std::boolalpha << params.run_until_time_limit << std::endl;
  str << "    Evolutionary Steps:          " << params.evolutionary_steps << std::endl;
  str << "    Solution Pool Size:          " << params.solution_pool_size << std::endl;
  str << "    Initial LP Iterations:       " << params.initial_lp_iterations << std::endl;
  str << "    Initial Node Swapper Iter.:  " << params.initial_node_swapper_iterations << std::endl;
  str << "    Intensivate LP Iterations:   " << params.intensivate_lp_iterations << std::endl;
  str << "    LP Iterations After Mutate:  " << params.lp_iterations_after_mutate << std::endl;
  str << "    Clique Remover Iterations:   " << params.clique_remover_iterations << std::endl;
  str << "    Clique Splitter Iterations:  " << params.clique_splitter_iterations << std::endl;
  str << "    Node Swapper Iterations:     " << params.node_swapper_iterations << std::endl;
  str << "    Node Swapper Max Cl. Size:   " << params.node_swapper_max_cluster_size << std::endl;
  str << "    Use Random Node Ordering:    " << std::boolalpha << params.use_random_node_ordering << std::endl;
  str << "    Enable All Mut. After Steps: " << params.enable_all_mutations_after_steps << std::endl;
  str << "    Enabled Mutations:           " << params.enabled_mutations << std::endl;
  str << "    Large Clique Threshold:      " << params.large_clique_threshold << std::endl;
  str << "    Min Clique Isolation Prob.:  " << params.min_clique_isolate_prob << std::endl;
  str << "    Max Clique Isolation Prob.:  " << params.max_clique_isolate_prob << std::endl;
  str << "    Min NeighClique Isolation P: " << params.min_neighbor_clique_isolate_prob << std::endl;
  str << "    Max NeighClique Isolation P: " << params.max_neighbor_clique_isolate_prob << std::endl;
  str << "    Min Node Isolation Prob.:    " << params.min_node_isolation_prob << std::endl;
  str << "    Max Node Isolation Prob.:    " << params.max_node_isolation_prob << std::endl;
  str << "    Min Node Move Prob.:         " << params.min_node_move_prob << std::endl;
  str << "    Max Node Move Prob.:         " << params.max_node_move_prob << std::endl;
  str << "    Min Clique Split Prob.:      " << params.min_clique_split_mutation_prob << std::endl;
  str << "    Max Clique Split Prob.:      " << params.max_clique_split_mutation_prob << std::endl;
  str << "    Min Test Mutation Prob.:     " << params.min_test_mutation_prob << std::endl;
  str << "    Max Test Mutation Prob.:     " << params.max_test_mutation_prob << std::endl;
  return str;
}

std::ostream & operator<< (std::ostream& str, const LocalizedEvolutionaryParameters& params) {
  str << "\n  Localized Evolutionary Parameters:" << std::endl;
  str << "    Localized Evo Steps:         " << params.steps << std::endl;
  str << "    Run until Time Limit:        " << std::boolalpha << params.run_until_time_limit << std::endl;
  str << "    Maximum LP Iterations:       " << params.max_lp_iterations << std::endl;
  str << "    Min. Mutation Nodes:         " << params.min_mutations_nodes << std::endl;
  str << "    Max. Mutation Nodes:         " << params.max_mutations_nodes << std::endl;
  str << "    Choose Adj. Mut. Node Prob.: " << params.choose_adjacent_mutation_node_prob << std::endl;
  str << "    Max. Distance to Mut. Node:  " << params.max_distance_to_mutation_node << std::endl;
  str << "    Degree Sampling Threshold:   " << params.degree_sampling_threshold << std::endl;
  return str;
}

std::ostream & operator<< (std::ostream& str, const LabelPropagationRefinerParameters& params) {
  str << "\n  Label Propagation Refiner Parameters:" << std::endl;
  str << "    Maximum LP Iterations:       " << params.maximum_lp_iterations << std::endl;
  str << "    Random Shuffle each Round:   " << std::boolalpha << params.random_shuffle_each_round << std::endl;
  str << "    Node Ordering:               " << params.node_order << std::endl;
  str << "    Rating Map Degree Threshold: " << params.rating_map_degree_threshold << std::endl;
  str << "    Min. Target Edit Distance:   " << params.min_target_edit_distance << std::endl;
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
  if ( params.use_localized_evo ) {
    str << params.localized_evo;
  }
  if ( params.use_evo || params.use_lp_refiner ) {
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
} // namespace cluster_editing