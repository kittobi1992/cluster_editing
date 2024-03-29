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

#include "command_line_options.h"
#include "cluster_editing/macros.h"

#include <boost/program_options.hpp>
#include <sys/ioctl.h>

#include <iostream>
#include <fstream>
#include <limits>

namespace po = boost::program_options;

namespace cluster_editing {

  namespace platform {
    int getTerminalWidth() {
      int columns = 0;
      struct winsize w = { };
      ioctl(0, TIOCGWINSZ, &w);
      columns = w.ws_col;
      return columns;
    }

    int getProcessID() {
      return getpid();
    }
  }  // namespace platform


  po::options_description createGeneralOptionsDescription(Context& context,
                                                          const int num_columns) {
    po::options_description options("Required Options", num_columns);
    options.add_options()
            ("help", "show help message")
            ("graph,g",
             po::value<std::string>(&context.general.graph_filename)->value_name("<string>"),
             "Graph filename")
            ("config,c",
             po::value<std::string>(&context.general.config_file)->value_name("<string>"),
             "Context Config File (see config directory):\n"
             " - <path-to-custom-ini-file>")
            ("enable-logging", po::value<bool>(&context.general.enable_logging)->value_name("<bool>"),
             "If true, logging is enabled")
            ("verbose,v", po::value<bool>(&context.general.verbose_output)->value_name("<bool>"),
             "Verbose main partitioning output")
            ("print-result-line", po::value<bool>(&context.general.print_result_line)->value_name("<bool>"),
             "Prints RESULT-line containing stats, timings and metrics")
            ("write-to-file", po::value<bool>(&context.general.write_to_file)->value_name("<bool>"),
             "If true, then solution is written to a file")
            ("read-from-file", po::value<bool>(&context.general.read_from_file)->value_name("<bool>"),
             "If true, then solution is read from a file")
            ("csv", po::value<bool>(&context.general.print_csv)->value_name("<bool>"),
             "Prints CSV output")
            ("seed", po::value<int>(&context.general.seed)->value_name("<int>"),
             "Random Seed");
    return options;
  }

  po::options_description createGenericOptionsDescription(Context& context,
                                                          const int num_columns) {
    po::options_description options("General Options", num_columns);
    options.add_options()
            ("num-repititions", po::value<int>(&context.general.num_repititions)->value_name("<int>"),
             "Number of Repititions")
            ("num-fruitless-repititions", po::value<int>(&context.general.num_fruitless_repititions)->value_name("<int>"),
             "If for specified number of repititions no improvement is found, then the algorithm terminates.")
            ("time-limit", po::value<double>(&context.general.time_limit)->value_name("<double>"),
             "Time limit.");
    return options;
  }

  po::options_description createRefinementOptionsDescription(Context& context,
                                                             const int num_columns) {
    po::options_description options("Refinement Options", num_columns);
    options.add_options()
            ("r-use-evo", po::value<bool>(&context.refinement.use_evo)->value_name("<bool>")->default_value(false),
             "If true, then evolutionary algorithm is used.")
            ("r-use-localized-evo", po::value<bool>(&context.refinement.use_localized_evo)->value_name("<bool>")->default_value(false),
             "If true, then localized evolutionary algorithm is used.")
            ("r-use-lp-refiner", po::value<bool>(&context.refinement.use_lp_refiner)->value_name("<bool>")->default_value(false),
             "If true, then label propagation is used to improve quality.")
            ("r-use-boundary-fm-refiner", po::value<bool>(&context.refinement.use_boundary_fm_refiner)->value_name("<bool>")->default_value(false),
             "If true, then Boundary FM is used to improve quality.")
            ("r-use-localized-fm-refiner", po::value<bool>(&context.refinement.use_localized_fm_refiner)->value_name("<bool>")->default_value(false),
             "If true, then Localized FM is used to improve quality.")
             // Evolutionary Parameters
            ("r-evo-enable-detailed-output", po::value<bool>(&context.refinement.evo.enable_detailed_output)->value_name("<bool>"),
             "If true, then detailed output is shown in evolutionary algorithm")
            ("r-evo-time-limit", po::value<double>(&context.refinement.evo.time_limit)->value_name("<double>"),
             "Time limit for evolutionary algorithm.")
            ("r-evo-pool-size", po::value<int>(&context.refinement.evo.solution_pool_size)->value_name("<int>"),
             "Number of different solution in evolutionary algorithm")
            ("r-evo-steps", po::value<int>(&context.refinement.evo.evolutionary_steps)->value_name("<int>"),
             "Steps made by evolutionary algorithm")
            ("r-evo-initial-lp-iterations", po::value<int>(&context.refinement.evo.initial_lp_iterations)->value_name("<int>"),
             "Number of LP iterations for initial solution in evolutionary algorithm")
            ("r-evo-initial-node-swapper-iterations", po::value<int>(&context.refinement.evo.initial_node_swapper_iterations)->value_name("<int>"),
             "Number of node swapper iterations for initial solution in evolutionary algorithm")
            ("r-evo-intensivate-lp-iterations", po::value<int>(&context.refinement.evo.intensivate_lp_iterations)->value_name("<int>"),
             "Number of LP iterations for intensivate operation in evolutionary algorithm")
            ("r-evo-lp-iterations-after-mutate", po::value<int>(&context.refinement.evo.lp_iterations_after_mutate)->value_name("<int>"),
             "Number of LP iterations after mutate operation in evolutionary algorithm")
            ("r-evo-clique-remover-iterations", po::value<int>(&context.refinement.evo.clique_remover_iterations)->value_name("<int>"),
             "Number of clique remover iterations after mutate operation in evolutionary algorithm")
            ("r-evo-clique-splitter-iterations", po::value<int>(&context.refinement.evo.clique_splitter_iterations)->value_name("<int>"),
             "Number of clique splitter iterations after mutate operation in evolutionary algorithm")
            ("r-evo-node-swapper-iterations", po::value<int>(&context.refinement.evo.node_swapper_iterations)->value_name("<int>"),
             "Number of node swapper iterations after mutate operation in evolutionary algorithm")
            ("r-evo-node-swapper-max-cluster-size", po::value<int>(&context.refinement.evo.node_swapper_max_cluster_size)->value_name("<int>"),
             "Maximum cluster size consider in node swapper refiner")
            ("r-evo-lp-use-random-node-ordering", po::value<bool>(&context.refinement.evo.use_random_node_ordering)->value_name("<bool>"),
             "If true, then we choose a random node order in LP refiner")
            ("r-evo-enable-all-mutations-after-step", po::value<int>(&context.refinement.evo.enable_all_mutations_after_steps)->value_name("<int>"),
             "All mutations are automatically enabled after this number of steps")
            ("r-evo-enabled-mutations", po::value<std::string>(&context.refinement.evo.enabled_mutations)->value_name("<string>"),
             "String that indicates which mutations operators are enabled (e.g., 1101)")
            ("r-evo-large-clique-threshold", po::value<size_t>(&context.refinement.evo.large_clique_threshold)->value_name("<size_t>"),
             "All cliques greater than this threshold are considered as large")
            ("r-evo-min-clique-isolation-prob", po::value<float>(&context.refinement.evo.min_clique_isolate_prob)->value_name("<float>"),
             "Minimum probability for mutation operator that isolate cliques")
            ("r-evo-min-neighbor-clique-isolation-prob", po::value<float>(&context.refinement.evo.min_neighbor_clique_isolate_prob)->value_name("<float>"),
             "Minimum probability for mutation operator that isolate cliques with neighbors")
            ("r-evo-min-node-isolation-prob", po::value<float>(&context.refinement.evo.min_node_isolation_prob)->value_name("<float>"),
             "Minimum probability for mutation operator to isolate a node")
            ("r-evo-min-node-move-prob", po::value<float>(&context.refinement.evo.min_node_move_prob)->value_name("<float>"),
             "Minimum probability for mutation operator to move a node")
            ("r-evo-min-clique-split-mutation-prob", po::value<float>(&context.refinement.evo.min_clique_split_mutation_prob)->value_name("<float>"),
             "Minimum probability for clique split mutation operator")
            ("r-evo-min-test-mutation-prob", po::value<float>(&context.refinement.evo.min_test_mutation_prob)->value_name("<float>"),
             "Minimum probability for test mutation operator")
            ("r-evo-max-clique-isolation-prob", po::value<float>(&context.refinement.evo.max_clique_isolate_prob)->value_name("<float>"),
             "Maximum probability for mutation operator that isolate cliques")
            ("r-evo-max-neighbor-clique-isolation-prob", po::value<float>(&context.refinement.evo.max_neighbor_clique_isolate_prob)->value_name("<float>"),
             "Maximum probability for mutation operator that isolate cliques with neighbors")
            ("r-evo-max-node-isolation-prob", po::value<float>(&context.refinement.evo.max_node_isolation_prob)->value_name("<float>"),
             "Maximum probability for mutation operator to isolate a node")
            ("r-evo-max-node-move-prob", po::value<float>(&context.refinement.evo.max_node_move_prob)->value_name("<float>"),
             "Maximum probability for mutation operator to move a node")
            ("r-evo-max-clique-split-mutation-prob", po::value<float>(&context.refinement.evo.max_clique_split_mutation_prob)->value_name("<float>"),
             "Maximum probability for clique split mutation operator")
            ("r-evo-max-test-mutation-prob", po::value<float>(&context.refinement.evo.max_test_mutation_prob)->value_name("<float>"),
             "Maximum probability for test mutation operator")
            // Localized Evolutionary Parameters
            ("r-localized-evo-steps", po::value<size_t>(&context.refinement.localized_evo.steps)->value_name("<size_t>"),
             "Number of steps made by localized evolutionary algorithm")
            ("r-localized-evo-lp-iterations", po::value<int>(&context.refinement.localized_evo.max_lp_iterations)->value_name("<int>"),
             "Maximum LP iterations in localized evolutionary algorithm")
            ("r-localized-evo-min-mutations-nodes", po::value<int>(&context.refinement.localized_evo.min_mutations_nodes)->value_name("<int>"),
             "Minimum number of mutated nodes in localized evolutionary algorithm")
            ("r-localized-evo-max-mutations-nodes", po::value<int>(&context.refinement.localized_evo.max_mutations_nodes)->value_name("<int>"),
             "Maximum number of mutated nodes in localized evolutionary algorithm")
            ("r-localized-evo-choose-adjacent-mutation-node-prob", po::value<float>(&context.refinement.localized_evo.choose_adjacent_mutation_node_prob)->value_name("<float>"),
             "Probability that we choose an additional adjacent mutation node")
            ("r-localized-evo-max-distance-to-mutation-node", po::value<int>(&context.refinement.localized_evo.max_distance_to_mutation_node)->value_name("<int>"),
             "The maximum distance of a refinement node to a mutation node")
            ("r-localized-evo-degree-sampling-threshold", po::value<int>(&context.refinement.localized_evo.degree_sampling_threshold)->value_name("<int>"),
             "Some weird stuff")
            // Label Propagation Parameters
            ("r-maximum-lp-iterations", po::value<int>(&context.refinement.lp.maximum_lp_iterations)->value_name("<int>"),
             "Maximum iterations made by the label propagation refiner")
            ("r-random-shuffle-each-round", po::value<bool>(&context.refinement.lp.random_shuffle_each_round)->value_name("<bool>")->default_value(false),
             "If true, label propagation random shuffles all nodes after each iteration")
            ("r-lp-node-ordering",
             po::value<std::string>()->value_name("<string>")->notifier(
                     [&](const std::string& order) {
                       context.refinement.lp.node_order = cluster_editing::nodeOrderingFromString(order);
                     }),
             "Node Orderings:\n"
             " - none\n"
             " - random_shuffle\n"
             " - degree_increasing\n"
             " - degree_decreasing\n")
            ("r-lp-min-improvement", po::value<int>(&context.refinement.lp.min_improvement)->value_name("<int>"),
             "Minimal improvement per early exit window, if FM refiner is also activated")
            ("r-lp-rating-map-degree-threshold", po::value<NodeID>(&context.refinement.lp.rating_map_degree_threshold)->value_name("<uint32_t>"),
             "For all nodes with degree greater than this threshold, rating are aggregated in hash map."
             "For all others, incident edges are sorted and rating are aggregated on-the-fly.")
            ("r-lp-min-target-edit-distance", po::value<int>(&context.refinement.lp.min_target_edit_distance)->value_name("<int>"),
             "Intensivate on solution, if current edits are not that far away from target edits (important for evo)")
            ("r-lp-early-exit-window", po::value<size_t>(&context.refinement.lp.early_exit_window)->value_name("<size_t>"),
             "If label propagation improvement is less than min_improvement in the last early_exit_window rounds"
             "then label propagation terminates")
            // FM Parameters
            ("r-maximum-fm-iterations", po::value<int>(&context.refinement.fm.maximum_fm_iterations)->value_name("<int>"),
             "Maximum iterations made by the FM refiner")
            ("r-fm-min-improvement", po::value<int>(&context.refinement.fm.min_improvement)->value_name("<int>"),
             "Minimal improvement per early exit window")
            ("r-fm-early-exit-window", po::value<size_t>(&context.refinement.fm.early_exit_window)->value_name("<size_t>"),
             "If FM improvement is less than min_improvement in the last early_exit_window rounds"
             "then FM terminates")
            ("r-fm-fraction-fruitless-moves", po::value<double>(&context.refinement.fm.fraction_of_fruitless_moves)->value_name("<double>"),
             "FM terminates, if it not finds an improvement after fruitless_moves * |V| moves")
            ("r-fm-seed-nodes", po::value<size_t>(&context.refinement.fm.num_seed_nodes)->value_name("<size_t>"),
             "Number of seed nodes used by the localized FM Refiner");
    return options;
  }

  void processCommandLineInput(Context& context, int argc, char *argv[]) {
    const int num_columns = platform::getTerminalWidth();

    po::options_description general_options =
      createGeneralOptionsDescription(context, num_columns);
    po::options_description generic_options =
      createGenericOptionsDescription(context, num_columns);
    po::options_description refinement_options =
      createRefinementOptionsDescription(context, num_columns);

    po::options_description cmd_line_options;
    cmd_line_options
      .add(general_options)
      .add(generic_options)
      .add(refinement_options);

    po::variables_map cmd_vm;
    po::store(po::parse_command_line(argc, argv, cmd_line_options), cmd_vm);

    // placing vm.count("help") here prevents required attributes raising an
    // error if only help was supplied
    if (cmd_vm.count("help") != 0) {
      std::cout << cmd_line_options << std::endl;
      exit(0);
    }

    po::notify(cmd_vm);

    if ( context.general.config_file == "" ) {
      context.general.config_file = "strong.ini";
    }

    std::ifstream file(context.general.config_file.c_str());
    if (!file) {
      ERROR("Could not load config file at: " + context.general.config_file);
    }

    po::options_description ini_line_options;
    ini_line_options
      .add(generic_options)
      .add(refinement_options);

    po::store(po::parse_config_file(file, ini_line_options, true), cmd_vm);
    po::notify(cmd_vm);

    context.general.output_file = context.general.graph_filename + ".sol";
  }

} // namespace cluster_editing