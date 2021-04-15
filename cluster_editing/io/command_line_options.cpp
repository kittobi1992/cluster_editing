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
             po::value<std::string>(&context.general.graph_filename)->value_name("<string>")->required(),
             "Graph filename")
            ("config,c",
             po::value<std::string>(&context.general.config_file)->value_name("<string>"),
             "Context Config File (see config directory):\n"
             " - <path-to-custom-ini-file>")
            ("verbose,v", po::value<bool>(&context.general.verbose_output)->value_name("<bool>")->default_value(true),
             "Verbose main partitioning output")
            ("print-result-line", po::value<bool>(&context.general.print_result_line)->value_name("<bool>")->default_value(true),
             "Prints RESULT-line containing stats, timings and metrics")
            ("csv", po::value<bool>(&context.general.print_csv)->value_name("<bool>")->default_value(false),
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
             "If for specified number of repititions no improvement is found, then the algorithm terminates.");
    return options;
  }

  po::options_description createRefinementOptionsDescription(Context& context,
                                                             const int num_columns) {
    po::options_description options("Refinement Options", num_columns);
    options.add_options()
            ("r-use-lp-refiner", po::value<bool>(&context.refinement.use_lp_refiner)->value_name("<bool>")->default_value(false),
             "If true, then label propagation is used to improve quality.")
            ("r-use-boundary-fm-refiner", po::value<bool>(&context.refinement.use_boundary_fm_refiner)->value_name("<bool>")->default_value(false),
             "If true, then Boundary FM is used to improve quality.")
            ("r-use-localized-fm-refiner", po::value<bool>(&context.refinement.use_localized_fm_refiner)->value_name("<bool>")->default_value(false),
             "If true, then Localized FM is used to improve quality.")
            ("r-maximum-lp-repititions", po::value<int>(&context.refinement.lp.maximum_lp_repititions)->value_name("<int>"),
             "Maximum number of restarts of LP refiner per run")
            ("r-maximum-lp-iterations", po::value<int>(&context.refinement.lp.maximum_lp_iterations)->value_name("<int>"),
             "Maximum iterations made by the label propagation refiner")
            ("r-maximum-fm-iterations", po::value<int>(&context.refinement.fm.maximum_fm_iterations)->value_name("<int>"),
             "Maximum iterations made by the FM refiner")
            ("r-activate-all-cliques-after-rounds", po::value<int>(&context.refinement.lp.activate_all_cliques_after_rounds)->value_name("<int>"),
             "Each #activate_all_cliques_after_rounds iterations, label propagation refiner reactivates all nodes again.")
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
            ("r-lp-early-exit-window", po::value<size_t>(&context.refinement.lp.early_exit_window)->value_name("<size_t>"),
             "If label propagation improvement is less than min_improvement in the last early_exit_window rounds"
             "then label propagation terminates")
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
    if (cmd_vm.count("help") != 0 || argc == 1) {
      std::cout << cmd_line_options << std::endl;
      exit(0);
    }

    po::notify(cmd_vm);

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