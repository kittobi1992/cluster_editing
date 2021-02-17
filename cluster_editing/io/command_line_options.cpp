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
            ("seed", po::value<int>(&context.general.seed)->value_name("<int>"),
             "Random Seed");
    return options;
  }

  po::options_description createGenericOptionsDescription(Context& context,
                                                          const int num_columns) {
    po::options_description options("General Options", num_columns);
    options.add_options()
            ("use-multilevel", po::value<bool>(&context.general.use_multilevel)->value_name("<bool>")->default_value(false),
             "If true, than multilevel paradigm is used.");
    return options;
  }

  po::options_description createCoarseningOptionsDescription(Context& context,
                                                             const int num_columns) {
    po::options_description options("Coarsening Options", num_columns);
    options.add_options()
            ("c-type",
             po::value<std::string>()->value_name("<string>")->notifier(
                     [&](const std::string& ctype) {
                       context.coarsening.algorithm = cluster_editing::coarseningAlgorithmFromString(ctype);
                     }),
             "Coarsening Algorithm:\n"
             " - do_nothing");
    return options;
  }

  po::options_description createRefinementOptionsDescription(Context& context,
                                                             const int num_columns) {
    po::options_description options("Refinement Options", num_columns);
    options.add_options()
            ("r-use-lp-refiner", po::value<bool>(&context.refinement.use_lp_refiner)->value_name("<bool>")->default_value(false),
             "If true, than label propagation is used to improve quality.")
            ("r-maximum-lp-iterations", po::value<int>(&context.refinement.lp.maximum_lp_iterations)->value_name("<int>"),
             "Maximum iterations made by the label propagation refiner");
    return options;
  }

  void processCommandLineInput(Context& context, int argc, char *argv[]) {
    const int num_columns = platform::getTerminalWidth();

    po::options_description general_options =
      createGeneralOptionsDescription(context, num_columns);
    po::options_description generic_options =
      createGenericOptionsDescription(context, num_columns);
    po::options_description coarsening_options =
      createCoarseningOptionsDescription(context, num_columns);
    po::options_description refinement_options =
      createRefinementOptionsDescription(context, num_columns);

    po::options_description cmd_line_options;
    cmd_line_options
      .add(general_options)
      .add(generic_options)
      .add(coarsening_options)
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
      .add(coarsening_options)
      .add(refinement_options);

    po::store(po::parse_config_file(file, ini_line_options, true), cmd_vm);
    po::notify(cmd_vm);

    context.general.output_file = context.general.graph_filename + ".sol";
  }

} // namespace cluster_editing