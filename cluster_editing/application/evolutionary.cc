#include "cluster_editing/io/command_line_options.h"
#include "cluster_editing/io/graph_io.h"
#include "cluster_editing/utils/randomize.h"
#include "cluster_editing/io/output.h"

#include "cluster_editing/evolutionary/evolutionary.h"

using namespace cluster_editing;

int main(int argc, char* argv[]) {
  Context context;
  processCommandLineInput(context, argc, argv);
  utils::Randomize::instance().setSeed(context.general.seed);
  Graph graph = io::readGraphFile(context.general.graph_filename);
  auto start = std::chrono::steady_clock::now();
  double time_limit = 10 /* minutes */ * 60 /* seconds */;
  time_limit = 15; // seconds

  // TODO lower bound?


  evolutionary::EvolutionaryAlgorithm evo(graph, context);
  evo.generate_initial_population();
  std::chrono::duration<double> elapsed_time { std::chrono::steady_clock::now() - start };
  size_t max_evo_steps = 1000;
  size_t steps = 0;
  while (steps++ < max_evo_steps && elapsed_time.count() < time_limit) {
    evo.evolution_step();
    elapsed_time = std::chrono::steady_clock::now() - start;
  }

  // Print Stats
  std::chrono::duration<double> elapsed_seconds(std::chrono::steady_clock::now() - start);

  if ( context.general.verbose_output || context.general.print_result_line ) {
    // write best solution to graph object
    evo.apply_best_solution();
  }

  if ( context.general.verbose_output ) {
    io::printClusterEditingResults(graph, context, elapsed_seconds);
  }

  if ( context.general.print_result_line ) {
    io::printResultLine(graph, context, elapsed_seconds);
  }

  if ( context.general.print_csv ) {
    std::cout << context.general.graph_filename.substr(context.general.graph_filename.find_last_of('/') + 1)
              << "," << evo.best().cost << "," << elapsed_seconds.count() << std::endl;
  }

  if ( context.general.print_edits ) {
    evo.print_best_solution();
  }
  return 0;
}
