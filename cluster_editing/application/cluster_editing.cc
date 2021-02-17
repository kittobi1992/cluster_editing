#include <iostream>

#include "cluster_editing/definitions.h"
#include "cluster_editing/preprocessing.h"
#include "cluster_editing/flat.h"
#include "cluster_editing/multilevel.h"
#include "cluster_editing/utils/timer.h"
#include "cluster_editing/utils/randomize.h"
#include "cluster_editing/metrics.h"
#include "cluster_editing/io/command_line_options.h"
#include "cluster_editing/io/output.h"
#include "cluster_editing/io/graph_io.h"

using namespace cluster_editing;

int main(int argc, char* argv[]) {
  Context context;
  processCommandLineInput(context, argc, argv);
  utils::Randomize::instance().setSeed(context.general.seed);

  if ( context.general.verbose_output ) {
    io::printBanner();
    // Print context description
    LOG << context;
  }

  utils::Timer::instance().start_timer("import_graph", "Import Graph");
  Graph graph = io::readGraphFile(context.general.graph_filename);
  utils::Timer::instance().stop_timer("import_graph");
  io::printInputInfo(graph, context);

  // Preprocessing
  io::printPreprocessingBanner(context);
  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  utils::Timer::instance().start_timer("preprocessing", "Preprocessing");
  Preprocessor preprocessor(graph, context);
  preprocessor.preprocess();
  utils::Timer::instance().stop_timer("preprocessing");

  // Multilevel Solver
  utils::Timer::instance().start_timer("solver", "Solver");
  std::vector<CliqueID> best_cliques(
    graph.numNodes(), INVALID_CLIQUE);
  size_t best_objective = std::numeric_limits<size_t>::max();
  for ( int i = 0; i < context.general.num_repititions; ++i ) {
    graph.reset();
    if ( context.general.use_multilevel ) {
      multilevel::solve(graph, context);
    } else {
      flat::solve(graph, context);
    }

    // Check if solution is better than best solution found so far
    const size_t current_objective =
      metrics::edge_deletions(graph) +
      metrics::edge_insertions(graph);
    if ( current_objective < best_objective ) {
      for ( const NodeID& u : graph.nodes() ) {
        best_cliques[u] = graph.clique(u);
      }
      best_objective = current_objective;
    }
  }

  // Apply best found solution
  for ( const NodeID& u : graph.nodes() ) {
    graph.setClique(u, best_cliques[u]);
  }

  utils::Timer::instance().stop_timer("solver");

  // Undo Preprocessing
  io::printUndoPreprocessingBanner(context);
  utils::Timer::instance().start_timer("undo_preprocessing", "Undo Preprocessing");
  preprocessor.undoPreprocessing();
  utils::Timer::instance().stop_timer("undo_preprocessing");
  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();

  // Print Stats
  std::chrono::duration<double> elapsed_seconds(end - start);
  io::printClusterEditingResults(graph, context, elapsed_seconds);

  // Print RESULT line
  if ( context.general.print_result_line ) {
    io::printResultLine(graph, context, elapsed_seconds);
  }

  return 0;
}
