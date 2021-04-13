#include <iostream>
#include <signal.h>

#include "cluster_editing/definitions.h"
#include "cluster_editing/preprocessing.h"
#include "cluster_editing/flat.h"
#include "cluster_editing/multilevel.h"
#include "cluster_editing/datastructures/spin_lock.h"
#include "cluster_editing/utils/timer.h"
#include "cluster_editing/utils/randomize.h"
#include "cluster_editing/metrics.h"
#include "cluster_editing/io/command_line_options.h"
#include "cluster_editing/io/output.h"
#include "cluster_editing/io/graph_io.h"

using namespace cluster_editing;

SpinLock terminate_lock;
Context context;
Graph graph;
HighResClockTimepoint start, end;

void printResult(Graph& best) {
  // Print Stats
  std::chrono::duration<double> elapsed_seconds(end - start);
  io::printClusterEditingResults(best, context, elapsed_seconds);

  // Print RESULT line
  if ( context.general.print_result_line ) {
    io::printResultLine(best, context, elapsed_seconds);
  }

  if ( context.general.print_csv ) {
    std::cout << context.general.graph_filename.substr(context.general.graph_filename.find_last_of('/') + 1)
              << "," << metrics::edits(best) << "," << elapsed_seconds.count() << std::endl;
  }
}

void terminate(int) {
  if ( terminate_lock.tryLock() ) {
    end = std::chrono::high_resolution_clock::now();
    // Copy graph since refinement algorithm can checks
    // cliques while we try to output solution
    Graph cpy_graph = graph.copyBestSolution();
    printResult(cpy_graph);
    std::exit(0);
  }
}

int main(int argc, char* argv[]) {
  // Register signal handler
  struct sigaction action;
  memset(&action, 0, sizeof(struct sigaction));
  action.sa_handler = terminate;
  sigaction(SIGTERM, &action, NULL);

  processCommandLineInput(context, argc, argv);
  utils::Randomize::instance().setSeed(context.general.seed);

  utils::Timer::instance().start_timer("import_graph", "Import Graph");
  graph = io::readGraphFile(context.general.graph_filename);
  utils::Timer::instance().stop_timer("import_graph");
  context.computeParameters(graph.numNodes());

  if ( context.general.verbose_output ) {
    io::printBanner();
    // Print context description
    LOG << context;
  }
  io::printInputInfo(graph, context);

  // Preprocessing
  io::printPreprocessingBanner(context);
  start = std::chrono::high_resolution_clock::now();
  utils::Timer::instance().start_timer("preprocessing", "Preprocessing");
  Preprocessor preprocessor(graph, context);
  preprocessor.preprocess();
  utils::Timer::instance().stop_timer("preprocessing");

  // Multilevel Solver
  utils::Timer::instance().start_timer("solver", "Solver");
  int fruitless_repititions = 0;
  std::vector<CliqueID> best_cliques(graph.numNodes(), INVALID_CLIQUE);
  size_t best_objective = std::numeric_limits<size_t>::max();
  for ( int i = 0;
        i < context.general.num_repititions &&
        fruitless_repititions < context.general.num_fruitless_repititions ; ++i ) {
    graph.reset();
    if ( context.general.use_multilevel ) {
      multilevel::solve(graph, context);
    } else {
      flat::solve(graph, context);
    }

    // Check if solution is better than best solution found so far
    const size_t current_objective = metrics::edits(graph);
    if ( current_objective < best_objective ) {
      if ( context.general.verbose_output ) {
        LOG << GREEN << "Improved best solution from"
            << best_objective << "to" << current_objective << END;
      }
      graph.checkpoint(current_objective);
      best_objective = current_objective;
      fruitless_repititions = 0;
    } else {
      ++fruitless_repititions;
    }
  }

  utils::Timer::instance().stop_timer("solver");

  // Undo Preprocessing
  io::printUndoPreprocessingBanner(context);
  utils::Timer::instance().start_timer("undo_preprocessing", "Undo Preprocessing");
  preprocessor.undoPreprocessing();
  utils::Timer::instance().stop_timer("undo_preprocessing");
  end = std::chrono::high_resolution_clock::now();

  if ( terminate_lock.tryLock() ) {
    graph.applyBestCliques();
    printResult(graph);
  }

  return 0;
}
