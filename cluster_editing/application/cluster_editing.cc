#include <iostream>

#include "cluster_editing/definitions.h"
#include "cluster_editing/preprocessing.h"
#include "cluster_editing/multilevel.h"
#include "cluster_editing/utils/timer.h"
#include "cluster_editing/io/command_line_options.h"
#include "cluster_editing/io/output.h"
#include "cluster_editing/io/graph_io.h"


int main(int argc, char* argv[]) {
  cluster_editing::Context context;
  cluster_editing::processCommandLineInput(context, argc, argv);

  if ( context.general.verbose_output ) {
    cluster_editing::io::printBanner();
    // Print context description
    LOG << context;
  }

  cluster_editing::utils::Timer::instance().start_timer("import_graph", "Import Graph");
  cluster_editing::Graph graph = cluster_editing::io::readGraphFile(context.general.graph_filename);
  cluster_editing::utils::Timer::instance().stop_timer("import_graph");
  cluster_editing::io::printInputInfo(graph, context);

  // Preprocessing
  cluster_editing::io::printPreprocessingBanner(context);
  cluster_editing::HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  cluster_editing::utils::Timer::instance().start_timer("preprocessing", "Preprocessing");
  cluster_editing::Preprocessor preprocessor(graph, context);
  preprocessor.preprocess();
  cluster_editing::utils::Timer::instance().stop_timer("preprocessing");

  // Multilevel Solver
  cluster_editing::utils::Timer::instance().start_timer("multilevel_solver", "Multilevel Solver");
  cluster_editing::multilevel::solve(graph, context);
  cluster_editing::utils::Timer::instance().stop_timer("multilevel_solver");

  // Undo Preprocessing
  cluster_editing::io::printUndoPreprocessingBanner(context);
  cluster_editing::utils::Timer::instance().start_timer("undo_preprocessing", "Undo Preprocessing");
  preprocessor.undoPreprocessing();
  cluster_editing::utils::Timer::instance().stop_timer("undo_preprocessing");
  cluster_editing::HighResClockTimepoint end = std::chrono::high_resolution_clock::now();

  // Print Stats
  std::chrono::duration<double> elapsed_seconds(end - start);
  cluster_editing::io::printClusterEditingResults(graph, context, elapsed_seconds);

  return 0;
}
