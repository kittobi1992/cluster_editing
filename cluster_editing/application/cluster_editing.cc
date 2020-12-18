#include <iostream>

#include "cluster_editing/io/command_line_options.h"
#include "cluster_editing/definitions.h"
#include "cluster_editing/preprocessing.h"
#include "cluster_editing/multilevel.h"
#include "cluster_editing/utils/timer.h"

int main(int argc, char* argv[]) {
  cluster_editing::Context context;
  cluster_editing::processCommandLineInput(context, argc, argv);

  if ( context.general.verbose_output ) {
    // Print context description
    LOG << context;
  }

  cluster_editing::Graph graph;

  cluster_editing::Preprocessor preprocessor(graph, context);
  cluster_editing::utils::Timer::instance().start_timer("preprocessing", "Preprocessing");
  preprocessor.preprocess();
  cluster_editing::utils::Timer::instance().stop_timer("preprocessing");

  cluster_editing::utils::Timer::instance().start_timer("multilevel_solver", "Multilevel Solver");
  cluster_editing::multilevel::solve(graph, context);
  cluster_editing::utils::Timer::instance().stop_timer("multilevel_solver");

  cluster_editing::utils::Timer::instance().start_timer("undo_preprocessing", "Undo Preprocessing");
  preprocessor.undoPreprocessing();
  cluster_editing::utils::Timer::instance().stop_timer("undo_preprocessing");

  if ( context.general.verbose_output ) {
    // Print timings
    LOG << cluster_editing::utils::Timer::instance(true);
  }

  return 0;
}
