#include <iostream>

#include "cluster_editing/io/command_line_options.h"
#include "cluster_editing/definitions.h"
#include "cluster_editing/preprocessing.h"
#include "cluster_editing/multilevel.h"
#include "cluster_editing/macros.h"

int main(int argc, char* argv[]) {
  cluster_editing::Context context;
  cluster_editing::processCommandLineInput(context, argc, argv);

  std::cout << context << std::endl;

  cluster_editing::Graph graph;

  cluster_editing::Preprocessor preprocessor(graph, context);
  preprocessor.preprocess();

  cluster_editing::multilevel::solve(graph, context);

  preprocessor.undoPreprocessing();

  return 0;
}
