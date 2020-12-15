#pragma once

#include "cluster_editing/context/context.h"
#include "cluster_editing/definitions.h"

namespace cluster_editing {

class Preprocessor {

 public:
  Preprocessor(Graph& graph, const Context& context) :
    _graph(graph),
    _context(context) { }

  void preprocess();
  void undoPreprocessing();

 private:
  Graph& _graph;
  const Context& _context;
};

} // namespace cluster_editing