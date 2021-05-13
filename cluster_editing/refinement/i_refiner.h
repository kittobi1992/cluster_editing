#pragma once

#include "cluster_editing/definitions.h"

namespace cluster_editing {

class IRefiner {

 public:
  IRefiner(const IRefiner&) = delete;
  IRefiner(IRefiner&&) = delete;
  IRefiner & operator= (const IRefiner &) = delete;
  IRefiner & operator= (IRefiner &&) = delete;

  virtual ~IRefiner() = default;

  void initialize(Graph& graph) {
    initializeImpl(graph);
  }

  EdgeWeight refine(Graph& graph,
                    const EdgeWeight current_edits,
                    const EdgeWeight target_edits = 0) {
    return refineImpl(graph, current_edits, target_edits);
  }

 protected:
  IRefiner() = default;

 private:
  virtual void initializeImpl(Graph& graph) = 0;
  virtual EdgeWeight refineImpl(Graph& graph,
                                const EdgeWeight current_edits,
                                const EdgeWeight target_edits) = 0;
};

} // namespace cluster_editing