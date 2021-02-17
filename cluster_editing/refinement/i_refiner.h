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

  bool refine(Graph& graph) {
    return refineImpl(graph);
  }

 protected:
  IRefiner() = default;

 private:
  virtual void initializeImpl(Graph& graph) = 0;
  virtual bool refineImpl(Graph& graph) = 0;
};

} // namespace cluster_editing