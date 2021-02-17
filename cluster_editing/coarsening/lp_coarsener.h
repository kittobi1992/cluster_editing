#pragma once

#include "cluster_editing/definitions.h"
#include "cluster_editing/context/context.h"
#include "cluster_editing/datastructures/sparse_map.h"
#include "cluster_editing/coarsening/i_coarsener.h"
#include "cluster_editing/refinement/lp_refiner.h"

namespace cluster_editing {

class LabelPropagationCoarsener final : public ICoarsener {

  static constexpr bool debug = false;

  struct Rating {
    CliqueID clique;
    EdgeWeight rating;
    EdgeWeight delta;
  };

  using Hierarchy = std::pair<Graph, std::vector<NodeID>>;

 public:
  LabelPropagationCoarsener(const LabelPropagationCoarsener&) = delete;
  LabelPropagationCoarsener(LabelPropagationCoarsener&&) = delete;
  LabelPropagationCoarsener & operator= (const LabelPropagationCoarsener &) = delete;
  LabelPropagationCoarsener & operator= (LabelPropagationCoarsener &&) = delete;

  explicit LabelPropagationCoarsener(Graph& graph, const Context& context) :
    _graph(graph),
    _context(context),
    _lp_refiner(graph, context),
    _hierarchies() { }


 private:
  void coarsenImpl() override final;

  void uncoarsenImpl(std::unique_ptr<IRefiner>& refiner) override final;

  Graph& currentGraph() {
    if ( _hierarchies.empty() ) {
      return _graph;
    } else {
      return _hierarchies.back().first;
    }
  }

  Graph& _graph;
  const Context& _context;
  LabelPropagationRefiner _lp_refiner;
  std::vector<Hierarchy> _hierarchies;
};
}  // namespace cluster_editing
