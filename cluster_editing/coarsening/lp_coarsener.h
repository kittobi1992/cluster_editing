#pragma once

#include "cluster_editing/definitions.h"
#include "cluster_editing/context/context.h"
#include "cluster_editing/datastructures/sparse_map.h"
#include "cluster_editing/coarsening/i_coarsener.h"

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
    _clique_weight(graph.numNodes()),
    _empty_cliques(),
    _rating(graph.numNodes()),
    _hierarchies() { }


 private:
  void coarsenImpl() override final;

  void uncoarsenImpl(std::unique_ptr<IRefiner>& refiner) override final;

  void moveVertex(Graph& graph, const NodeID u, const CliqueID to);

  Rating computeBestClique(Graph& graph, const NodeID u);

  Graph& currentGraph() {
    if ( _hierarchies.empty() ) {
      return _graph;
    } else {
      return _hierarchies.back().first;
    }
  }

  Graph& _graph;
  const Context& _context;
  std::vector<NodeWeight> _clique_weight;
  std::vector<CliqueID> _empty_cliques;
  ds::SparseMap<CliqueID, EdgeWeight> _rating;
  std::vector<Hierarchy> _hierarchies;
};
}  // namespace cluster_editing
