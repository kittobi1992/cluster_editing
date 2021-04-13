#pragma once

#include "cluster_editing/refinement/i_refiner.h"
#include "cluster_editing/context/context.h"
#include "cluster_editing/datastructures/sparse_map.h"

#include "cluster_editing/datastructures/pq.h"

namespace cluster_editing {

class FMRefiner final : public IRefiner {
private:
  static constexpr bool debug = false;

  struct Rating {
    CliqueID clique;
    EdgeWeight rating;
    EdgeWeight delta;
    EdgeWeight weight_to_clique;
  };

  struct NodeData {
    CliqueID desired_target = INVALID_CLIQUE;
    EdgeWeight weight_to_target_clique, weight_to_current_clique = 0;
    NodeID index_in_target_clique = INVALID_NODE, index_in_current_clique = INVALID_NODE;
    uint32_t num_skips = 0,  last_move_round = 0;
  };

  struct Move {
    NodeID node;
    CliqueID from, to;
  };

  CliqueID ISOLATE_CLIQUE = INVALID_CLIQUE - 1;

public:
  explicit FMRefiner(const Graph& graph, const Context& context) :
      _context(context),
      _nodes(),
      _clique_weight(graph.numNodes()),
      _empty_cliques(),
      edge_weight_to_clique(graph.numNodes()),
      n(graph.numNodes()),
      target_cliques(graph.numNodes()),
      current_cliques(graph.numNodes()),
      pq(graph.numNodes())
      { }

  size_t movedVertices() const {
    return _moved_vertices;
  }

private:

  void initializeImpl(Graph& graph) final;

  bool refineImpl(Graph& graph) final ;

  void moveVertex(Graph& graph, NodeID u, CliqueID to, bool manage_empty_cliques = true);

  Rating computeBestClique(Graph& graph, const NodeID u);

  void updateTargetClique(NodeID u, Rating& r) {
    removeFromTargetClique(u);
    insertIntoTargetClique(u, r);
  }

  void removeFromTargetClique(NodeID u) {
    CliqueID old_target = n[u].desired_target;
    if (old_target != ISOLATE_CLIQUE) {
      target_cliques[old_target][n[u].index_in_target_clique] = target_cliques[old_target].back();
      n[target_cliques[old_target].back()].index_in_target_clique = n[u].index_in_target_clique;
      target_cliques[old_target].pop_back();
      n[u].index_in_target_clique = std::numeric_limits<uint32_t>::max();
    }
  }

  void insertIntoTargetClique(NodeID u, Rating& r) {
    if (r.clique == ISOLATE_CLIQUE) assert(r.weight_to_clique == 0);
    CliqueID target = r.clique;
    n[u].desired_target = target;
    n[u].weight_to_target_clique = r.weight_to_clique;
    if (target != ISOLATE_CLIQUE) {
      n[u].index_in_target_clique = target_cliques[target].size();
      target_cliques[target].push_back(u);
    }
  }

  void removeFromCurrentClique(NodeID u, CliqueID current) {
    current_cliques[current][n[u].index_in_current_clique] = current_cliques[current].back();
    n[current_cliques[current].back()].index_in_current_clique = n[u].index_in_current_clique;
    current_cliques[current].pop_back();
    n[u].index_in_current_clique = std::numeric_limits<uint32_t>::max();
  }

  void insertIntoCurrentClique(NodeID u, CliqueID to, EdgeWeight weight_to_clique) {
    n[u].weight_to_current_clique = weight_to_clique;
    n[u].index_in_current_clique = current_cliques[to].size();
    current_cliques[to].push_back(u);
  }


  EdgeWeight insertions(NodeWeight u_weight, NodeWeight target_clique_weight, EdgeWeight edge_weight_to_target_clique) const {
    return u_weight * target_clique_weight - edge_weight_to_target_clique;
  }

  EdgeWeight deletions(EdgeWeight weighted_degree, EdgeWeight edge_weight_to_target_clique) const {
    return weighted_degree - edge_weight_to_target_clique;
  }

  EdgeWeight gain(const Graph& graph, NodeID u, CliqueID from, CliqueID target, EdgeWeight weight_to_current_clique, EdgeWeight weight_to_target_clique) const {
    const NodeWeight target_weight = target == ISOLATE_CLIQUE ? 0 : _clique_weight[target];
    return (target_weight - _clique_weight[from] + graph.nodeWeight(u)) * graph.nodeWeight(u)
            + 2 * (weight_to_current_clique - weight_to_target_clique);
  }

  // ! Only for debug
  void checkPQGains(const Graph& graph);

  // ! Only for debug
  void checkCliqueWeights(const Graph& graph);

  const Context& _context;
  size_t _moved_vertices;
  std::vector<NodeID> _nodes;
  std::vector<NodeWeight> _clique_weight;
  std::vector<CliqueID> _empty_cliques;
  ds::SparseMap<CliqueID, EdgeWeight> edge_weight_to_clique;

  std::vector<NodeData> n;
  std::vector<std::vector<NodeID>> target_cliques, current_cliques;
  mt_kahypar::ds::MinHeap<EdgeWeight, NodeID> pq;
  std::vector<Move> moves;

};
}  // namespace cluster_editing
