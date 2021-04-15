#pragma once

#include "cluster_editing/refinement/i_refiner.h"
#include "cluster_editing/context/context.h"
#include "cluster_editing/datastructures/sparse_map.h"

#include "cluster_editing/datastructures/pq.h"
#include "cluster_editing/datastructures/fast_reset_flag_array.h"

namespace cluster_editing {

template<typename StoppingRule>
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
  explicit FMRefiner(const Graph& graph, const Context& context, const FMType type) :
      _context(context),
      _type(type),
      _nodes(),
      _clique_weight(graph.numNodes()),
      _empty_cliques(),
      edge_weight_to_clique(graph.numNodes()),
      n(graph.numNodes()),
      target_cliques(graph.numNodes()),
      current_cliques(graph.numNodes()),
      pq(graph.numNodes()),
      _moved_nodes(graph.numNodes()),
      _window_improvement(0),
      _round_improvements() { }

  size_t movedVertices() const {
    return _moved_vertices;
  }

private:

  void initializeImpl(Graph& graph) final;

  EdgeWeight refineImpl(Graph& graph) final ;

  EdgeWeight boundaryFM(Graph& graph, const EdgeWeight current_metric);

  EdgeWeight localizedFM(Graph& graph, const EdgeWeight current_metric);

  EdgeWeight localizedFMSearch(Graph& graph,
                               const EdgeWeight current_metric,
                               const std::vector<NodeID>& refinement_nodes);

  void moveVertex(Graph& graph, NodeID u, CliqueID to, bool manage_empty_cliques = true);

  void deltaGainUpdates(const Graph& graph, const NodeID u, const CliqueID from, const CliqueID to);

  void insertIntoPQ(const Graph& graph, const NodeID u);

  void clearPQ();

  void rollback(Graph& graph);

  Rating computeBestClique(const Graph& graph, const NodeID u);

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

  CliqueID getEmptyClique() {
    CliqueID empty_clique = INVALID_CLIQUE;
    while ( !_empty_cliques.empty() ) {
      CliqueID id = _empty_cliques.back();
      if ( _clique_weight[id] == 0 ) {
        empty_clique = id;
        break;
      }
      _empty_cliques.pop_back();
    }

    if ( empty_clique == INVALID_CLIQUE ) {
      for ( CliqueID id = 0; id < _clique_weight.size(); ++id ) {
        if ( _clique_weight[id] == 0 ) {
          _empty_cliques.push_back(id);
        }
      }
      ASSERT(!_empty_cliques.empty());
      empty_clique = _empty_cliques.back();
    }

    return empty_clique;
  }

  EdgeWeight insertions(NodeWeight target_clique_weight, EdgeWeight edge_weight_to_target_clique) const {
    return target_clique_weight - edge_weight_to_target_clique;
  }

  EdgeWeight deletions(EdgeWeight degree, EdgeWeight edge_weight_to_target_clique) const {
    return degree - edge_weight_to_target_clique;
  }

  EdgeWeight gain(CliqueID from, CliqueID target, EdgeWeight weight_to_current_clique, EdgeWeight weight_to_target_clique) const {
    const NodeWeight target_weight = target == ISOLATE_CLIQUE ? 0 : _clique_weight[target];
    return (target_weight - _clique_weight[from] + 1)
            + 2 * (weight_to_current_clique - weight_to_target_clique);
  }

  // ! Only for debug
  void checkPQGains(const Graph& graph);

  // ! Only for debug
  void checkCliqueWeights(const Graph& graph);

  const Context& _context;
  const FMType _type;
  size_t _moved_vertices;
  std::vector<NodeID> _nodes;
  std::vector<NodeWeight> _clique_weight;
  std::vector<CliqueID> _empty_cliques;

  ds::SparseMap<CliqueID, EdgeWeight> edge_weight_to_clique;
  std::vector<NodeData> n;
  std::vector<std::vector<NodeID>> target_cliques, current_cliques;
  mt_kahypar::ds::MinHeap<EdgeWeight, NodeID> pq;
  std::vector<Move> moves;
  ds::FastResetFlagArray<> _moved_nodes;
  EdgeWeight _window_improvement;
  std::vector<EdgeWeight> _round_improvements;
};
}  // namespace cluster_editing
