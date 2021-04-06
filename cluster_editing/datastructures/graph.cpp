#include "graph.h"

#include <vector>

#include "cluster_editing/macros.h"
#include "cluster_editing/datastructures/sparse_map.h"
#include "cluster_editing/utils/timer.h"

namespace cluster_editing {
namespace ds {

  std::pair<Graph, std::vector<NodeID>> Graph::contract() const {
    Graph c_graph;
    std::vector<NodeID> clique_to_coarse_node(numNodes(), INVALID_NODE);

    // Create clique to coarse node mapping
    utils::Timer::instance().start_timer("create_mapping", "Create Mapping");
    std::vector<std::vector<NodeID>> nodes_in_clique(numNodes());
    NodeID num_nodes = 0;
    for ( const NodeID& u : nodes() ) {
      const CliqueID u_clique = clique(u);
      NodeID coarse_node = clique_to_coarse_node[u_clique];
      if ( coarse_node == INVALID_NODE ) {
        coarse_node = num_nodes++;
      }
      ASSERT(coarse_node != INVALID_NODE);
      clique_to_coarse_node[u_clique] = coarse_node;
      nodes_in_clique[coarse_node].push_back(u);
    }
    c_graph._num_nodes = num_nodes;
    c_graph._total_weight = _total_weight;
    c_graph._nodes.resize(num_nodes + 1);
    nodes_in_clique.resize(num_nodes);
    utils::Timer::instance().stop_timer("create_mapping");

    // Create edges of coarse graph
    utils::Timer::instance().start_timer("contract_edges", "Contract Edges");
    SparseMap<NodeID, EdgeWeight> aggregated_edge_weight(num_nodes);
    for ( CliqueID c = 0; c < num_nodes; ++c ) {
      aggregated_edge_weight.clear();
      c_graph.node(c).setClique(c);
      c_graph.node(c).setFirstEntry(c_graph._edges.size());
      NodeWeight contracted_weight = 0;
      EdgeWeight selfloop_weight = 0;
      for ( const NodeID& u : nodes_in_clique[c] ) {
        contracted_weight += nodeWeight(u);
        selfloop_weight += selfloopWeight(u);
        for ( const Neighbor& v : neighbors(u) ) {
          const CliqueID c_v = clique_to_coarse_node[clique(v.target)];
          if ( c != c_v ) {
            aggregated_edge_weight[c_v] += edgeWeight(v.id);
          } else if ( u < v.target ) {
            // Only aggregate weight of selfloops once
            selfloop_weight += edgeWeight(v.id);
          }
        }
      }
      c_graph.node(c).setWeight(contracted_weight);
      c_graph.node(c).setSelfloopWeight(selfloop_weight);

      EdgeWeight weighted_degree = 0;
      for ( const auto& entry : aggregated_edge_weight ) {
        const NodeID v = entry.key;
        const EdgeWeight weight = entry.value;
        weighted_degree += weight;
        ASSERT(c != v, "Selfloop contained in aggregated weights");
        c_graph._edges.emplace_back();
        c_graph._edges.back().setSource(c);
        c_graph._edges.back().setTarget(v);
        c_graph._edges.back().setWeight(weight);
      }
      c_graph.node(c).setWeightedDegree(weighted_degree);
    }
    c_graph.node(num_nodes).setFirstEntry(c_graph._edges.size());
    c_graph._num_edges = static_cast<EdgeID>(c_graph._edges.size());
    utils::Timer::instance().stop_timer("contract_edges");

    if ( c_graph._num_edges == 0 ) {
      // Add dummy edge
      c_graph._edges.emplace_back();
    }

    return std::make_pair(c_graph, clique_to_coarse_node);
  }

  void Graph::reset() {
    for ( const NodeID& u : nodes() ) {
      _nodes[u].setClique(u);
    }
  }
}
} // namespace cluster_editing::ds