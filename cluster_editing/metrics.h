#pragma once

#include "cluster_editing/definitions.h"

namespace cluster_editing::metrics {

inline EdgeWeight edge_insertions(const Graph& graph) {
  EdgeWeight insertions = 0;
  std::vector<NodeID> cluster_sizes(graph.numNodes(), 0);
  std::vector<NodeID> internal_edges(graph.numNodes(), 0);
  for ( const NodeID& u : graph.nodes() ) {
    const CliqueID u_id = graph.clique(u);
    cluster_sizes[u_id] += graph.nodeWeight(u);
    internal_edges[u_id] += graph.selfloopWeight(u);
    for ( const Neighbor& n : graph.neighbors(u) ) {
      const NodeID v = n.target;
      const CliqueID v_id = graph.clique(v);
      if ( u_id == v_id && u < v /* count each internal edge only once */ ) {
        internal_edges[u_id] += graph.edgeWeight(n.id);
      }
    }
  }

  for ( const NodeID& u : graph.nodes() ) {
    insertions += ( (cluster_sizes[u] * ( cluster_sizes[u] - 1 )) / 2 - internal_edges[u] );
  }
  return insertions;
}

inline EdgeWeight edge_deletions(const Graph& graph) {
  EdgeWeight deletions = 0;
  for ( const NodeID& u : graph.nodes() ) {
    for ( const Neighbor& n : graph.neighbors(u) ) {
      const NodeID v = n.target;
      // Only count deletions once
      if ( graph.clique(u) < graph.clique(v) ) {
        deletions += graph.edgeWeight(n.id);
      }
    }
  }
  return deletions;
}

inline EdgeWeight edits(const Graph& graph) {
	return edge_deletions(graph) + edge_insertions(graph);
}

} // namespace metrics