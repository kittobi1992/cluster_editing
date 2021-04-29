#pragma once

#include "cluster_editing/definitions.h"
#include "cluster_editing/utils/common_operations.h"

namespace cluster_editing::metrics {

inline EdgeWeight edge_insertions(const Graph& graph) {
  EdgeWeight insertions = 0;
  utils::CommonOperations::instance(graph).computeClusterSizesAndInternalEdges(graph);
  const std::vector<NodeID>& cluster_sizes =
    utils::CommonOperations::instance(graph)._cluster_sizes;
  const std::vector<NodeID>& internal_edges =
    utils::CommonOperations::instance(graph)._internal_edges;
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
        ++deletions;
      }
    }
  }
  return deletions;
}

inline EdgeWeight edits(const Graph& graph) {
	return edge_deletions(graph) + edge_insertions(graph);
}

} // namespace metrics