/*******************************************************************************
 * This file is part of KaPoCE.
 *
 * KaPoCE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaPoCE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaPoCE.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

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
    for ( const NodeID& v : graph.neighbors(u) ) {
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