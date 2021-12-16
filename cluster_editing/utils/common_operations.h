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

#include <algorithm>
#include <cstring>
#include <limits>
#include <vector>

#include "cluster_editing/definitions.h"
#include "cluster_editing/datastructures/sparse_map.h"

namespace cluster_editing {
namespace utils {
class CommonOperations {
 public:
  CommonOperations(const CommonOperations&) = delete;
  CommonOperations(CommonOperations&&) = delete;
  CommonOperations& operator= (const CommonOperations&) = delete;
  CommonOperations& operator= (CommonOperations&&) = delete;

  static bool dirty;

  static CommonOperations & instance(const Graph& graph) {
    static auto* instance = new CommonOperations(graph);
    if(dirty) {
        delete instance;
        instance = new CommonOperations(graph);
        dirty = false;
    }
    return *instance;
  }

  void computeClusterSizes(const Graph& graph) {
    reset(_cluster_sizes);
    for ( const NodeID& u : graph.nodes() ) {
      ++_cluster_sizes[graph.clique(u)];
    }
  }

  void computeClusterSizesAndInternalEdges(const Graph& graph) {
    reset(_cluster_sizes);
    reset(_internal_edges);
    for ( const NodeID& u : graph.nodes() ) {
      const CliqueID u_id = graph.clique(u);
      ++_cluster_sizes[u_id];
      for ( const NodeID& v : graph.neighbors(u) ) {
        const CliqueID v_id = graph.clique(v);
        if ( u_id == v_id && u < v /* count each internal edge only once */ ) {
          ++_internal_edges[u_id];
        }
      }
    }
  }

  void computeEmptyCliques(const Graph& graph) {
    computeClusterSizes(graph);
    _empty_cliques.clear();
    for ( const CliqueID& c : graph.nodes() ) {
      if ( _cluster_sizes[c] == 0 ) {
        _empty_cliques.push_back(c);
      }
    }
  }

  void computeNodesOfCliqueWithEmptyCliques(const Graph& graph) {
    for ( const NodeID& u : graph.nodes() ) {
      _cliques[u].clear();
    }
    for ( const NodeID& u : graph.nodes() ) {
      _cliques[graph.clique(u)].push_back(u);
    }
    _empty_cliques.clear();
    for ( const CliqueID& c : graph.nodes() ) {
      if ( _cliques[c].size() == 0 ) {
        _empty_cliques.push_back(c);
      }
    }
  }

 private:
  CommonOperations(const Graph& graph) :
    _cluster_sizes(graph.numNodes(), 0),
    _internal_edges(graph.numNodes(), 0),
    _empty_cliques(),
    _cliques(graph.numNodes()),
    _rating(3 * graph.maxDegree(), 0),
    _lp_aborted_flag(false),
    _is_special_instance(false) { }

  ~CommonOperations() = default;

  template<typename T>
  void reset(std::vector<T>& vec) {
    memset(vec.data(), 0, vec.size() * sizeof(T));
  }

 public:
  std::vector<NodeID> _cluster_sizes;
  std::vector<NodeID> _internal_edges;
  std::vector<CliqueID> _empty_cliques;
  std::vector<std::vector<NodeID>> _cliques;
  ds::FixedSizeSparseMap<CliqueID, EdgeWeight> _rating;
  bool _lp_aborted_flag;
  bool _is_special_instance;

};
}  // namespace utils
}  // namespace cluster_editing