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

#include "cluster_editing/heuristic/i_refiner.h"
#include "cluster_editing/context/context.h"
#include "cluster_editing/utils/common_operations.h"
#include "cluster_editing/datastructures/fast_reset_flag_array.h"

namespace cluster_editing {

class NodeSwapper final : public IRefiner {
 private:

  static constexpr bool debug = false;

  struct Rating {
    CliqueID clique;
    EdgeWeight rating;
    EdgeWeight delta;
  };

 public:
  explicit NodeSwapper(const Graph& graph,
                       const Context& context) :
    _context(context),
    _cluster_sizes(utils::CommonOperations::instance(graph)._cluster_sizes),
    _empty_cliques(utils::CommonOperations::instance(graph)._empty_cliques),
    _cliques(utils::CommonOperations::instance(graph)._cliques),
    _rating(utils::CommonOperations::instance(graph)._rating),
    _cliques_with_same_rating(),
    _prefer_isolation(false)  { }

  NodeSwapper(const NodeSwapper&) = delete;
  NodeSwapper(NodeSwapper&&) = delete;

  NodeSwapper & operator= (const NodeSwapper &) = delete;
  NodeSwapper & operator= (NodeSwapper &&) = delete;

 private:

  void initializeImpl(Graph& graph) final;

  EdgeWeight refineImpl(Graph& graph,
                        const EdgeWeight current_edits,
                        const EdgeWeight target_edits) final ;

  Rating computeBestTargetClique(const Graph& graph,
                                 const NodeID u,
                                 const bool restrict_max_cluster_size,
                                 const bool consider_isolating_vertex);

  const Context& _context;
  std::vector<NodeID>& _cluster_sizes;
  std::vector<CliqueID>& _empty_cliques;
  std::vector<std::vector<NodeID>>& _cliques;
  ds::FixedSizeSparseMap<CliqueID, EdgeWeight>& _rating;
  std::vector<CliqueID> _cliques_with_same_rating;
  bool _prefer_isolation;
};
}  // namespace cluster_editing
