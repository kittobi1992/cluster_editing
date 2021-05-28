/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include "cluster_editing/heuristic/i_refiner.h"
#include "cluster_editing/context/context.h"
#include "cluster_editing/datastructures/fast_reset_flag_array.h"
#include "cluster_editing/utils/common_operations.h"

namespace cluster_editing {

class LabelPropagationRefiner final : public IRefiner {
 private:

  static constexpr bool debug = false;

  struct Rating {
    CliqueID clique;
    EdgeWeight rating;
    EdgeWeight delta;
  };

 public:
  explicit LabelPropagationRefiner(const Graph& graph,
                                   const Context& context) :
    _context(context),
    _nodes(),
    _clique_sizes(utils::CommonOperations::instance(graph)._cluster_sizes),
    _empty_cliques(utils::CommonOperations::instance(graph)._empty_cliques),
    _rating(utils::CommonOperations::instance(graph)._rating),
    _window_improvement(0),
    _round_improvements(),
    _cliques_with_same_rating(),
    _prefer_isolation(false) { }

  LabelPropagationRefiner(const LabelPropagationRefiner&) = delete;
  LabelPropagationRefiner(LabelPropagationRefiner&&) = delete;

  LabelPropagationRefiner & operator= (const LabelPropagationRefiner &) = delete;
  LabelPropagationRefiner & operator= (LabelPropagationRefiner &&) = delete;

  size_t movedVertices() const {
    return _moved_vertices;
  }

 private:

  void initializeImpl(Graph& graph) final;

  EdgeWeight refineImpl(Graph& graph,
                        const EdgeWeight current_edits,
                        const EdgeWeight target_edits) final ;

  void moveVertex(Graph& graph, const NodeID u, const CliqueID to);

  Rating computeBestTargetCliqueWithRatingMap(Graph& graph, const NodeID u);

  Rating computeBestTargetCliqueWithSorting(Graph& graph, const NodeID u);

  const Context& _context;
  size_t _moved_vertices;
  std::vector<NodeID> _nodes;
  std::vector<NodeID>& _clique_sizes;
  std::vector<CliqueID>& _empty_cliques;
  ds::FixedSizeSparseMap<CliqueID, EdgeWeight>& _rating;
  EdgeWeight _window_improvement;
  std::vector<EdgeWeight> _round_improvements;
  std::vector<CliqueID> _cliques_with_same_rating;
  bool _prefer_isolation;
};
}  // namespace cluster_editing
