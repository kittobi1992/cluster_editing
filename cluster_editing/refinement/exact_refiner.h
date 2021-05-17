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

#include <numeric>

#include "cluster_editing/refinement/i_refiner.h"
#include "cluster_editing/context/context.h"
#include "cluster_editing/datastructures/fast_reset_flag_array.h"
#include "cluster_editing/datastructures/sparse_map.h"

namespace cluster_editing {

class ExactRefiner final : public IRefiner {
 private:

  static constexpr bool debug = false;

 public:
  explicit ExactRefiner(const Graph& graph, const Context& context) :
    _context(context),
    _visited(graph.numNodes()),
    _in_queue(graph.numNodes()),
    _nodes(graph.numNodes()),
    _mapping(graph.numNodes()) {
    std::iota(_nodes.begin(), _nodes.end(), 0);
  }

  ExactRefiner(const ExactRefiner&) = delete;
  ExactRefiner(ExactRefiner&&) = delete;

  ExactRefiner & operator= (const ExactRefiner &) = delete;
  ExactRefiner & operator= (ExactRefiner &&) = delete;

 private:

  void initializeImpl(Graph& graph) final;

  EdgeWeight refineImpl(Graph& graph,
                        const EdgeWeight current_edits,
                        const EdgeWeight target_edits) final ;

  std::vector<NodeID> growSubgraph(const Graph& graph,
                                   const NodeID seed,
                                   const std::vector<std::vector<NodeID>>& cliques);

  const Context& _context;
  ds::FastResetFlagArray<> _visited;
  ds::FastResetFlagArray<> _in_queue;
  std::vector<NodeID> _nodes;
  ds::SparseMap<NodeID, NodeID> _mapping;
};
}  // namespace cluster_editing
