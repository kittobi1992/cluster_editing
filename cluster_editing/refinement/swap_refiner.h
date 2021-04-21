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

#include "cluster_editing/refinement/i_refiner.h"
#include "cluster_editing/context/context.h"
#include "cluster_editing/datastructures/sparse_map.h"
#include "cluster_editing/datastructures/fast_reset_flag_array.h"

namespace cluster_editing {

class SwapRefiner final : public IRefiner {
 private:

  static constexpr bool debug = false;

  struct Rating {
    CliqueID to;
    NodeID swap_node;
    EdgeWeight delta;
    EdgeWeight weight_to_target_clique;
  };

  struct CliqueNode {
    NodeID u;
    EdgeWeight weight_to_current_clique;
  };

 public:
  explicit SwapRefiner(const Graph& graph,
                       const Context& context) :
    _context(context),
    _nodes(),
    _cliques(),
    _clique_weight(graph.numNodes()),
    _empty_cliques(),
    _rating(graph.numNodes()),
    _adjacent_nodes(graph.numNodes()),
    _window_improvement(0),
    _round_improvements() { }

  SwapRefiner(const SwapRefiner&) = delete;
  SwapRefiner(SwapRefiner&&) = delete;

  SwapRefiner & operator= (const SwapRefiner &) = delete;
  SwapRefiner & operator= (SwapRefiner &&) = delete;

  size_t movedVertices() const {
    return _moved_vertices;
  }

 private:

  void initializeImpl(Graph& graph) final;

  EdgeWeight refineImpl(Graph& graph) final ;

  bool moveVertex(Graph& graph, const NodeID u, const Rating& rating);

  Rating computeBestSwap(Graph& graph, const NodeID u);

  const Context& _context;
  size_t _moved_vertices;
  std::vector<NodeID> _nodes;
  std::vector<std::vector<CliqueNode>> _cliques;
  std::vector<NodeWeight> _clique_weight;
  std::vector<CliqueID> _empty_cliques;
  ds::SparseMap<CliqueID, EdgeWeight> _rating;
  ds::FastResetFlagArray<> _adjacent_nodes;
  EdgeWeight _window_improvement;
  std::vector<EdgeWeight> _round_improvements;
};
}  // namespace cluster_editing
