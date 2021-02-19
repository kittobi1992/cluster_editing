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

namespace cluster_editing {

class LabelPropagationRefiner final : public IRefiner {
 private:

  static constexpr bool debug = true;

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
    _clique_weight(graph.numNodes()),
    _empty_cliques(),
    _rating(graph.numNodes()) { }

  LabelPropagationRefiner(const LabelPropagationRefiner&) = delete;
  LabelPropagationRefiner(LabelPropagationRefiner&&) = delete;

  LabelPropagationRefiner & operator= (const LabelPropagationRefiner &) = delete;
  LabelPropagationRefiner & operator= (LabelPropagationRefiner &&) = delete;

  size_t movedVertices() const {
    return _moved_vertices;
  }

 private:

  void initializeImpl(Graph& graph) final;

  bool refineImpl(Graph& graph) final ;

  void moveVertex(Graph& graph, const NodeID u, const CliqueID to);

  Rating computeBetTargetClique(Graph& graph, const NodeID u);

  const Context& _context;
  size_t _moved_vertices;
  std::vector<NodeID> _nodes;
  std::vector<NodeWeight> _clique_weight;
  std::vector<CliqueID> _empty_cliques;
  ds::SparseMap<CliqueID, EdgeWeight> _rating;
};
}  // namespace cluster_editing
