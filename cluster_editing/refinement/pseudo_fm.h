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
#include "cluster_editing/utils/common_operations.h"
#include "cluster_editing/datastructures/fast_reset_flag_array.h"

namespace cluster_editing {

class PseudoFM final : public IRefiner {
 private:

  static constexpr bool debug = false;

  struct Rating {
    NodeID u;
    CliqueID from;
    CliqueID to;
    EdgeWeight rating;
    EdgeWeight delta;
  };

 public:
  explicit PseudoFM(const Graph& graph,
                       const Context& context) :
    _context(context),
    _cluster_sizes(utils::CommonOperations::instance(graph)._cluster_sizes),
    _cliques(utils::CommonOperations::instance(graph)._cliques),
    _rating(utils::CommonOperations::instance(graph)._rating),
    _cliques_with_same_rating(),
    _best_ratings()  { }

  PseudoFM(const PseudoFM&) = delete;
  PseudoFM(PseudoFM&&) = delete;

  PseudoFM & operator= (const PseudoFM &) = delete;
  PseudoFM & operator= (PseudoFM &&) = delete;

 private:

  void moveVertex(Graph& graph, const NodeID u, const CliqueID to);

  void initializeImpl(Graph& graph) final;

  EdgeWeight refineImpl(Graph& graph,
                        const EdgeWeight current_edits,
                        const EdgeWeight target_edits) final ;

  Rating computeBestMoveOfClique(const Graph& graph,
                                 const CliqueID from,
                                 const NodeID forbidden_u);

  Rating computeBestTargetClique(const Graph& graph,
                                 const NodeID u);

  const Context& _context;
  std::vector<NodeID>& _cluster_sizes;
  std::vector<std::vector<NodeID>>& _cliques;
  ds::FixedSizeSparseMap<CliqueID, EdgeWeight>& _rating;
  std::vector<CliqueID> _cliques_with_same_rating;
  std::vector<Rating> _best_ratings;
};
}  // namespace cluster_editing
