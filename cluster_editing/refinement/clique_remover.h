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

namespace cluster_editing {

class CliqueRemover final : public IRefiner {
 private:

  static constexpr bool debug = false;

 public:
  explicit CliqueRemover(const Graph& graph,
                                   const Context& context) :
    _context(context),
    _cluster_sizes(utils::CommonOperations::instance(graph)._cluster_sizes),
    _cliques(utils::CommonOperations::instance(graph)._cliques),
    _rating(utils::CommonOperations::instance(graph)._rating),
    _best_cliques() { }

  CliqueRemover(const CliqueRemover&) = delete;
  CliqueRemover(CliqueRemover&&) = delete;

  CliqueRemover & operator= (const CliqueRemover &) = delete;
  CliqueRemover & operator= (CliqueRemover &&) = delete;

 private:

  void initializeImpl(Graph& graph) final;

  EdgeWeight refineImpl(Graph& graph,
                        const EdgeWeight current_edits,
                        const EdgeWeight target_edits) final ;

  const Context& _context;
  std::vector<NodeID>& _cluster_sizes;
  std::vector<std::vector<NodeID>>& _cliques;
  ds::FixedSizeSparseMap<CliqueID, EdgeWeight>& _rating;
  std::vector<CliqueID> _best_cliques;
};
}  // namespace cluster_editing
