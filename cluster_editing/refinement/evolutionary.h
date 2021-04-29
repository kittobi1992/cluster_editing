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
#include "cluster_editing/refinement/lp_refiner.h"
#include "cluster_editing/context/context.h"
#include "cluster_editing/datastructures/sparse_map.h"
#include "cluster_editing/datastructures/fast_reset_flag_array.h"
#include "cluster_editing/utils/randomize.h"
#include "cluster_editing/refinement/mutation.h"
#include "cluster_editing/refinement/action_selector.h"

namespace cluster_editing {

class Evolutionary final : public IRefiner {
 private:

  static constexpr bool debug = false;

  using Solution = std::vector<CliqueID>;

  struct SolutionStats {
    EdgeWeight edits;
    Solution* solution;
  };

  using Population = std::vector<SolutionStats>;

 public:
  explicit Evolutionary(const Graph& graph,
                        const Context& context) :
    _context(context),
    _original_context(context),
    _show_detailed_output(context.general.verbose_output && context.refinement.evo.enable_detailed_output),
    _population(context.refinement.evo.solution_pool_size, SolutionStats {
      static_cast<EdgeWeight>(graph.numEdges()), nullptr}),
    _solutions(context.refinement.evo.solution_pool_size),
    _lp_refiner(graph, _context),
    _mutator(_context),
    _evo_action_selector( { EvoAction::INTESIVATE, EvoAction::MUTATE } ) { }

  Evolutionary(const Evolutionary&) = delete;
  Evolutionary(Evolutionary&&) = delete;

  Evolutionary & operator= (const Evolutionary &) = delete;
  Evolutionary & operator= (Evolutionary &&) = delete;

 private:

  void initializeImpl(Graph& graph) final;

  EdgeWeight refineImpl(Graph& graph, const EdgeWeight current_edits) final ;

  void createInitialPopulation(Graph& graph, const EdgeWeight current_edits);

  void evolutionaryStep(Graph& graph);

  EdgeWeight intensivate(Graph& graph, SolutionStats& stats);

  EdgeWeight mutate(Graph& graph, SolutionStats& stats);

  EdgeWeight refineSolution(Graph& graph,
                            const EdgeWeight current_edits,
                            const int lp_iterations,
                            const bool use_random_node_order,
                            const bool show_detailed_output);

  void evictSolution(const Graph& graph,
                     const EdgeWeight edits) {
    sortSolutions();
    const EdgeWeight worst_edits = _population.back().edits;
    if ( edits < worst_edits ) {
      _population.back().edits = edits;
      storeSolution(graph, *_population.back().solution);
      if ( _show_detailed_output) {
        LOG << "Evict worst solution in pool with" << worst_edits
            << "edits and store new solution with" << edits << "edits";
      }
    } else if ( _show_detailed_output ) {
      LOG << "New solution with" << edits << "edits is not"
          << "better than worst solution in the pool with"
          << worst_edits << "edits";
    }
  }

  void storeSolution(const Graph& graph,
                     Solution& solution) {
    ASSERT(graph.numNodes() == solution.size());
    for ( const NodeID& u : graph.nodes() ) {
      solution[u] = graph.clique(u);
    }
  }

  void applySolution(Graph& graph,
                     const Solution& solution) {
    ASSERT(graph.numNodes() == solution.size());
    for ( const NodeID& u : graph.nodes() ) {
      ASSERT(solution[u] != INVALID_CLIQUE);
      graph.setClique(u, solution[u]);
    }
  }

  void sortSolutions() {
    std::sort(_population.begin(), _population.end(),
      [&](const SolutionStats& lhs, const SolutionStats& rhs) {
        return lhs.edits < rhs.edits;
      });
  }

  Context _context;
  const Context& _original_context;
  const bool _show_detailed_output;
  Population _population;
  std::vector<Solution> _solutions;
  LabelPropagationRefiner _lp_refiner;
  Mutator _mutator;
  ActionSelector<EvoAction> _evo_action_selector;
};
}  // namespace cluster_editing
