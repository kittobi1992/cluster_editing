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

namespace cluster_editing {

class Evolutionary final : public IRefiner {
 private:

  enum class Action : uint8_t {
    INTESIVATE = 0,
    MUTATE = 1,
    COMBINE = 2
  };

  enum class Mutation : uint8_t {
    LARGE_CLIQUE_ISOLATOR = 0,
    LARGE_CLIQUE_WITH_NEIGHBOR_ISOLATOR = 1,
    RANDOM_NODE_ISOLATOR = 2,
    RANDOM_NODE_MOVER = 3,
    NUM_MUTATIONS = 4
  };

  static constexpr bool debug = false;

  using Solution = std::vector<CliqueID>;

  struct SolutionStats {
    EdgeWeight edits;
    Solution* solution;
    bool no_intensivate = false;

    void reset() {
      no_intensivate = false;
    }
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
    _active_mutations() {
    ASSERT(_context.refinement.evo.enabled_mutations.size() ==
           static_cast<size_t>(Mutation::NUM_MUTATIONS));
    for ( size_t i = 0; i < _context.refinement.evo.enabled_mutations.size(); ++i ) {
      if ( _context.refinement.evo.enabled_mutations[i] == '1' ) {
        _active_mutations.push_back(static_cast<Mutation>(i));
      }
    }
  }

  Evolutionary(const Evolutionary&) = delete;
  Evolutionary(Evolutionary&&) = delete;

  Evolutionary & operator= (const Evolutionary &) = delete;
  Evolutionary & operator= (Evolutionary &&) = delete;

 private:

  void initializeImpl(Graph& graph) final;

  EdgeWeight refineImpl(Graph& graph) final ;

  void createInitialPopulation(Graph& graph);

  void evolutionaryStep(Graph& graph);

  void intensivate(Graph& graph, SolutionStats& stats);

  void mutate(Graph& graph, SolutionStats& stats);

  EdgeWeight refineSolution(Graph& graph,
                            const int lp_iterations,
                            const bool use_random_node_order,
                            const bool show_detailed_output);

  Action selectAction(const bool no_intensivate) const {
    const float p =
      utils::Randomize::instance().getRandomFloat(
        ( no_intensivate ? _context.refinement.evo.intensivate_prob : 0.0 ), 1.0);
    if ( p < _context.refinement.evo.intensivate_prob ) {
      return Action::INTESIVATE;
    } else if ( p < _context.refinement.evo.intensivate_prob +
                    _context.refinement.evo.mutate_prob ) {
      return Action::MUTATE;
    } else {
      return Action::COMBINE;
    }
  }

  void evictSolution(const Graph& graph,
                     const EdgeWeight edits) {
    sortSolutions();
    const EdgeWeight worst_edits = _population.back().edits;
    if ( edits < worst_edits ) {
      _population.back().edits = edits;
      _population.back().reset();
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
  std::vector<Mutation> _active_mutations;
};
}  // namespace cluster_editing
