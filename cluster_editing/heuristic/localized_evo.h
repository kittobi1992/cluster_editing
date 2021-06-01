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
#include "cluster_editing/datastructures/fast_reset_flag_array.h"
#include "cluster_editing/utils/common_operations.h"

namespace cluster_editing {

class LocalizedEvolutionary final : public IRefiner {
 private:

  static constexpr bool debug = false;

  #define HIGH_DEGREE_THRESHOLD 1000

  struct Rating {
    NodeID u;
    CliqueID from;
    CliqueID to;
    EdgeWeight rating;
    EdgeWeight delta;
  };

 public:
  explicit LocalizedEvolutionary(const Graph& graph,
                                 const Context& context) :
    _context(context),
    _moves(),
    _clique_sizes(utils::CommonOperations::instance(graph)._cluster_sizes),
    _empty_cliques(utils::CommonOperations::instance(graph)._empty_cliques),
    _rating(utils::CommonOperations::instance(graph)._rating),
    _mutation_nodes(),
    _refinement_nodes(),
    _cliques_with_same_rating(),
    _marked(graph.numNodes()),
    _max_distance(context.refinement.localized_evo.max_distance_to_mutation_node),
    _max_mutation_nodes(context.refinement.localized_evo.max_mutations_nodes),
    _step(0),
    _max_steps(context.refinement.localized_evo.steps),
    _prng(420),
    _prefer_isolation(false) {
    if ( utils::CommonOperations::instance(graph)._is_special_instance ) {
      _max_distance = std::max(5, _max_distance);
      _max_mutation_nodes = std::min(5, _max_mutation_nodes);
    }
    if ( _context.refinement.localized_evo.run_until_time_limit ) {
      _max_steps = std::numeric_limits<size_t>::max();
    }
  }

  LocalizedEvolutionary(const LocalizedEvolutionary&) = delete;
  LocalizedEvolutionary(LocalizedEvolutionary&&) = delete;

  LocalizedEvolutionary & operator= (const LocalizedEvolutionary &) = delete;
  LocalizedEvolutionary & operator= (LocalizedEvolutionary &&) = delete;

  EdgeWeight performTimeLimitedEvoSteps(Graph& graph, double time_limit, EdgeWeight current_edits);

  bool done() const {
    return _step >= _max_steps;
  }

  void setDone() {
    _step = _max_steps;
  }

 private:

  void initializeImpl(Graph& graph) final;

  EdgeWeight refineImpl(Graph& graph,
                        const EdgeWeight current_edits,
                        const EdgeWeight target_edits) final ;

  void mutate(Graph& graph, EdgeWeight& delta);

  Rating isolateVertex(Graph& graph,
                       const NodeID u,
                       const EdgeWeight from_rating);

  Rating moveToRandomTarget(Graph& graph,
                            const NodeID u,
                            const std::vector<CliqueID>& targets,
                            const EdgeWeight from_rating);

  void findRefinementNodes(const Graph& graph);

  NodeID selectNonMarkedNeighbor(const Graph& graph,
                                 const NodeID u,
                                 const bool adjacent);

  void moveVertex(Graph& graph, const NodeID u, const CliqueID to);

  Rating computeBestTargetCliqueWithRatingMap(Graph& graph, const NodeID u);

  Rating computeBestTargetCliqueWithSorting(Graph& graph, const NodeID u);

  const Context& _context;
  std::vector<Rating> _moves;
  std::vector<NodeID>& _clique_sizes;
  std::vector<CliqueID>& _empty_cliques;
  ds::FixedSizeSparseMap<CliqueID, EdgeWeight>& _rating;
  std::vector<NodeID> _mutation_nodes;
  std::vector<NodeID> _refinement_nodes;
  std::vector<CliqueID> _cliques_with_same_rating;
  ds::FastResetFlagArray<> _marked;
  int _max_distance;
  int _max_mutation_nodes;
  size_t _step;
  size_t _max_steps;
  std::mt19937 _prng;
  bool _prefer_isolation;
};
}  // namespace cluster_editing
