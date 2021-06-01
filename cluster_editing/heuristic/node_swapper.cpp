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

#include "node_swapper.h"

#include "cluster_editing/io/output.h"
#include "cluster_editing/utils/randomize.h"
#include "cluster_editing/utils/timer.h"
#include "cluster_editing/utils/progress_bar.h"
#include "cluster_editing/metrics.h"

namespace cluster_editing {


void NodeSwapper::initializeImpl(Graph& graph) {
  utils::CommonOperations::instance(graph).computeNodesOfCliqueWithEmptyCliques(graph);
  utils::CommonOperations::instance(graph).computeClusterSizes(graph);
}

EdgeWeight NodeSwapper::refineImpl(Graph& graph,
                                   const EdgeWeight current_edits,
                                   const EdgeWeight) {
  utils::Timer::instance().start_timer("node_swapper", "Node Swapper");
  EdgeWeight edits = current_edits;
  _prefer_isolation = utils::CommonOperations::instance(graph)._is_special_instance;

  if ( _context.isTimeLimitReached() ) {
    return current_edits;
  }

  utils::ProgressBar node_swapper_progress(
    _context.refinement.evo.node_swapper_iterations, edits,
    _context.general.verbose_output && !debug);
  for ( int i = 0; i < _context.refinement.evo.node_swapper_iterations; ++i ) {

    for ( const NodeID& u : graph.nodes() ) {
      // Find best target clique for vertex u different from its current clique
      // => Can worsen solution quality
      const CliqueID from = graph.clique(u);
      Rating rating = computeBestTargetClique(graph, u, true, false);
      EdgeWeight delta = rating.delta;
      const CliqueID to = rating.clique;
      if ( to != INVALID_CLIQUE ) {
        graph.setClique(u, to);
        --_cluster_sizes[from];
        ++_cluster_sizes[to];

        // Try to move nodes from u desired clique out in order to improve solution quality
        std::random_shuffle(_cliques[to].begin(), _cliques[to].end());
        size_t to_size = _cliques[to].size();
        for ( size_t i = 0; i < to_size; ++i ) {
          const NodeID v = _cliques[to][i];
          rating = computeBestTargetClique(graph, v, false, true);
          if ( rating.clique != INVALID_CLIQUE && ( rating.delta < 0 ||
               ( rating.delta == 0 && utils::Randomize::instance().flipCoin() ) ) ) {
            delta += rating.delta;
            const CliqueID v_to = rating.clique;
            if ( v_to == _empty_cliques.back() ) {
              _empty_cliques.pop_back();
            }
            graph.setClique(v, v_to);
            --_cluster_sizes[to];
            ++_cluster_sizes[v_to];
            std::swap(_cliques[to][i--], _cliques[to][--to_size]);
          }
        }


        if ( delta <= 0 ) {
          // If solution quality improved, then apply moves to data structures ...
          edits += delta;
          for ( int i = _cliques[to].size() - 1; i >= static_cast<int>(to_size); --i ) {
            const NodeID v = _cliques[to][i];
            const CliqueID v_to = graph.clique(v);
            _cliques[v_to].push_back(v);
            _cliques[to].pop_back();
          }

          _cliques[to].push_back(u);
          for ( size_t i = 0; i < _cliques[from].size(); ++i ) {
            if ( _cliques[from][i] == u ) {
              std::swap(_cliques[from][i], _cliques[from].back());
              _cliques[from].pop_back();
              break;
            }
          }

          if ( _cluster_sizes[from] == 0 ) {
            _empty_cliques.push_back(from);
          }
        } else {
          // ... otherwise revert move
          graph.setClique(u, from);
          --_cluster_sizes[to];
          ++_cluster_sizes[from];

          for ( size_t i = to_size; i < _cliques[to].size(); ++i ) {
            const NodeID v = _cliques[to][i];
            const CliqueID v_to = graph.clique(v);
            graph.setClique(v, to);
            ++_cluster_sizes[to];
            --_cluster_sizes[v_to];
            if ( _cluster_sizes[v_to] == 0 ) {
              _empty_cliques.push_back(v_to);
            }
          }
        }
      }
    }

    node_swapper_progress.setObjective(edits);
    node_swapper_progress += 1;
  }
  node_swapper_progress +=
    (_context.refinement.evo.clique_remover_iterations - node_swapper_progress.count());

  utils::Timer::instance().stop_timer("node_swapper");
  return edits;
}


namespace {

  ATTRIBUTE_ALWAYS_INLINE EdgeWeight insertions(const NodeWeight target_clique_weight,
                                                const EdgeWeight incident_edge_weight_to_target_clique) {
    return target_clique_weight - incident_edge_weight_to_target_clique;
  }

  ATTRIBUTE_ALWAYS_INLINE EdgeWeight deletions(const EdgeWeight u_degree,
                                               const EdgeWeight incident_edge_weight_to_target_clique) {
    return u_degree - incident_edge_weight_to_target_clique;
  }

}

NodeSwapper::Rating NodeSwapper::computeBestTargetClique(const Graph& graph,
                                                         const NodeID u,
                                                         const bool restrict_max_cluster_size,
                                                         const bool consider_isolating_vertex) {
  _rating.clear();
  Rating best_rating;
  const CliqueID from = graph.clique(u);
  _rating[from] = 0;

  for ( const NodeID& v : graph.neighbors(u) ) {
    const CliqueID to = graph.clique(v);
    ++_rating[to];
  }

  const EdgeWeight u_degree = graph.degree(u);
  const EdgeWeight from_rating =
    insertions(_cluster_sizes[from] - 1 /* assumes u is not part of clique 'from'*/, _rating[from]) +
    deletions(u_degree, _rating[from]);
  best_rating.clique = INVALID_CLIQUE;
  best_rating.rating = std::numeric_limits<EdgeWeight>::max();
  best_rating.delta = std::numeric_limits<EdgeWeight>::max();
  _cliques_with_same_rating.clear();
  for ( const auto& entry : _rating ) {
    const CliqueID to = entry.key;
    if ( to != from && ( !restrict_max_cluster_size ||
        static_cast<int>(_cluster_sizes[to]) <= _context.refinement.evo.node_swapper_max_cluster_size ) ) {
      const EdgeWeight to_rating =
        insertions(_cluster_sizes[to], entry.value) +
        deletions(u_degree, entry.value);

      // It looks like that tie breaking is very important to achieve better quality
      if ( to_rating < best_rating.rating ) {
        best_rating.clique = to;
        best_rating.rating = to_rating;
        best_rating.delta = to_rating - from_rating;
        _cliques_with_same_rating.clear();
        _cliques_with_same_rating.push_back(to);
      } else if ( to_rating == best_rating.rating ) {
        _cliques_with_same_rating.push_back(to);
      }
    }
  }

  if ( consider_isolating_vertex ) {
  // Check if it is beneficial to isolate the vertex
    if ( !_empty_cliques.empty() && ( u_degree < best_rating.rating ||
        ( u_degree == best_rating.rating && _prefer_isolation ) ) ) {
      best_rating.clique = _empty_cliques.back();
      best_rating.rating = u_degree;
      best_rating.delta = u_degree - from_rating;
      _cliques_with_same_rating.clear();
      _cliques_with_same_rating.push_back(_empty_cliques.back());
    } else if ( !_empty_cliques.empty() && u_degree == best_rating.rating ) {
      _cliques_with_same_rating.push_back(_empty_cliques.back());
    }
  }

  // Random tie breaking
  if ( _cliques_with_same_rating.size() > 1 ) {
    const size_t tie_breaking_idx = utils::Randomize::instance().getRandomInt(
      0, _cliques_with_same_rating.size() - 1);
    best_rating.clique = _cliques_with_same_rating[tie_breaking_idx];
  }

  return best_rating;
}

} // namespace cluster_editing