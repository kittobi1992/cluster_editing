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

#include "clique_remover.h"

#include "cluster_editing/io/output.h"
#include "cluster_editing/utils/randomize.h"
#include "cluster_editing/utils/timer.h"
#include "cluster_editing/utils/progress_bar.h"
#include "cluster_editing/metrics.h"

namespace cluster_editing {


void CliqueRemover::initializeImpl(Graph& graph) {
  utils::CommonOperations::instance(graph).computeNodesOfCliqueWithEmptyCliques(graph);
  utils::CommonOperations::instance(graph).computeClusterSizes(graph);
}

EdgeWeight CliqueRemover::refineImpl(Graph& graph,
                                     const EdgeWeight current_edits,
                                     const EdgeWeight) {
  utils::Timer::instance().start_timer("clique_remover", "Clique Remover");
  EdgeWeight edits = current_edits;

  if ( _context.isTimeLimitReached() ) {
    return current_edits;
  }

  utils::ProgressBar clique_remover_progress(
    _context.refinement.evo.clique_remover_iterations, edits,
    _context.general.verbose_output && !debug);
  for ( int i = 0; i < _context.refinement.evo.clique_remover_iterations; ++i ) {
    for ( const CliqueID& from : graph.nodes() ) {
      if ( _cliques[from].size() > 1 ) {
        EdgeWeight delta = 0;
        std::random_shuffle(_cliques[from].begin(), _cliques[from].end());
        for ( const NodeID& u : _cliques[from] ) {
          // Compute number of incident edges to all cliques
          _rating.clear();
          for ( const NodeID& v : graph.neighbors(u) ) {
            const CliqueID to = graph.clique(v);
            ++_rating[to];
          }

          // Compute best target clique (different to from)
          _best_cliques.clear();
          const EdgeWeight from_rating = _cluster_sizes[from] - 1 +
            graph.degree(u) - 2 * _rating[from];
          EdgeWeight best_rating = std::numeric_limits<EdgeWeight>::max();
          for ( const auto& entry : _rating ) {
            const CliqueID to = entry.key;
            if ( to != from ) {
              const EdgeWeight to_rating = _cluster_sizes[to] + graph.degree(u) - 2 * entry.value;
              if ( to_rating < best_rating ) {
                best_rating = to_rating;
                _best_cliques.clear();
                _best_cliques.push_back(to);
              } else if ( to_rating == best_rating ) {
                _best_cliques.push_back(to);
              }
            }
          }

          // Random tie breaking
          if ( !_best_cliques.empty() ) {
            const size_t tie_breaking_idx = utils::Randomize::instance().getRandomInt(
              0, _best_cliques.size() - 1);
            const CliqueID to = _best_cliques[tie_breaking_idx];
            delta += (best_rating - from_rating);
            graph.setClique(u, to);
            --_cluster_sizes[from];
            ++_cluster_sizes[to];
          }
        }

        if ( delta <= 0 ) {
          edits += delta;
          // Add nodes to target cliques
          size_t clique_size = _cliques[from].size();
          for ( size_t i = 0; i < clique_size; ++i ) {
            const NodeID u = _cliques[from][i];
            const CliqueID to = graph.clique(u);
            if ( from != to ) {
              std::swap(_cliques[from][i--], _cliques[from][--clique_size]);
              _cliques[from].pop_back();
              _cliques[to].push_back(u);
            }
          }
        } else {
          // If no improvement revert
          for ( const NodeID& u : _cliques[from] ) {
            ++_cluster_sizes[from];
            --_cluster_sizes[graph.clique(u)];
            graph.setClique(u, from);
          }
        }

        if ( _context.isTimeLimitReached() ) {
          break;
        }
      }
    }
    clique_remover_progress.setObjective(edits);
    clique_remover_progress += 1;
  }
  clique_remover_progress +=
    (_context.refinement.evo.clique_remover_iterations - clique_remover_progress.count());


  utils::Timer::instance().stop_timer("clique_remover");
  return edits;
}

} // namespace cluster_editing