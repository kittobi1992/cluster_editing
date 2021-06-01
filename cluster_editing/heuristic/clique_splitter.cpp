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

#include "clique_splitter.h"

#include "cluster_editing/io/output.h"
#include "cluster_editing/utils/randomize.h"
#include "cluster_editing/utils/timer.h"
#include "cluster_editing/utils/progress_bar.h"
#include "cluster_editing/metrics.h"

namespace cluster_editing {


void CliqueSplitter::initializeImpl(Graph& graph) {
  utils::CommonOperations::instance(graph).computeNodesOfCliqueWithEmptyCliques(graph);
  utils::CommonOperations::instance(graph).computeClusterSizes(graph);
}

EdgeWeight CliqueSplitter::refineImpl(Graph& graph,
                                     const EdgeWeight current_edits,
                                     const EdgeWeight) {
  utils::Timer::instance().start_timer("clique_splitter", "Clique Splitter");
  EdgeWeight edits = current_edits;

  if ( _context.isTimeLimitReached() ) {
    return current_edits;
  }

  utils::ProgressBar clique_splitter_progress(
    _context.refinement.evo.clique_remover_iterations, edits,
    _context.general.verbose_output && !debug);
  for ( int i = 0; i < _context.refinement.evo.clique_remover_iterations; ++i ) {
    for ( const CliqueID& from : graph.nodes() ) {
      if ( _cliques[from].size() > 2 ) {
        std::random_shuffle(_cliques[from].begin(), _cliques[from].end());

        // Find highest rated non-edge
        _non_edges.clear();
        for ( const NodeID& u : _cliques[from] ) {
          _marked.reset();
          _rating_to_from.clear();
          for ( const NodeID& v : graph.neighbors(u) ) {
            _marked.set(v, true);
            if ( from == graph.clique(v) ) {
              ++_rating_to_from[u];
            }
          }

          for ( const NodeID& v : _cliques[from] ) {
            if ( u < v && !_marked[v] ) {
              _non_edges.push_back(NonEdge { u, v });
            }
          }
        }

        // In case we found two non-edges, we try to split the clique
        if ( !_non_edges.empty() ) {
          std::sort(_non_edges.begin(), _non_edges.end(),
            [&](const NonEdge& lhs, const NonEdge& rhs) {
              return _rating_to_from[lhs.u] + _rating_to_from[rhs.v] <
                _rating_to_from[rhs.u] + _rating_to_from[rhs.v];
            });
          const NodeID seed_1 = _non_edges.back().u;
          const NodeID seed_2 = _non_edges.back().v;
          EdgeWeight delta = isolateAllVertices(graph, from, _cliques[from]);
          const CliqueID from_1 = graph.clique(seed_1);
          const CliqueID from_2 = graph.clique(seed_2);
          for ( const NodeID& u : _cliques[from] ) {
            if ( u != seed_1 && u != seed_2 ) {
              EdgeWeight edges_to_from_1 = 0;
              EdgeWeight edges_to_from_2 = 0;
              for ( const NodeID& v : graph.neighbors(u) ) {
                const CliqueID v_from = graph.clique(v);
                if ( from_1 == v_from ) {
                  ++edges_to_from_1;
                } else if ( from_2 == v_from ) {
                  ++edges_to_from_2;
                }
              }

              const CliqueID u_from = graph.clique(u);
              const NodeID u_degree = graph.degree(u);
              const EdgeWeight from_rating = u_degree;
              const EdgeWeight from_1_rating = _cluster_sizes[from_1] + u_degree - 2 * edges_to_from_1;
              const EdgeWeight from_2_rating = _cluster_sizes[from_2] + u_degree - 2 * edges_to_from_2;
              CliqueID to = INVALID_CLIQUE;
              if ( from_1_rating < from_2_rating ) {
                to = from_1;
                delta += (from_1_rating - from_rating);
              } else if ( from_1_rating > from_2_rating ) {
                to = from_2;
                delta += (from_2_rating - from_rating);
              } else {
                if ( utils::Randomize::instance().flipCoin() ) {
                  to = from_1;
                  delta += (from_1_rating - from_rating);
                } else {
                  to = from_2;
                  delta += (from_2_rating - from_rating);
                }
              }

              ASSERT(to != INVALID_CLIQUE);
              graph.setClique(u, to);
              ++_cluster_sizes[to];
              --_cluster_sizes[u_from];
              if ( _cluster_sizes[u_from] == 0 && u_from != from ) {
                _empty_cliques.push_back(u_from);
              }
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
            if ( _cluster_sizes[from] == 0 ) {
              _empty_cliques.push_back(from);
            }
          } else {
            // If no improvement revert
            for ( const NodeID& u : _cliques[from] ) {
              const CliqueID to = graph.clique(u);
              ++_cluster_sizes[from];
              --_cluster_sizes[to];
              if ( _cluster_sizes[to] == 0 ) {
                _empty_cliques.push_back(to);
              }
              graph.setClique(u, from);
            }
          }
        }

        if ( _context.isTimeLimitReached() ) {
          break;
        }
      }
    }
    clique_splitter_progress.setObjective(edits);
    clique_splitter_progress += 1;
  }
  clique_splitter_progress +=
    (_context.refinement.evo.clique_remover_iterations - clique_splitter_progress.count());


  utils::Timer::instance().stop_timer("clique_splitter");
  return edits;
}

EdgeWeight CliqueSplitter::isolateAllVertices(Graph& graph,
                                              const CliqueID from,
                                              const std::vector<NodeID>& clique) {
  EdgeWeight isolate_delta = 0;
  for ( const NodeID& u : clique ) {
    if ( _cluster_sizes[from] > 1 ) {
      const NodeID u_degree = graph.degree(u);
      EdgeWeight edges_to_from = 0;
      for ( const NodeID& v : graph.neighbors(u) ) {
        if ( from == graph.clique(v) ) {
          ++edges_to_from;
        }
      }

      const EdgeWeight from_rating = _cluster_sizes[from] - 1 + u_degree - 2 * edges_to_from;
      isolate_delta += (u_degree - from_rating);
      const CliqueID to = _empty_cliques.back();
      _empty_cliques.pop_back();
      graph.setClique(u, to);
      --_cluster_sizes[from];
      ++_cluster_sizes[to];
    }
  }
  return isolate_delta;
}

} // namespace cluster_editing