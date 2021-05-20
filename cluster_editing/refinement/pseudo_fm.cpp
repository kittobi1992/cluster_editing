#include "pseudo_fm.h"

#include "cluster_editing/io/output.h"
#include "cluster_editing/utils/randomize.h"
#include "cluster_editing/utils/timer.h"
#include "cluster_editing/utils/progress_bar.h"
#include "cluster_editing/metrics.h"

namespace cluster_editing {


void PseudoFM::initializeImpl(Graph& graph) {
  utils::CommonOperations::instance(graph).computeNodesOfCliqueWithEmptyCliques(graph);
  utils::CommonOperations::instance(graph).computeClusterSizes(graph);
}

EdgeWeight PseudoFM::refineImpl(Graph& graph,
                                const EdgeWeight current_edits,
                                const EdgeWeight) {
  utils::Timer::instance().start_timer("pseudo_fm", "Pseudo FM");
  EdgeWeight edits = current_edits;

  if ( _context.isTimeLimitReached() ) {
    return current_edits;
  }

  utils::ProgressBar pseudo_fm_progress(
    _context.refinement.evo.node_swapper_iterations, edits,
    _context.general.verbose_output && !debug);
  std::vector<Rating> moves;
  for ( int i = 0; i < _context.refinement.evo.node_swapper_iterations; ++i ) {

    for ( const CliqueID& from : graph.nodes() ) {
      if ( _cluster_sizes[from] > 1 ) {
        moves.clear();
        EdgeWeight delta = 0;

        Rating initial_rating = computeBestMoveOfClique(graph, from, INVALID_NODE);
        if ( initial_rating.to != INVALID_CLIQUE ) {
          moveVertex(graph, initial_rating.u, initial_rating.to);
          moves.push_back(initial_rating);
          delta += initial_rating.delta;

          CliqueID next_clique = initial_rating.to;
          for ( size_t i = 0; i < 5; ++i ) {
            if ( delta < 0 ) {
              break;
            }

            Rating rating = computeBestMoveOfClique(graph, next_clique, moves.back().u);
            if ( rating.to != INVALID_CLIQUE ) {
              moves.push_back(rating);
              delta += rating.delta;
              moveVertex(graph, rating.u, rating.to);
              next_clique = rating.to;
            } else {
              break;
            }
          }

          if ( delta <= 0 ) {
            edits += delta;
          } else {
            for ( int i = moves.size() - 1; i >= 0; --i ) {
              moveVertex(graph, moves[i].u, moves[i].from);
            }
          }
        }
      }
    }

    pseudo_fm_progress.setObjective(edits);
    pseudo_fm_progress += 1;
  }
  pseudo_fm_progress +=
    (_context.refinement.evo.clique_remover_iterations - pseudo_fm_progress.count());

  utils::Timer::instance().stop_timer("pseudo_fm");
  return edits;
}


void PseudoFM::moveVertex(Graph& graph, const NodeID u, const CliqueID to) {
  const CliqueID from = graph.clique(u);
  graph.setClique(u, to);
  --_cluster_sizes[from];
  ++_cluster_sizes[to];
  for ( size_t i = 0; i < _cliques[from].size(); ++i ) {
    if ( _cliques[from][i] == u ) {
      std::swap(_cliques[from][i], _cliques[from].back());
      _cliques[from].pop_back();
      break;
    }
  }
  _cliques[to].push_back(u);
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

PseudoFM::Rating PseudoFM::computeBestMoveOfClique(const Graph& graph,
                                                   const CliqueID from,
                                                   const NodeID forbidden_u) {
  _best_ratings.clear();
  for ( const NodeID& u : _cliques[from] ) {
    if ( u != forbidden_u ) {
      const Rating rating = computeBestTargetClique(graph, u);
      if ( _best_ratings.empty() ) {
        _best_ratings.push_back(rating);
      } else if ( rating.delta < _best_ratings.back().delta ) {
        _best_ratings.clear();
        _best_ratings.push_back(rating);
      } else if ( rating.delta == _best_ratings.back().delta ) {
        _best_ratings.push_back(rating);
      }
    }
  }

  Rating best_rating;
  best_rating.to = INVALID_CLIQUE;
  if ( !_best_ratings.empty() ) {
    const size_t tie_breaking_idx = utils::Randomize::instance().getRandomInt(
      0, _best_ratings.size() - 1);
    best_rating = _best_ratings[tie_breaking_idx];
  }

  return best_rating;
}

PseudoFM::Rating PseudoFM::computeBestTargetClique(const Graph& graph,
                                                   const NodeID u) {
  _rating.clear();
  Rating best_rating;
  const CliqueID from = graph.clique(u);
  best_rating.u = u;
  best_rating.from = from;
  _rating[from] = 0;

  for ( const NodeID& v : graph.neighbors(u) ) {
    const CliqueID to = graph.clique(v);
    ++_rating[to];
  }

  const EdgeWeight u_degree = graph.degree(u);
  const EdgeWeight from_rating =
    insertions(_cluster_sizes[from] - 1 /* assumes u is not part of clique 'from'*/, _rating[from]) +
    deletions(u_degree, _rating[from]);
  best_rating.to = INVALID_CLIQUE;
  best_rating.rating = std::numeric_limits<EdgeWeight>::max();
  best_rating.delta = std::numeric_limits<EdgeWeight>::max();
  _cliques_with_same_rating.clear();
  for ( const auto& entry : _rating ) {
    const CliqueID to = entry.key;
    if ( to != from && static_cast<int>(_cluster_sizes[to]) <= _context.refinement.evo.node_swapper_max_cluster_size ) {
      const EdgeWeight to_rating =
        insertions(_cluster_sizes[to], entry.value) +
        deletions(u_degree, entry.value);

      // It looks like that tie breaking is very important to achieve better quality
      if ( to_rating < best_rating.rating ) {
        best_rating.to = to;
        best_rating.rating = to_rating;
        best_rating.delta = to_rating - from_rating;
        _cliques_with_same_rating.clear();
        _cliques_with_same_rating.push_back(to);
      } else if ( to_rating == best_rating.rating ) {
        _cliques_with_same_rating.push_back(to);
      }
    }
  }

  // Random tie breaking
  if ( _cliques_with_same_rating.size() > 1 ) {
    const size_t tie_breaking_idx = utils::Randomize::instance().getRandomInt(
      0, _cliques_with_same_rating.size() - 1);
    best_rating.to = _cliques_with_same_rating[tie_breaking_idx];
  }

  return best_rating;
}

} // namespace cluster_editing