#include "lp_refiner.h"

#include "cluster_editing/io/output.h"
#include "cluster_editing/utils/randomize.h"
#include "cluster_editing/utils/timer.h"
#include "cluster_editing/utils/progress_bar.h"
#include "cluster_editing/metrics.h"

namespace cluster_editing {

/**
 * Lessons learned:
 *  - LP Refiner converges slowly for large sparse graphs (still good improvements
 *    later iterations)
 *  - A move of a vertex can change the gain of a non-adjacent vertex, since
 *    number of insertions depends on cluster size
 *  - Tie breaking is important to achieve higher quality
 *  - Often finds exact solution for small instances
 */

void LabelPropagationRefiner::initializeImpl(Graph& graph) {
  _moved_vertices = 0;
  _clique_weight.assign(graph.numNodes(), 0);
  _empty_cliques.clear();
  _nodes.clear();
  _window_improvement = 0;
  _round_improvements.clear();
  for ( const NodeID& u : graph.nodes() ) {
    _nodes.push_back(u);
    const CliqueID from = graph.clique(u);
    ++_clique_weight[from];
  }

  for ( const NodeID& u : graph.nodes() ) {
    if ( _clique_weight[u] == 0 ) {
      _empty_cliques.push_back(u);
    }
  }
}

EdgeWeight LabelPropagationRefiner::refineImpl(Graph& graph) {
  utils::Timer::instance().start_timer("lp", "Label Propagation");
  bool converged = false;
  EdgeWeight start_metric =
    metrics::edge_deletions(graph) + metrics::edge_insertions(graph);
  EdgeWeight current_metric = start_metric;
  utils::ProgressBar lp_progress(
    _context.refinement.lp.maximum_lp_iterations, start_metric,
    _context.general.verbose_output && !debug);

  if ( _context.isTimeLimitReached() ) {
    return current_metric;
  }

  if ( _context.refinement.lp.node_order == NodeOrdering::random_shuffle ) {
    utils::Randomize::instance().shuffleVector(_nodes, _nodes.size());
  } else if ( _context.refinement.lp.node_order == NodeOrdering::degree_increasing ) {
    std::sort(_nodes.begin(), _nodes.end(), [&](const NodeID& u, const NodeID& v) {
      return graph.degree(u) < graph.degree(v);
    });
  } else if ( _context.refinement.lp.node_order == NodeOrdering::degree_decreasing ) {
    std::sort(_nodes.begin(), _nodes.end(), [&](const NodeID& u, const NodeID& v) {
      return graph.degree(u) > graph.degree(v);
    });
  }

  // enable early exit on large graphs, if FM refiner is used afterwards
  for ( int i = 0; i < _context.refinement.lp.maximum_lp_iterations && !converged; ++i ) {
    utils::Timer::instance().start_timer("local_moving", "Local Moving");
    converged = true;
    const EdgeWeight initial_metric = current_metric;
    for ( const NodeID& u : _nodes ) {
      const CliqueID from = graph.clique(u);
      Rating rating = computeBestTargetClique(graph, u, false);
      const CliqueID to = rating.clique;
      if ( to != from ) {
        moveVertex(graph, u, to);
        ++_moved_vertices;
        converged = false;

        // Verify, if rating is correct
        ASSERT(current_metric + rating.delta ==
          metrics::edge_deletions(graph) + metrics::edge_insertions(graph),
          "Rating is wrong. Expected:" <<
          (metrics::edge_deletions(graph) + metrics::edge_insertions(graph)) <<
          "but is" << (current_metric + rating.delta) << "(" <<
          V(current_metric) << "," << V(rating.delta) << ")");
        current_metric += rating.delta;
      }
    }
    utils::Timer::instance().stop_timer("local_moving");

    lp_progress.setObjective(current_metric);
    lp_progress += 1;
    DBG << "Pass Nr." << (i + 1) << "improved metric from"
        << initial_metric << "to" << current_metric
        << "( Moved Vertices:" << _moved_vertices << ")";

    if ( _context.refinement.lp.random_shuffle_each_round &&
         _context.refinement.lp.node_order == NodeOrdering::random_shuffle ) {
      utils::Timer::instance().start_timer("random_shuffle", "Random Shuffle");
      utils::Randomize::instance().shuffleVector(_nodes, _nodes.size());
      utils::Timer::instance().stop_timer("random_shuffle");
    }

    const EdgeWeight round_delta = initial_metric - current_metric;
    if ( round_delta > 0 ) {
      utils::Timer::instance().start_timer("checkpoint", "Checkpoint");
      graph.checkpoint(current_metric);
      utils::Timer::instance().stop_timer("checkpoint");
    }

    _window_improvement += round_delta;
    _round_improvements.push_back(round_delta);
    if ( _round_improvements.size() >= _context.refinement.lp.early_exit_window) {
      _window_improvement -= _round_improvements[
        _round_improvements.size() - _context.refinement.lp.early_exit_window];
      if ( _window_improvement <= _context.refinement.lp.min_improvement ) {
        break;
      }
    }

    if ( _context.isTimeLimitReached() ) {
      break;
    }
  }
  lp_progress += (_context.refinement.lp.maximum_lp_iterations - lp_progress.count());
  utils::Timer::instance().stop_timer("lp");
  return current_metric;
}

void LabelPropagationRefiner::moveVertex(Graph& graph, const NodeID u, const CliqueID to) {
  const CliqueID from = graph.clique(u);
  ASSERT(from != to);
  --_clique_weight[from];
  const bool from_becomes_empty = _clique_weight[from] == 0;
  const bool to_becomes_non_empty = _clique_weight[to] == 0;
  ++_clique_weight[to];
  graph.setClique(u, to);


  if ( to_becomes_non_empty ) {
    ASSERT(_empty_cliques.back() == to);
    _empty_cliques.pop_back();
  }
  if ( from_becomes_empty ) {
    _empty_cliques.push_back(from);
  }
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

LabelPropagationRefiner::Rating LabelPropagationRefiner::computeBestTargetClique(Graph& graph,
                                                                                 const NodeID u,
                                                                                 const bool force_isolation) {
  _rating.clear();
  Rating best_rating;
  const CliqueID from = graph.clique(u);
  best_rating.clique = from;
  _rating[from] = 0;

  for ( const Neighbor& n : graph.neighbors(u) ) {
    const CliqueID to = graph.clique(n.target);
    ++_rating[to];
  }

  const EdgeWeight u_degree = graph.degree(u);
  const EdgeWeight from_rating =
    insertions(_clique_weight[from] - 1 /* assumes u is not part of clique 'from'*/, _rating[from]) +
    deletions(u_degree, _rating[from]);
  best_rating.rating = from_rating;
  best_rating.delta = 0;
  std::vector<CliqueID> cliques_with_same_rating = { from };
  for ( const auto& entry : _rating ) {
    const CliqueID to = entry.key;
    if ( to != from ) {
      const EdgeWeight to_rating =
        insertions(_clique_weight[to], entry.value) +
        deletions(u_degree, entry.value);

      // It looks like that tie breaking is very important to achieve better quality
      if ( to_rating < best_rating.rating ) {
        best_rating.clique = to;
        best_rating.rating = to_rating;
        best_rating.delta = to_rating - from_rating;
        cliques_with_same_rating.clear();
        cliques_with_same_rating.push_back(to);
      } else if ( to_rating == best_rating.rating ) {
        cliques_with_same_rating.push_back(to);
      }
    }
  }

  // Check if it is beneficial to isolate the vertex again
  if ( !_empty_cliques.empty() && ( u_degree < best_rating.rating ||
       ( u_degree == best_rating.rating && force_isolation ) ) ) {
    best_rating.clique = _empty_cliques.back();
    best_rating.rating = u_degree;
    best_rating.delta = u_degree - from_rating;
    cliques_with_same_rating.clear();
    cliques_with_same_rating.push_back(_empty_cliques.back());
  } else if ( !_empty_cliques.empty() && u_degree == best_rating.rating ) {
    cliques_with_same_rating.push_back(_empty_cliques.back());
  }

  // Random tie breaking
  if ( cliques_with_same_rating.size() > 1 ) {
    const size_t tie_breaking_idx = utils::Randomize::instance().getRandomInt(
      0, cliques_with_same_rating.size() - 1);
    best_rating.clique = cliques_with_same_rating[tie_breaking_idx];
  }

  return best_rating;
}

} // namespace cluster_editing