#include "lp_refiner.h"

#include "cluster_editing/io/output.h"
#include "cluster_editing/utils/randomize.h"
#include "cluster_editing/utils/timer.h"
#include "cluster_editing/utils/progress_bar.h"
#include "cluster_editing/metrics.h"

namespace cluster_editing {


void LabelPropagationRefiner::initializeImpl(Graph& graph) {
  _moved_vertices = 0;
  _nodes.clear();
  _window_improvement = 0;
  _round_improvements.clear();
  utils::CommonOperations::instance(graph).computeEmptyCliques(graph);
  if ( _context.refinement.lp.node_order != NodeOrdering::none ) {
    for ( const NodeID u : graph.nodes() ) {
      _nodes.push_back(u);
    }
  }
}

EdgeWeight LabelPropagationRefiner::refineImpl(Graph& graph,
                                               const EdgeWeight current_edits,
                                               const EdgeWeight target_edits) {
  utils::Timer::instance().start_timer("lp", "Label Propagation");
  utils::CommonOperations::instance(graph)._lp_aborted_flag = false;
  _prefer_isolation = utils::CommonOperations::instance(graph)._is_special_instance;
  bool converged = false;
  EdgeWeight start_metric = current_edits;
  EdgeWeight current_metric = current_edits;
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

  // Sorting-based rating is beneficial if the number of nodes is greater
  // than 150000
  const NodeID rating_map_degree_threshold = graph.numNodes() > 150000 ?
    _context.refinement.lp.rating_map_degree_threshold : 0;
  auto process_vertex = [&](const NodeID u) {
    const CliqueID from = graph.clique(u);
    const NodeID u_degree = graph.degree(u);
    if ( u_degree > 0 ) {
      Rating rating;
      if ( u_degree < rating_map_degree_threshold ) {
        rating = computeBestTargetCliqueWithSorting(graph, u);
      } else {
        rating = computeBestTargetCliqueWithRatingMap(graph, u);
      }
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
  };

  // enable early exit on large graphs, if FM refiner is used afterwards
  int max_lp_iterations = _context.refinement.lp.maximum_lp_iterations;
  for ( int i = 0; i < max_lp_iterations && !converged; ++i ) {
    utils::Timer::instance().start_timer("local_moving", "Local Moving");
    converged = true;
    const EdgeWeight initial_metric = current_metric;
    if ( _context.refinement.lp.node_order == NodeOrdering::none ) {
      for ( const NodeID& u : graph.nodes() ) {
        process_vertex(u);
      }
    } else {
      for ( const NodeID& u : _nodes ) {
        process_vertex(u);
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

    if ( i + 1 == _context.refinement.lp.maximum_lp_iterations &&
         current_metric - 25 < target_edits ) {
      max_lp_iterations *= 2;
      lp_progress.reset();
    }

    const EdgeWeight remaining_rounds = (max_lp_iterations - (i + 1));
    if ( target_edits != 0 && remaining_rounds > 10 &&
         current_metric - std::max(3, round_delta) * remaining_rounds > target_edits ) {
      utils::CommonOperations::instance(graph)._lp_aborted_flag = true;
      break;
    }
  }
  lp_progress += (_context.refinement.lp.maximum_lp_iterations - lp_progress.count());

  const EdgeWeight delta = start_metric - current_metric;
  if ( delta > 0 ) {
    utils::Timer::instance().start_timer("checkpoint", "Checkpoint");
    graph.checkpoint(current_metric);
    utils::Timer::instance().stop_timer("checkpoint");
  }
  utils::Timer::instance().stop_timer("lp");

  return current_metric;
}

void LabelPropagationRefiner::moveVertex(Graph& graph, const NodeID u, const CliqueID to) {
  const CliqueID from = graph.clique(u);
  ASSERT(from != to);
  --_clique_sizes[from];
  const bool from_becomes_empty = _clique_sizes[from] == 0;
  const bool to_becomes_non_empty = _clique_sizes[to] == 0;
  ++_clique_sizes[to];
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

LabelPropagationRefiner::Rating LabelPropagationRefiner::computeBestTargetCliqueWithRatingMap(Graph& graph,
                                                                                              const NodeID u) {
  _rating.clear();
  Rating best_rating;
  const CliqueID from = graph.clique(u);
  best_rating.clique = from;
  _rating[from] = 0;

  for ( const NodeID& v : graph.neighbors(u) ) {
    const CliqueID to = graph.clique(v);
    ++_rating[to];
  }

  const EdgeWeight u_degree = graph.degree(u);
  const EdgeWeight from_rating =
    insertions(_clique_sizes[from] - 1 /* assumes u is not part of clique 'from'*/, _rating[from]) +
    deletions(u_degree, _rating[from]);
  best_rating.rating = from_rating;
  best_rating.delta = 0;
  _cliques_with_same_rating.clear();
  _cliques_with_same_rating.push_back(from);
  for ( const auto& entry : _rating ) {
    const CliqueID to = entry.key;
    if ( to != from ) {
      const EdgeWeight to_rating =
        insertions(_clique_sizes[to], entry.value) +
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

  // Random tie breaking
  if ( _cliques_with_same_rating.size() > 1 ) {
    const size_t tie_breaking_idx = utils::Randomize::instance().getRandomInt(
      0, _cliques_with_same_rating.size() - 1);
    best_rating.clique = _cliques_with_same_rating[tie_breaking_idx];
  }

  return best_rating;
}

LabelPropagationRefiner::Rating LabelPropagationRefiner::computeBestTargetCliqueWithSorting(Graph& graph,
                                                                                            const NodeID u) {
  // Sort incident edges according to their clique id
  graph.sortNeighborsByCliqueID(u);

  Rating best_rating;
  best_rating.clique = INVALID_CLIQUE;
  best_rating.rating = std::numeric_limits<EdgeWeight>::max();
  Rating current_rating;
  current_rating.clique = graph.clique(*graph.neighbors(u).begin());
  current_rating.rating = 0;

  _cliques_with_same_rating.clear();
  const EdgeWeight u_degree = graph.degree(u);
  const CliqueID from = graph.clique(u);
  EdgeWeight from_rating = std::numeric_limits<EdgeWeight>::max();
  auto computeRating = [&]() {
    const CliqueID to = current_rating.clique;
    const EdgeWeight to_rating =
      insertions(_clique_sizes[to] - (to == from ? 1 : 0) , current_rating.rating) +
      deletions(u_degree, current_rating.rating);

    if ( to_rating < best_rating.rating ) {
      best_rating.clique = to;
      best_rating.rating = to_rating;
      _cliques_with_same_rating.clear();
      _cliques_with_same_rating.push_back(to);
    } else if ( to_rating == best_rating.rating ) {
      _cliques_with_same_rating.push_back(to);
    }

    // From rating is required to calculate delta
    // for the number of edits if we apply the move
    if ( to == from ) {
      from_rating = to_rating;
    }
  };

  // All neighbors with same clique id are now in a
  // consecutive range within the adjacency list
  for ( const NodeID& v : graph.neighbors(u) ) {
    const CliqueID to = graph.clique(v);
    if ( current_rating.clique != to ) {
      computeRating();
      current_rating.clique = to;
      current_rating.rating = 0;
    }
    ++current_rating.rating;
  }
  computeRating();

  // From rating is required to calculate delta
  // for the number of edits if we apply the move
  // => fallback if u is not adjacent to from clique
  if ( from_rating == std::numeric_limits<EdgeWeight>::max() ) {
    from_rating = insertions(_clique_sizes[from] - 1, 0) + deletions(u_degree, 0);
  }

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

  // Random tie breaking
  if ( _cliques_with_same_rating.size() > 1 ) {
    const size_t tie_breaking_idx = utils::Randomize::instance().getRandomInt(
      0, _cliques_with_same_rating.size() - 1);
    best_rating.clique = _cliques_with_same_rating[tie_breaking_idx];
  }

  best_rating.delta = best_rating.rating - from_rating;
  return best_rating;
}

} // namespace cluster_editing