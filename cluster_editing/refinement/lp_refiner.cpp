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
 */

void LabelPropagationRefiner::initializeImpl(Graph& graph) {
  _moved_vertices = 0;
  _clique_weight.assign(graph.numNodes(), 0);
  _empty_cliques.clear();
  _nodes.clear();
  for ( const NodeID& u : graph.nodes() ) {
    _nodes.push_back(u);
    _clique_weight[graph.clique(u)] += graph.nodeWeight(u);
  }
}

bool LabelPropagationRefiner::refineImpl(Graph& graph) {

  bool converged = false;
  EdgeWeight start_metric =
    metrics::edge_deletions(graph) + metrics::edge_insertions(graph);
  EdgeWeight current_metric = start_metric;
  utils::ProgressBar lp_progress(
    _context.refinement.lp.maximum_lp_iterations, start_metric,
    _context.general.verbose_output && !debug);
  for ( int i = 0; i < _context.refinement.lp.maximum_lp_iterations && !converged; ++i ) {

    utils::Timer::instance().start_timer("random_shuffle", "Random Shuffle");
    utils::Randomize::instance().shuffleVector(_nodes, _nodes.size());
    utils::Timer::instance().stop_timer("random_shuffle");

    utils::Timer::instance().start_timer("local_moving", "Local Moving");
    converged = true;
    const EdgeWeight initial_metric = current_metric;
    for ( const NodeID& u : _nodes ) {
      Rating rating = computeBetTargetClique(graph, u);
      if ( rating.clique != graph.clique(u) ) {
        moveVertex(graph, u, rating.clique);
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
  }
  lp_progress += (_context.refinement.lp.maximum_lp_iterations - lp_progress.count());

  return current_metric < start_metric;
}

void LabelPropagationRefiner::moveVertex(Graph& graph, const NodeID u, const CliqueID to) {
  const CliqueID from = graph.clique(u);
  ASSERT(from != to);
  _clique_weight[from] -= graph.nodeWeight(u);
  const bool from_becomes_empty = _clique_weight[from] == 0;
  const bool to_becomes_non_empty = _clique_weight[to] == 0;
  _clique_weight[to] += graph.nodeWeight(u);
  graph.setClique(u, to);

  if ( from_becomes_empty ) {
    _empty_cliques.push_back(from);
  }
  if ( to_becomes_non_empty ) {
    ASSERT(_empty_cliques.back() == to);
    _empty_cliques.pop_back();
  }
}

namespace {

  ATTRIBUTE_ALWAYS_INLINE EdgeWeight insertions(const NodeWeight u_weight,
                                                const NodeWeight target_clique_weight,
                                                const EdgeWeight incident_edge_weight_to_target_clique) {
    return u_weight * target_clique_weight - incident_edge_weight_to_target_clique;
  }

  ATTRIBUTE_ALWAYS_INLINE EdgeWeight deletions(const EdgeWeight u_weight_degree,
                                               const EdgeWeight incident_edge_weight_to_target_clique) {
    return u_weight_degree - incident_edge_weight_to_target_clique;
  }

}

LabelPropagationRefiner::Rating LabelPropagationRefiner::computeBetTargetClique(Graph& graph, const NodeID u) {
  _rating.clear();
  Rating best_rating;
  const CliqueID from = graph.clique(u);
  best_rating.clique = from;
  _rating[from] = 0;

  for ( const Neighbor& n : graph.neighbors(u) ) {
    const CliqueID to = graph.clique(n.target);
    _rating[to] += graph.edgeWeight(n.id);
  }

  const EdgeWeight u_weighted_degree = graph.weightedDegree(u);
  const NodeWeight u_weight = graph.nodeWeight(u);
  const EdgeWeight from_rating =
    insertions(u_weight, _clique_weight[from] - u_weight /* assumes u is not part of clique 'from'*/, _rating[from]) +
    deletions(u_weighted_degree, _rating[from]);
  best_rating.rating = from_rating;
  best_rating.delta = 0;
  for ( const auto& entry : _rating ) {
    const CliqueID to = entry.key;
    if ( to != from ) {
      const EdgeWeight to_rating =
        insertions(u_weight, _clique_weight[to], entry.value) +
        deletions(u_weighted_degree, entry.value);

      // It looks like that tie breaking is very important to achieve better quality
      if (   to_rating < best_rating.rating ||
           ( to_rating == best_rating.rating && utils::Randomize::instance().flipCoin() ) ) {
        best_rating.clique = to;
        best_rating.rating = to_rating;
        best_rating.delta = to_rating - from_rating;
      }
    }
  }

  // Check if it is beneficial to isolate the vertex again
  if ( !_empty_cliques.empty() && u_weighted_degree < best_rating.rating ) {
    best_rating.clique = _empty_cliques.back();
    best_rating.rating = u_weighted_degree;
    best_rating.delta = u_weighted_degree - from_rating;
  }

  return best_rating;
}

} // namespace cluster_editing