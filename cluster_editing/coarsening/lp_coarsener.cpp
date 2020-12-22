#include "lp_coarsener.h"

#include "cluster_editing/io/output.h"
#include "cluster_editing/utils/randomize.h"
#include "cluster_editing/utils/timer.h"
#include "cluster_editing/metrics.h"

namespace cluster_editing {

void LabelPropagationCoarsener::coarsenImpl() {

  while ( true ) {
    Graph& current_graph = currentGraph();
    std::vector<NodeID> nodes;
    _clique_weight.assign(current_graph.numNodes(), 0);
    for ( const NodeID& u : current_graph.nodes() ) {
      nodes.push_back(u);
      _clique_weight[current_graph.clique(u)] += current_graph.nodeWeight(u);
    }

    size_t moved_vertices = 0;
    bool converged = false;
    EdgeWeight current_metric =
      metrics::edge_deletions(current_graph) + metrics::edge_insertions(current_graph);
    utils::Timer::instance().start_timer("clustering", "Clustering");
    for ( int i = 0; i < _context.coarsening.maximum_lp_iterations && !converged; ++i ) {

      utils::Timer::instance().start_timer("random_shuffle", "Random Shuffle");
      utils::Randomize::instance().shuffleVector(nodes, nodes.size());
      utils::Timer::instance().stop_timer("random_shuffle");

      utils::Timer::instance().start_timer("local_moving", "Local Moving");
      converged = true;
      const EdgeWeight initial_metric = current_metric;
      for ( const NodeID& u : nodes ) {
        Rating rating = computeBestClique(current_graph, u);
        if ( rating.clique != current_graph.clique(u) ) {
          moveVertex(current_graph, u, rating.clique);
          ++moved_vertices;
          converged = false;

          // Verify, if rating is correct
          ASSERT(current_metric + rating.delta ==
            metrics::edge_deletions(current_graph) + metrics::edge_insertions(current_graph),
            "Rating is wrong. Expected:" <<
            (metrics::edge_deletions(current_graph) + metrics::edge_insertions(current_graph)) <<
            "but is" << (current_metric + rating.delta) << "(" <<
            V(current_metric) << "," << V(rating.delta) << ")");
          current_metric += rating.delta;
        }
      }
      utils::Timer::instance().stop_timer("local_moving");

      DBG << "Pass Nr." << (i + 1) << "improved metric from"
          << initial_metric << "to" << current_metric;
    }
    utils::Timer::instance().stop_timer("clustering");

    if ( moved_vertices > 0 ) {
      utils::Timer::instance().start_timer("contraction", "Contraction");
      _hierarchies.emplace_back(current_graph.contract());
      utils::Timer::instance().stop_timer("contraction");

      if ( debug ) {
        io::printGraphInfo(_hierarchies.back().first, _context, "Contracted Graph");
      }
    } else {
      break;
    }
  }

  if ( _context.general.verbose_output ) {
    const EdgeWeight edge_insertions = metrics::edge_insertions(currentGraph());
    const EdgeWeight edge_deletions = metrics::edge_deletions(currentGraph());
    LOG << "Initial Solution:";
    LOG << "Number of Hierarchies:" << _hierarchies.size();
    LOG << "Edge Insertions:" << edge_insertions << "-"
        << "Edge Deletions:" << edge_deletions << "-"
        << "Edge Modifications:" << (edge_deletions + edge_insertions);
    io::printGraphInfo(_hierarchies.back().first, _context, "Contracted Graph");
  }
}

namespace {

  void refine(Graph& graph, std::unique_ptr<IRefiner>& refiner) {
    utils::Timer::instance().start_timer("initialize_refiner", "Initialize Refiner");
    refiner->initialize(graph);
    utils::Timer::instance().stop_timer("initialize_refiner");

    utils::Timer::instance().start_timer("refine", "Refine");
    refiner->refine(graph);
    utils::Timer::instance().stop_timer("refine");
  }

  void project(Graph& original_graph,
               const Graph& coarse_graph,
               const std::vector<NodeID>& clique_to_coarse) {
    utils::Timer::instance().start_timer("project_clique", "Projecting Cliques");
    for ( const NodeID& u : original_graph.nodes() ) {
      const CliqueID& u_c = original_graph.clique(u);
      original_graph.setClique(u, coarse_graph.clique(clique_to_coarse[u_c]));
    }
    utils::Timer::instance().stop_timer("project_clique");
  }

} // namespace

void LabelPropagationCoarsener::uncoarsenImpl(std::unique_ptr<IRefiner>& refiner) {
  // Refine coarsest graph
  refine(currentGraph(), refiner);

  while ( true ) {
    const Graph& coarse_graph = currentGraph();
    if ( &coarse_graph != &_graph ) {
      // Project partition to next level finer graph
      const size_t h_size = _hierarchies.size();
      if ( h_size > 1 ) {
        project(_hierarchies[h_size - 2].first, coarse_graph, _hierarchies.back().second);
      } else {
        project(_graph, coarse_graph, _hierarchies.back().second);
      }
      _hierarchies.pop_back();

      // Refine current graph
      refine(currentGraph(), refiner);
    } else {
      break;
    }
  }
}


namespace {

  ATTRIBUTE_ALWAYS_INLINE EdgeWeight insertions(const NodeWeight u_weight,
                                                const NodeWeight clique_weight,
                                                const EdgeWeight incident_edge_weight) {
    return u_weight * clique_weight - incident_edge_weight;
  }

}

void LabelPropagationCoarsener::moveVertex(Graph& graph, const NodeID u, const CliqueID to) {
  const CliqueID from = graph.clique(u);
  ASSERT(from != to);
  _clique_weight[from] -= graph.nodeWeight(u);
  _clique_weight[to] += graph.nodeWeight(u);
  graph.setClique(u, to);
}

LabelPropagationCoarsener::Rating LabelPropagationCoarsener::computeBestClique(Graph& graph, const NodeID u) {
  _rating.clear();
  Rating best_rating;
  const CliqueID u_c = graph.clique(u);
  best_rating.clique = u_c;
  _rating[u_c] = 0;

  for ( const Neighbor& n : graph.neighbors(u) ) {
    const CliqueID v_c = graph.clique(n.target);
    _rating[v_c] += graph.edgeWeight(n.id);
  }

  const EdgeWeight u_weighted_degree = graph.weightedDegree(u);
  const NodeWeight u_weight = graph.nodeWeight(u);
  const EdgeWeight current_rating = insertions(u_weight, _clique_weight[u_c] -
    u_weight, _rating[u_c]) + (u_weighted_degree - _rating[u_c]) /* deletions */;
  best_rating.rating = current_rating;
  best_rating.delta = 0;
  for ( const auto& entry : _rating ) {
    const CliqueID target_c = entry.key;
    if ( target_c != u_c ) {
      const EdgeWeight rating = insertions(u_weight, _clique_weight[target_c],
        entry.value) + (u_weighted_degree - entry.value) /* deletions */;
      if (   rating < best_rating.rating ||
           ( rating == best_rating.rating && utils::Randomize::instance().flipCoin() )) {
        best_rating.clique = target_c;
        best_rating.rating = rating;
        best_rating.delta = rating - current_rating;
      }
    }
  }

  return best_rating;
}

} // namespace cluster_editing