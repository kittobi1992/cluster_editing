#include "swap_refiner.h"

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

void SwapRefiner::initializeImpl(Graph& graph) {
  _moved_vertices = 0;
  _cliques.clear();
  _cliques.resize(graph.numNodes());
  _clique_weight.assign(graph.numNodes(), 0);
  _empty_cliques.clear();
  _nodes.clear();
  for ( const NodeID& u : graph.nodes() ) {
    _nodes.push_back(u);
    const CliqueID from = graph.clique(u);
    ++_clique_weight[from];
    CliqueNode node { u, 0 };
    for ( const Neighbor& n : graph.neighbors(u) ) {
      if ( graph.clique(n.target) == from ) {
        ++node.weight_to_current_clique;
      }
    }
    _cliques[from].push_back(node);
  }

  for ( const NodeID& u : graph.nodes() ) {
    if ( _clique_weight[u] == 0 ) {
      _empty_cliques.push_back(u);
    }
  }
}

EdgeWeight SwapRefiner::refineImpl(Graph& graph) {
  utils::Timer::instance().start_timer("swap_refiner", "Swap Refiner");
  bool converged = false;
  EdgeWeight start_metric =
    metrics::edge_deletions(graph) + metrics::edge_insertions(graph);
  EdgeWeight current_metric = start_metric;
  utils::ProgressBar swap_progress(
    _context.refinement.lp.maximum_lp_iterations, start_metric,
    _context.general.verbose_output && !debug);

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
  const bool enable_early_exit = _context.refinement.use_boundary_fm_refiner || _context.refinement.use_localized_fm_refiner;
  for ( int i = 0; i < _context.refinement.lp.maximum_lp_iterations && !converged; ++i ) {
    utils::Timer::instance().start_timer("local_moving", "Local Moving");
    converged = true;
    const EdgeWeight initial_metric = current_metric;
    for ( const NodeID& u : _nodes ) {
      const CliqueID from = graph.clique(u);
      Rating rating = computeBestSwap(graph, u);
      if ( rating.to != from && moveVertex(graph, u, rating) ) {
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
        // LOG << V(current_metric) << V(metrics::edits(graph));
      }
    }
    utils::Timer::instance().stop_timer("local_moving");

    swap_progress.setObjective(current_metric);
    swap_progress += 1;
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
      if ( enable_early_exit && _window_improvement <= _context.refinement.lp.min_improvement ) {
        break;
      }
    }
  }
  swap_progress += (_context.refinement.lp.maximum_lp_iterations - swap_progress.count());
  utils::Timer::instance().stop_timer("swap_refiner");
  return current_metric;
}


bool SwapRefiner::moveVertex(Graph& graph, const NodeID u, const Rating& rating) {
  if ( !_empty_cliques.empty() ) {
    const CliqueID from = graph.clique(u);
    const CliqueID to = rating.to;
    ASSERT(from != to);
    ASSERT(graph.clique(rating.swap_node) == to);

    // Isolate swap node
    const CliqueID swap_node_target = _empty_cliques.back();
    _empty_cliques.pop_back();
    graph.setClique(rating.swap_node, swap_node_target);
    _cliques[swap_node_target].push_back(CliqueNode { rating.swap_node, 0 });
    ++_clique_weight[swap_node_target];

    for ( size_t i = 0; i < _cliques[to].size(); ++i ) {
      if ( rating.swap_node == _cliques[to][i].u ) {
        std::swap(_cliques[to][i], _cliques[to].back());
        _cliques[to].pop_back();
        break;
      }
    }
    _adjacent_nodes.reset();
    for ( const Neighbor& n : graph.neighbors(rating.swap_node) ) {
      _adjacent_nodes.set(n.target, true);
    }
    for ( CliqueNode& c : _cliques[to] ) {
      if ( _adjacent_nodes[c.u] ) {
        --c.weight_to_current_clique;
      }
    }

    // Move u to target clique
    const bool from_becomes_empty = _clique_weight[from] == 1;
    --_clique_weight[from];
    graph.setClique(u, to);
    _cliques[to].push_back(CliqueNode { u, rating.weight_to_target_clique });
    if ( from_becomes_empty ) {
      _empty_cliques.push_back(from);
    }

    for ( size_t i = 0; i < _cliques[from].size(); ++i ) {
      if ( u == _cliques[from][i].u ) {
        std::swap(_cliques[from][i], _cliques[from].back());
        _cliques[from].pop_back();
        break;
      }
    }

    _adjacent_nodes.reset();
    for ( const Neighbor& n : graph.neighbors(u) ) {
      _adjacent_nodes.set(n.target, true);
    }
    for ( CliqueNode& c : _cliques[from] ) {
      if ( _adjacent_nodes[c.u] && c.u != u ) {
        --c.weight_to_current_clique;
      }
    }
    for ( CliqueNode& c : _cliques[to] ) {
      if ( _adjacent_nodes[c.u] && c.u != u ) {
        ++c.weight_to_current_clique;
      }
    }

    return true;
  }
  return false;
}

namespace {

  ATTRIBUTE_ALWAYS_INLINE EdgeWeight insertions(const NodeWeight target_clique_weight,
                                                const EdgeWeight incident_edge_weight_to_target_clique) {
    return std::max(target_clique_weight - incident_edge_weight_to_target_clique, 0);
  }

  ATTRIBUTE_ALWAYS_INLINE EdgeWeight deletions(const EdgeWeight u_degree,
                                               const EdgeWeight incident_edge_weight_to_target_clique) {
    return u_degree - incident_edge_weight_to_target_clique;
  }

}


SwapRefiner::Rating SwapRefiner::computeBestSwap(Graph& graph, const NodeID u) {
  _rating.clear();
  _adjacent_nodes.reset();
  Rating best_rating;
  const CliqueID from = graph.clique(u);
  _rating[from] = 0;

  for ( const Neighbor& n : graph.neighbors(u) ) {
    const CliqueID to = graph.clique(n.target);
    ++_rating[to];
    _adjacent_nodes.set(n.target, true);
  }

  const EdgeWeight u_degree = graph.degree(u);
  const EdgeWeight from_rating =
    insertions(_clique_weight[from] - 1 /* assumes u is not part of clique 'from'*/, _rating[from]) +
    deletions(u_degree, _rating[from]);
  best_rating.to = from;
  best_rating.swap_node = INVALID_NODE;
  best_rating.delta = 0;
  for ( const auto& entry : _rating ) {
    const CliqueID to = entry.key;
    if ( to != from && _clique_weight[to] > 2 ) {
      // Assume that we can potentially remove one vertex from the target clique
      EdgeWeight weight_to_target_clique = entry.value;
      EdgeWeight to_rating =
        insertions(_clique_weight[to] - 1, weight_to_target_clique) +
        deletions(u_degree, weight_to_target_clique);

      if ( to_rating < from_rating ) {
        for ( const CliqueNode& v : _cliques[to] ) {
          const bool is_adjacent_to_u = _adjacent_nodes[v.u];
          const EdgeWeight isolation_delta =
            2 * v.weight_to_current_clique - (_clique_weight[to] - 1);
          weight_to_target_clique -= is_adjacent_to_u;
          to_rating =
            insertions(_clique_weight[to] - 1, weight_to_target_clique) +
            deletions(u_degree, weight_to_target_clique);
          const EdgeWeight swap_delta =
            ( to_rating - from_rating ) + isolation_delta;
          if ( swap_delta < 0 ) {
            // LOG << V(from_rating) << V(to_rating) << V(weight_to_target_clique) << V(swap_delta);
            // LOG << V(v.weight_to_current_clique) << V(_clique_weight[to]) << V(_cliques[to].size()) << V(isolation_delta);
            return Rating { to, v.u, swap_delta, weight_to_target_clique };
          } else if ( swap_delta == 0 ) {
            best_rating = Rating { to, v.u, swap_delta, weight_to_target_clique };
          }
        }
      }
    }
  }


  return best_rating;
}

} // namespace cluster_editing