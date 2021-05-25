#include "localized_evo.h"

#include <queue>

#include "cluster_editing/io/output.h"
#include "cluster_editing/utils/randomize.h"
#include "cluster_editing/utils/timer.h"
#include "cluster_editing/utils/progress_bar.h"
#include "cluster_editing/metrics.h"

namespace cluster_editing {


void LocalizedEvolutionary::initializeImpl(Graph& graph) {
  utils::CommonOperations::instance(graph).computeEmptyCliques(graph);
}

EdgeWeight LocalizedEvolutionary::refineImpl(Graph& graph,
                                             const EdgeWeight current_edits,
                                             const EdgeWeight) {
  utils::Timer::instance().start_timer("localized_evo", "Localized Evolutionary");
  utils::CommonOperations::instance(graph)._lp_aborted_flag = false;
  EdgeWeight start_metric = current_edits;
  EdgeWeight current_metric = current_edits;
  if ( _context.isTimeLimitReached() ) {
    return current_metric;
  }

  size_t steps = _context.refinement.localized_evo.run_until_time_limit ?
    std::numeric_limits<size_t>::max() : _context.refinement.localized_evo.steps;
  utils::ProgressBar lp_progress(
    steps, start_metric, _context.general.verbose_output && !debug);
  // Sorting-based rating is beneficial if the number of nodes is greater
  // than 150000
  const NodeID rating_map_degree_threshold = graph.numNodes() > 150000 ?
    _context.refinement.lp.rating_map_degree_threshold : 0;
  // enable early exit on large graphs, if FM refiner is used afterwards
  for ( size_t i = 0; i < steps; ++i ) {
    EdgeWeight delta = 0;
    // Mutate a small number of the vertices
    mutate(graph, delta);
    findRefinementNodes(graph);
    // Perform highly-localized label propagation around the mutated vertices
    for ( int j = 0; j < _context.refinement.localized_evo.max_lp_iterations; ++j ) {
      std::random_shuffle(_refinement_nodes.begin(), _refinement_nodes.end());
      for ( const NodeID& u : _refinement_nodes ) {
        const CliqueID from = graph.clique(u);
        const NodeID u_degree = graph.degree(u);
        if ( u_degree > 0 ) {
          Rating rating;
          if ( u_degree < rating_map_degree_threshold ) {
            rating = computeBestTargetCliqueWithSorting(graph, u);
          } else {
            rating = computeBestTargetCliqueWithRatingMap(graph, u);
          }
          const CliqueID to = rating.to;
          if ( to != from ) {
            moveVertex(graph, u, to);
            delta += rating.delta;
            _moves.push_back(rating);
          }
        }
      }
    }

    if ( delta <= 0 ) {
      // If improvement found -> update metric
      current_metric += delta;
    } else {
      // Otherwise rollback
      for ( int i = _moves.size() - 1; i >= 0; --i ) {
        const Rating& rating = _moves[i];
        moveVertex(graph, rating.u, rating.from);
      }
    }
    _moves.clear();
    lp_progress.setObjective(current_metric);
    lp_progress += 1;

    if ( _context.isTimeLimitReached() ) {
      break;
    }

    if ( i % 100 == 0 ) {
      const EdgeWeight delta = start_metric - current_metric;
      if ( delta > 0 ) {
        utils::Timer::instance().start_timer("checkpoint", "Checkpoint");
        graph.checkpoint(current_metric);
        utils::Timer::instance().stop_timer("checkpoint");
      }
    }
  }
  lp_progress += (_context.refinement.localized_evo.steps - lp_progress.count());
  utils::Timer::instance().stop_timer("localized_evo");

  return current_metric;
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

void LocalizedEvolutionary::mutate(Graph& graph,
                                   EdgeWeight& delta) {
  // Select vertices that we mutate
  _marked.reset();
  _mutation_nodes.clear();
  const int num_mutation_nodes = std::min(
    _context.refinement.localized_evo.num_mutations_nodes, static_cast<int>(graph.numNodes() / 4));
  utils::Randomize& rnd = utils::Randomize::instance();
  for ( int i = 0; i < num_mutation_nodes; ++i ) {
    NodeID u = rnd.getRandomInt(0, graph.numNodes() - 1);
    CliqueID from = graph.clique(u);
    if ( !_marked[u] && _clique_sizes[from] > 1 ) {
      _marked.set(u, true);
      _mutation_nodes.push_back(u);
      while ( rnd.getRandomFloat(0.0f, 1.0f) <=
           _context.refinement.localized_evo.choose_adjacent_mutation_node_prob ) {
        const NodeID v = graph.randomNeighbor(u);
        from = graph.clique(v);
        if ( v != INVALID_NODE && !_marked[v] && _clique_sizes[from] > 1 ) {
          u = v;
          _marked.set(u, true);
          _mutation_nodes.push_back(u);
          ++i;
        }
      }
    } else {
      --i;
    }
  }

  // Mutate selected vertices
  std::vector<NodeID> targets;
  for ( const NodeID& u : _mutation_nodes ) {
    _rating.clear();
    targets.clear();
    const CliqueID from = graph.clique(u);
    const NodeID u_degree = graph.degree(u);
    for ( const NodeID& v : graph.neighbors(u) ) {
      const CliqueID v_from = graph.clique(v);
      ++_rating[v_from];
      if ( v_from != from && _rating[v_from] == 1 ) {
        targets.push_back(v_from);
      }
    }

    const EdgeWeight from_rating =
      insertions(_clique_sizes[from] - 1, _rating[from]) +
      deletions(u_degree, _rating[from]);
    Rating rating;
    rating.to = INVALID_CLIQUE;
    if ( _clique_sizes[from] > 1 &&
        ( rnd.flipCoin() || targets.empty() ) ) {
      rating = isolateVertex(graph, u, from_rating);
    } else if ( !targets.empty() ) {
      if ( rnd.flipCoin() ) {
        targets.clear();
        EdgeWeight best_rating = std::numeric_limits<EdgeWeight>::max();
        for ( const auto& entry : _rating ) {
          const CliqueID to = entry.key;
          if ( to != from ) {
            const EdgeWeight to_rating =
              insertions(_clique_sizes[to], entry.value) +
              deletions(u_degree, entry.value);

            if ( to_rating < best_rating ) {
              best_rating = to_rating;
              targets.clear();
              targets.push_back(to);
            } else if ( to_rating == best_rating ) {
              targets.push_back(to);
            }
          }
        }
      }
      rating = moveToRandomTarget(graph, u, targets, from_rating);
    }

    if ( rating.to != INVALID_CLIQUE ) {
      _moves.emplace_back(rating);
      delta += rating.delta;
    }
  }
}

LocalizedEvolutionary::Rating LocalizedEvolutionary::isolateVertex(Graph& graph,
                                                                   const NodeID u,
                                                                   const EdgeWeight from_rating) {
  const NodeID u_degree = graph.degree(u);
  Rating rating;
  rating.u = u;
  rating.from = graph.clique(u);
  rating.to = _empty_cliques.back();
  rating.rating = u_degree;
  rating.delta = u_degree - from_rating;
  moveVertex(graph, u, rating.to);
  return rating;
}

LocalizedEvolutionary::Rating LocalizedEvolutionary::moveToRandomTarget(Graph& graph,
                                                                        const NodeID u,
                                                                        const std::vector<CliqueID>& targets,
                                                                        const EdgeWeight from_rating) {
  const NodeID u_degree = graph.degree(u);
  const size_t tie_breaking_idx = utils::Randomize::instance().getRandomInt(
    0, targets.size() - 1);
  Rating rating;
  rating.u = u;
  rating.from = graph.clique(u);
  rating.to = targets[tie_breaking_idx];
  rating.rating =
    insertions(_clique_sizes[rating.to], _rating[rating.to]) +
    deletions(u_degree, _rating[rating.to]);
  rating.delta = rating.rating - from_rating;
  moveVertex(graph, u, rating.to);
  return rating;
}

void LocalizedEvolutionary::findRefinementNodes(const Graph& graph) {
  _refinement_nodes.clear();
  _marked.reset();

  std::queue<NodeID> q;
  std::queue<NodeID> next_q;
  for ( const NodeID& u : _mutation_nodes ) {
    q.push(u);
    _marked.set(u, true);
  }

  utils::Randomize& rnd = utils::Randomize::instance();
  int distance = 0;
  while ( !q.empty() ) {
    const NodeID u = q.front();
    q.pop();
    const NodeID u_degree = graph.degree(u);
    _refinement_nodes.push_back(u);

    if ( distance < _context.refinement.localized_evo.max_distance_to_mutation_node ) {
      if ( static_cast<int>(u_degree) <= _context.refinement.localized_evo.degree_sampling_threshold ) {
        for ( const NodeID& v : graph.neighbors(u) ) {
          if ( !_marked[v] ) {
            next_q.push(v);
            _marked.set(v, true);
          }
        }
      } else {
        const float prob = static_cast<float>(_context.refinement.localized_evo.degree_sampling_threshold) / u_degree;
        for ( const NodeID& v : graph.neighbors(u) ) {
          const float p = rnd.getRandomFloat(0.0, 1.0);
          if ( p <= prob && !_marked[v] ) {
            next_q.push(v);
            _marked.set(v, true);
          }
        }
      }

      if ( q.empty() ) {
        std::swap(q, next_q);
        ++distance;
      }
    }
  }
}

void LocalizedEvolutionary::moveVertex(Graph& graph, const NodeID u, const CliqueID to) {
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

LocalizedEvolutionary::Rating LocalizedEvolutionary::computeBestTargetCliqueWithRatingMap(Graph& graph,
                                                                                          const NodeID u) {
  _rating.clear();
  Rating best_rating;
  const CliqueID from = graph.clique(u);
  best_rating.to = from;
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

  // Check if it is beneficial to isolate the vertex
  if ( !_empty_cliques.empty() &&  u_degree < best_rating.rating ) {
    best_rating.to = _empty_cliques.back();
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
    best_rating.to = _cliques_with_same_rating[tie_breaking_idx];
  }

  best_rating.u = u;
  best_rating.from = from;
  return best_rating;
}

LocalizedEvolutionary::Rating LocalizedEvolutionary::computeBestTargetCliqueWithSorting(Graph& graph,
                                                                                        const NodeID u) {
  // Sort incident edges according to their clique id
  graph.sortNeighborsByCliqueID(u);

  Rating best_rating;
  best_rating.to = INVALID_CLIQUE;
  best_rating.rating = std::numeric_limits<EdgeWeight>::max();
  Rating current_rating;
  current_rating.to = graph.clique(*graph.neighbors(u).begin());
  current_rating.rating = 0;

  _cliques_with_same_rating.clear();
  const EdgeWeight u_degree = graph.degree(u);
  const CliqueID from = graph.clique(u);
  EdgeWeight from_rating = std::numeric_limits<EdgeWeight>::max();
  auto computeRating = [&]() {
    const CliqueID to = current_rating.to;
    const EdgeWeight to_rating =
      insertions(_clique_sizes[to] - (to == from ? 1 : 0) , current_rating.rating) +
      deletions(u_degree, current_rating.rating);

    if ( to_rating < best_rating.rating ) {
      best_rating.to = to;
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
    if ( current_rating.to != to ) {
      computeRating();
      current_rating.to = to;
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
  if ( !_empty_cliques.empty() &&  u_degree < best_rating.rating ) {
    best_rating.to = _empty_cliques.back();
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
    best_rating.to = _cliques_with_same_rating[tie_breaking_idx];
  }

  best_rating.u = u;
  best_rating.from = from;
  best_rating.delta = best_rating.rating - from_rating;
  return best_rating;
}

} // namespace cluster_editing