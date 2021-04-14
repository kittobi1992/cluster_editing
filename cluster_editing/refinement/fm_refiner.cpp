#include "fm_refiner.h"

#include "cluster_editing/utils/randomize.h"
#include "cluster_editing/metrics.h"
#include "cluster_editing/utils/progress_bar.h"
#include "cluster_editing/utils/timer.h"

namespace cluster_editing {


void FMRefiner::initializeImpl(Graph& graph) {
  _moved_vertices = 0;
  _clique_weight.assign(graph.numNodes(), 0);

  _nodes.clear();
  for ( const NodeID& u : graph.nodes() ) {
    ++_clique_weight[graph.clique(u)];

    EdgeWeight weight = 0;
    const CliqueID clique = graph.clique(u);
    // TODO accelerate case with singleton init?
    for (const Neighbor& nb : graph.neighbors(u)) {
      if (graph.clique(nb.target) == clique) {
        ++weight;
      }
    }
    insertIntoCurrentClique(u, clique, weight);
  }
}

bool FMRefiner::refineImpl(Graph& graph) {
  EdgeWeight start_metric = metrics::edits(graph);
  EdgeWeight current_metric = start_metric;
  EdgeWeight round_delta = -1;

  utils::ProgressBar fm_progress(
    _context.refinement.fm.maximum_fm_iterations, start_metric,
    _context.general.verbose_output && !debug);

  utils::Timer::instance().start_timer("boundary_fm", "Boundary FM");
  for ( int round = 0; round < _context.refinement.fm.maximum_fm_iterations; ++round ) {
    round_delta = 0;
    EdgeWeight best_delta = 0;

    // init PQ and empty cliques
    _empty_cliques.clear();
    for (const NodeID u : graph.nodes()) {
      Rating rating = computeBestClique(graph, u);
      if (rating.clique != INVALID_CLIQUE) {
        pq.insert(u, rating.delta);
        insertIntoTargetClique(u, rating);
      }
      if (_clique_weight[u] == 0) {
        _empty_cliques.push_back(u);
      }
    }

    // perform moves
    size_t num_fruitless_moves = 0;
    while (!pq.empty() && num_fruitless_moves <= _context.refinement.fm.max_fruitless_moves) {
      EdgeWeight estimated_gain = pq.topKey();
      NodeID u = pq.top();

      Rating rating = computeBestClique(graph, u);
      if (rating.delta != estimated_gain) {
        // retry! 		or do KaHiP/Metis approach and just apply? only if improvement?
        pq.adjustKey(u, rating.delta);
        if (rating.clique != graph.clique(u)) {
          updateTargetClique(u, rating);
        }
        continue;
      }

      pq.deleteTop();
      removeFromTargetClique(u);
      const CliqueID from = graph.clique(u), to = rating.clique;
      if ( _clique_weight[from] == 1 && to == ISOLATE_CLIQUE ) {
        continue;
      }

      moveVertex(graph, u, to);
      round_delta += rating.delta;
      if (round_delta < best_delta) {
        best_delta = round_delta;
        // permanently keep all moves up to here
        moves.clear();
        num_fruitless_moves = 0;
      } else {
        // store for rollback later
        moves.push_back({ u, from, to });
        ++num_fruitless_moves;
      }
      assert(from != graph.clique(u));

      ASSERT(current_metric + round_delta == metrics::edits(graph),
        "Rating is wrong. Expected:" << metrics::edits(graph)
          << "but is" << (current_metric + round_delta));

      deltaGainUpdates(graph, u, from, to);

      #ifndef NDEBUG
      checkPQGains(graph);
      #endif
    }

    #ifndef NDEBUG
    checkCliqueWeights(graph);
    #endif

    clearPQ();
    rollback(graph);

    current_metric += best_delta;
    ASSERT(current_metric == metrics::edits(graph));

    if ( best_delta < 0 ) {
      utils::Timer::instance().start_timer("checkpoint", "Checkpoint");
      graph.checkpoint(current_metric);
      utils::Timer::instance().stop_timer("checkpoint");
    }

    #ifndef NDEBUG
    checkCliqueWeights(graph);
    #endif

    DBG << "Pass Nr." << (round + 1) << "improved metric from"
      << start_metric << "to" << current_metric;

    fm_progress.setObjective(current_metric);
    fm_progress += 1;
  }
  fm_progress += (_context.refinement.fm.maximum_fm_iterations - fm_progress.count());
  utils::Timer::instance().stop_timer("boundary_fm");

  return current_metric < start_metric;
}

void FMRefiner::moveVertex(Graph& graph, NodeID u, CliqueID to, bool manage_empty_cliques) {
  const CliqueID from = graph.clique(u);
  assert(from != to);

  if (manage_empty_cliques && to == ISOLATE_CLIQUE) {
    assert(!_empty_cliques.empty());
    to = _empty_cliques.back();
    assert(_clique_weight[to] == 0);
  }
  assert(manage_empty_cliques || to != ISOLATE_CLIQUE);

  --_clique_weight[from];
  const bool from_becomes_empty = _clique_weight[from] == 0;
  const bool to_becomes_non_empty = _clique_weight[to] == 0;
  ++_clique_weight[to];
  graph.setClique(u, to);

  removeFromCurrentClique(u, from);
  EdgeWeight weight_to_target_clique = 0;
  for (const Neighbor& nb : graph.neighbors(u)) {
    if (graph.clique(nb.target) == to) {
      ++weight_to_target_clique;
    }
  }
  insertIntoCurrentClique(u, to, weight_to_target_clique);

  if (manage_empty_cliques && to_becomes_non_empty) {
    assert(_empty_cliques.back() == to);
    _empty_cliques.pop_back();
  }
  if (manage_empty_cliques && from_becomes_empty) {
    _empty_cliques.push_back(from);
  }
  ++_moved_vertices;
}


void FMRefiner::deltaGainUpdates(const Graph& graph, const NodeID u, const CliqueID from, const CliqueID to) {
  // updates due to clique weight changes
  // these are a LOT of updates!
  for (NodeID v : target_cliques[from]) {
    assert(pq.contains(v));
    pq.adjustKey(v, pq.getKey(v) - 1);
  }
  if (to != ISOLATE_CLIQUE) {
    for (NodeID v : target_cliques[to]) {
      assert(pq.contains(v));
      pq.adjustKey(v, pq.getKey(v) + 1);
    }
  }

  for (NodeID v : current_cliques[from]) {
    if (pq.contains(v)) {
      pq.adjustKey(v, pq.getKey(v) + 1);
    }
  }
  if (to != ISOLATE_CLIQUE) {
    for (NodeID v : current_cliques[to]) {
      if (pq.contains(v)) {
        pq.adjustKey(v, pq.getKey(v) - 1);
      }
    }
  }

  static constexpr uint32_t SKIP_THRESHOLD = 20;

  const CliqueID actual_target = graph.clique(u);

  // update neighbors -- in original graph
  for (const Neighbor& nb : graph.neighbors(u)) {
    const NodeID v = nb.target;

    EdgeWeight delta = 0;
    if (graph.clique(v) == from) {
      --n[v].weight_to_current_clique;
      delta -= 2;
    } else if (graph.clique(v) == actual_target) {
      ++n[v].weight_to_current_clique;
      delta += 2;
    } else {
      // TODO could also sum up edge weight changes.
      // if weight changes allow for target cluster to change
      // (multiply by some treshold) --> recompute
      if (pq.contains(v) && ++n[v].num_skips > SKIP_THRESHOLD) {
        n[v].num_skips = 0;
        //LOG << "recalc" << V(v);
        Rating rv = computeBestClique(graph, v);
        pq.adjustKey(v, rv.delta);
        updateTargetClique(v, rv);
        continue;
      }
    }

    if (!pq.contains(v)) continue;

    if (n[v].desired_target == from) {
      --n[v].weight_to_target_clique;
      delta += 2;
    } else if (n[v].desired_target == actual_target) {
      ++n[v].weight_to_target_clique;
      delta -= 2;
    } else {
      if (++n[v].num_skips > SKIP_THRESHOLD) {
        n[v].num_skips = 0;
        //LOG << "recalc" << V(v);
        Rating rv = computeBestClique(graph, v);
        pq.adjustKey(v, rv.delta);
        updateTargetClique(v, rv);
        continue;
      }
    }

    if (delta != 0) {
      pq.adjustKey(v, pq.getKey(v) + delta);
    }
  }

  // special case for empty clique left behind
  if (_clique_weight[from] == 0) {
    for (NodeID v : target_cliques[from]) {
      assert(pq.contains(v));
      Rating rv = computeBestClique(graph, v);
      updateTargetClique(v, rv);
      pq.adjustKey(v, rv.delta);
      n[v].num_skips = 0;
    }
  }

  // special case for moving into isolated clique. see if
  // any neighbor wants to join
  if (to == ISOLATE_CLIQUE) {
    const CliqueID actual_to = graph.clique(u);
    for (const Neighbor& nb : graph.neighbors(u)) {
      const NodeID v = nb.target;
      if (!pq.contains(v)) continue;
      const EdgeWeight delta = gain(graph.clique(v), actual_to, n[v].weight_to_current_clique, 1);
      if (delta < pq.keyOf(v)) {
        pq.adjustKey(v, delta);
        Rating rv = { actual_to, 0, delta, 1 };
        updateTargetClique(v, rv);
        n[v].num_skips = 0;
      }
    }
  }
}

void FMRefiner::clearPQ() {
  while ( !pq.empty() ) {
    const NodeID u = pq.top();
    pq.deleteTop();
    removeFromTargetClique(u);
  }
}

void FMRefiner::rollback(Graph& graph) {
  for (Move& m : moves) {
    CliqueID from = graph.clique(m.node), to = m.from;
    for (const Neighbor& nb : graph.neighbors(m.node)) {
      const NodeID v = nb.target;
      if (graph.clique(v) == from) {
        --n[v].weight_to_current_clique;
      } else if (graph.clique(v) == to) {
        ++n[v].weight_to_current_clique;
      }
    }
    moveVertex(graph, m.node, to, false);
  }
  moves.clear();
}

FMRefiner::Rating FMRefiner::computeBestClique(const Graph& graph, const NodeID u) {
  edge_weight_to_clique.clear();
  for ( const Neighbor& nb : graph.neighbors(u) ) {
    const CliqueID v_c = graph.clique(nb.target);
    ++edge_weight_to_clique[v_c];
  }

  const EdgeWeight u_degree = graph.degree(u);
  const CliqueID from = graph.clique(u);
  const EdgeWeight from_rating =  insertions(_clique_weight[from] - 1, edge_weight_to_clique[from])
                                  + deletions(u_degree, edge_weight_to_clique[from]);

  // ignore the zero gain move of keeping u in its clique
  Rating best_rating = { INVALID_CLIQUE, std::numeric_limits<EdgeWeight>::max(),
                         std::numeric_limits<EdgeWeight>::max(), std::numeric_limits<EdgeWeight>::max() };

  for ( const auto& entry : edge_weight_to_clique ) {
    const CliqueID to = entry.key;
    if (to != from ) {
      const EdgeWeight to_rating =  insertions(_clique_weight[to], entry.value)
                                    + deletions(u_degree, entry.value);
      if (to_rating < best_rating.rating
          || (to_rating == best_rating.rating && utils::Randomize::instance().flipCoin())) {
        best_rating = { to, to_rating, to_rating - from_rating, entry.value };
      }
    }
  }

  // Check if it is beneficial to isolate the vertex again. if u is already isolated then this is not triggered
  if ( !_empty_cliques.empty() && u_degree < best_rating.rating ) {
    best_rating = { ISOLATE_CLIQUE, u_degree, u_degree - from_rating, 0 };
  }

  return best_rating;
}

void FMRefiner::checkPQGains(const Graph& graph) {
  for (size_t j = 0; j < pq.size(); ++j) {
    NodeID v = pq.at(j);
    EdgeWeight gain_in_pq = pq.keyAtPos(j);
    EdgeWeight gain_recalculated = gain(graph.clique(v), n[v].desired_target,
                                        n[v].weight_to_current_clique, n[v].weight_to_target_clique);

    edge_weight_to_clique.clear();
    for (const Neighbor& nb : graph.neighbors(v)) {
      ++edge_weight_to_clique[graph.clique(nb.target)];
    }

    if (n[v].desired_target == ISOLATE_CLIQUE) {
      assert(n[v].weight_to_target_clique == 0);
    } else {
      assert(n[v].weight_to_target_clique == edge_weight_to_clique[n[v].desired_target]);
    }
    assert(n[v].weight_to_current_clique == edge_weight_to_clique[graph.clique(v)]);

    if (gain_in_pq != gain_recalculated) {
      LOG << V(j) << V(v) << V(graph.clique(v)) << V(n[v].desired_target)
          << V(n[v].weight_to_current_clique) << V(n[v].weight_to_target_clique);
      LOG << V(gain_in_pq) << V(gain_recalculated);
      std::cout << "cliques: ";
      for (NodeID w : graph.nodes()) {
        std::cout << graph.clique(w) << " ";
      }
      std::cout << "\nneighbors: ";
      for (const Neighbor& nb : graph.neighbors(v)) {
        std::cout << nb.target << " " << graph.clique(nb.target) << " | ";
      }
      std::cout << std::endl;
    }
    assert(gain_in_pq == gain_recalculated);
  }
}

void FMRefiner::checkCliqueWeights(const Graph& graph) {
  for (NodeID u : graph.nodes()) {
    edge_weight_to_clique.clear();
    for (const Neighbor& nb : graph.neighbors(u)) {
      ++edge_weight_to_clique[graph.clique(nb.target)];
    }
    assert(n[u].weight_to_current_clique == edge_weight_to_clique[graph.clique(u)]);
  }
}

}