#include "fm_refiner.h"


#include "cluster_editing/utils/randomize.h"
#include "cluster_editing/metrics.h"
#include "cluster_editing/utils/timer.h"

namespace cluster_editing {


void FMRefiner::initializeImpl(Graph& graph) {
	_moved_vertices = 0;
	_clique_weight.assign(graph.numNodes(), 0);
	_empty_cliques.clear();
	_nodes.clear();
	for ( const NodeID& u : graph.nodes() ) {
		_clique_weight[graph.clique(u)] += graph.nodeWeight(u);

		EdgeWeight weight = 0;
		const CliqueID clique = graph.clique(u);
		// TODO accelerate case with singleton init?
		for (const Neighbor& nb : graph.neighbors(u)) {
			if (graph.clique(nb.target) == clique) {
				weight += graph.edgeWeight(nb.id);
			}
		}
		insertIntoCurrentClique(u, clique, weight);
	}
	for (NodeID u : graph.nodes()) {
		if (_clique_weight[u] == 0) {
			_empty_cliques.push_back(u);
		}
	}
}

bool FMRefiner::refineImpl(Graph& graph) {
	EdgeWeight start_metric = metrics::edits(graph);
	EdgeWeight current_metric = start_metric;
	EdgeWeight round_delta = -1;

	for ( int i = 0; i < _context.refinement.lp.maximum_lp_iterations && round_delta < 0; ++i ) {
		round_delta = 0;
		EdgeWeight best_delta = 0;

		// init PQ
		for (const NodeID u : graph.nodes()) {
			Rating rating = computeBestClique(graph, u);
			pq.insert(u, rating.delta);
			insertIntoTargetClique(u, rating);
		}

		// perform moves
		while (!pq.empty()) {
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

			const CliqueID from = graph.clique(u), to = rating.clique;
			round_delta += rating.delta;
			if (round_delta < best_delta) {
				// permanently apply all moves
				moves.clear();
			} else {
				// store for rollback later
				moves.push_back({ u, from, to });
			}

			removeFromTargetClique(u);
			moveVertex(graph, u, to);
			ASSERT(current_metric + round_delta == metrics::edits(graph), "Rating is wrong. Expected:" << metrics::edits(graph) << "but is" << (current_metric + round_delta));

			static constexpr uint32_t SKIP_THRESHOLD = 200;

			// update neighbors -- in original graph
			for (const Neighbor& nb : graph.neighbors(u)) {
				const NodeID v = nb.target;
				const EdgeWeight we = graph.edgeWeight(nb.id);
				if (graph.clique(v) == from) {
					n[v].weight_to_current_clique -= we;
					if (pq.contains(v)) {
						pq.adjustKey(v, pq.getKey(v) - 2*we);
					}
				} else if (graph.clique(v) == to) {
					n[v].weight_to_current_clique += we;
					if (pq.contains(v)) {
						pq.adjustKey(v, pq.getKey(v) + 2*we);
					}
				} else {
					// TODO could also sum up edge weight changes. if weight changes allow for target cluster to change (multiply by some treshold) --> recompute
					if (++n[v].num_skips > SKIP_THRESHOLD) {
						n[v].num_skips = 0;
						Rating rv = computeBestClique(graph, v);
						if (rv.clique != n[v].desired_target) {
							updateTargetClique(v, rv);
						}
					}
				}

				if (n[v].desired_target == from) {
					n[v].weight_to_target_clique -= we;
					if (pq.contains(v)) {
						pq.adjustKey(v, pq.getKey(v) + 2*we);
					}
				} else if (n[v].desired_target == to) {
					n[v].weight_to_target_clique += we;
					if (pq.contains(v)) {
						pq.adjustKey(v, pq.getKey(v) - 2*we);
					}
				} else {
					if (++n[v].num_skips > SKIP_THRESHOLD) {
						n[v].num_skips = 0;
						Rating rv = computeBestClique(graph, v);
						if (rv.clique != n[v].desired_target) {
							updateTargetClique(v, rv);
						}
					}
				}
			}

			// updates due to clique weight changes
			const NodeWeight wu = graph.nodeWeight(u);
			for (NodeID v : target_cliques[from]) {
				if (pq.contains(v)) {
					pq.adjustKey(v, pq.getKey(v) - wu * graph.nodeWeight(v));
				}
			}
			for (NodeID v : target_cliques[to]) {
				if (pq.contains(v)) {
					pq.adjustKey(v, pq.getKey(v) + wu * graph.nodeWeight(v));
				}
			}

			for (NodeID v : current_cliques[from]) {
				if (pq.contains(v)) {
					pq.adjustKey(v, pq.getKey(v) + wu * graph.nodeWeight(v));
				}
			}
			for (NodeID v : current_cliques[to]) {
				if (pq.contains(v)) {
					pq.adjustKey(v, pq.getKey(v) - wu * graph.nodeWeight(v));
				}
			}

			ASSERT([&] {
				for (size_t i = 0; i < pq.size(); ++i) {
					NodeID u = pq.at(i);
					if (pq.keyAtPos(i) != gain(graph, u, graph.clique(u), n[u].desired_target, n[u].weight_to_current_clique, n[u].weight_to_target_clique)) {
						return false;
					}
				}
				return true;
			}());

		}

		// revert leftovers
		for (Move& m : moves) {
			moveVertex(graph, m.node, m.from);
		}
		current_metric += best_delta;
		ASSERT(current_metric == metrics::edge_deletions(graph));


		DBG << "Pass Nr." << (i + 1) << "improved metric from"
			<< start_metric << "to" << current_metric;
	}

	return current_metric < start_metric;
}

void FMRefiner::moveVertex(Graph& graph, NodeID u, CliqueID to, EdgeWeight weight_to_clique) {
	const CliqueID from = graph.clique(u);
	ASSERT(from != to);

	if (to == ISOLATE_CLIQUE) {
		ASSERT(_empty_cliques.empty());
		to = _empty_cliques.back();
		ASSERT(_clique_weight[to] == 0);
	}
	_clique_weight[from] -= graph.nodeWeight(u);
	const bool from_becomes_empty = _clique_weight[from] == 0;
	const bool to_becomes_non_empty = _clique_weight[to] == 0;
	_clique_weight[to] += graph.nodeWeight(u);
	graph.setClique(u, to);

	removeFromCurrentClique(u, from);
	if (weight_to_clique == std::numeric_limits<EdgeWeight>::max()) {
		weight_to_clique = 0;
		for (const Neighbor& nb : graph.neighbors(u)) {
			if (graph.clique(nb.target) == to) {
				weight_to_clique += graph.edgeWeight(nb.id);
			}
		}
	}
	insertIntoCurrentClique(u, to, weight_to_clique);

	if ( from_becomes_empty ) {
		_empty_cliques.push_back(from);
	}
	if ( to_becomes_non_empty ) {
		ASSERT(_empty_cliques.back() == to);
		_empty_cliques.pop_back();
	}
	++_moved_vertices;
}

FMRefiner::Rating FMRefiner::computeBestClique(Graph& graph, const NodeID u) {
	edge_weight_to_clique.clear();
	const CliqueID from = graph.clique(u);

	for ( const Neighbor& nb : graph.neighbors(u) ) {
		const CliqueID v_c = graph.clique(nb.target);
		edge_weight_to_clique[v_c] += graph.edgeWeight(nb.id);
	}

	const EdgeWeight u_weighted_degree = graph.weightedDegree(u);
	const NodeWeight u_weight = graph.nodeWeight(u);
	const EdgeWeight from_rating = insertions(u_weight, _clique_weight[from] - u_weight, edge_weight_to_clique[from])
									+ deletions(u_weighted_degree, edge_weight_to_clique[from]);

	// ignore the zero gain move of keeping u in its clique
	Rating best_rating = { INVALID_CLIQUE, std::numeric_limits<EdgeWeight>::max(), std::numeric_limits<EdgeWeight>::max(), std::numeric_limits<EdgeWeight>::max() };

	for ( const auto& entry : edge_weight_to_clique ) {
		const CliqueID to = entry.key;
		if (to != from ) {
			const EdgeWeight to_rating = insertions(u_weight, _clique_weight[to], entry.value)
											+ deletions(u_weighted_degree, entry.value);
			if (to_rating < best_rating.rating || (to_rating == best_rating.rating && utils::Randomize::instance().flipCoin())) {
				best_rating = { to, to_rating, to_rating - from_rating, entry.value };
			}
		}
	}

	// Check if it is beneficial to isolate the vertex again. if u is already isolated then this is not triggered
	if ( !_empty_cliques.empty() && u_weighted_degree < best_rating.rating ) {
		best_rating = { ISOLATE_CLIQUE, u_weighted_degree, u_weighted_degree - from_rating, 0 };
	}

	return best_rating;
}

}