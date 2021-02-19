#include "fm_refiner.h"


#include "cluster_editing/utils/randomize.h"
#include "cluster_editing/metrics.h"
#include "cluster_editing/utils/timer.h"

namespace cluster_editing {



namespace {
ATTRIBUTE_ALWAYS_INLINE EdgeWeight insertions(const NodeWeight u_weight,
											  const NodeWeight clique_weight,
											  const EdgeWeight incident_edge_weight) {
	return u_weight * clique_weight - incident_edge_weight;
}
}

void FMRefiner::initializeImpl(Graph& graph) {
	_moved_vertices = 0;
	_clique_weight.assign(graph.numNodes(), 0);
	_empty_cliques.clear();
	_nodes.clear();
	for ( const NodeID& u : graph.nodes() ) {
		_nodes.push_back(u);
		_clique_weight[graph.clique(u)] += graph.nodeWeight(u);
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
		for ( const NodeID& u : _nodes ) {
			Rating rating = computeBestClique(graph, u);
			pq.insert(u, rating.delta);
			insertIntoLookupDatastructure(u, rating.clique);
		}

		// perform moves
		while (!pq.empty()) {
			EdgeWeight estimated_gain = pq.topKey();
			NodeID u = pq.top();

			Rating rating = computeBestClique(graph, u);
			if (rating.delta != estimated_gain) {
				// retry! 		or do KaHiP approach and just apply? only if improvement?
				pq.adjustKey(u, rating.delta);
				updateLookupDatastructure(u, rating.clique);
				continue;
			}

			round_delta += rating.delta;
			if (round_delta < best_delta) {
				// permanently apply all moves
				moves.clear();
			} else {
				// store for rollback later
				moves.push_back({ u, graph.clique(u), rating.clique });
			}

			removeFromLookupDataStructure(u);
			moveVertex(graph, u, rating.clique);
			ASSERT(current_metric + round_delta == metrics::edits(graph), "Rating is wrong. Expected:" << metrics::edits(graph) << "but is" << (current_metric + round_delta));

			// update neighbors
			

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

void FMRefiner::moveVertex(Graph& graph, const NodeID u, const CliqueID to) {
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
	++_moved_vertices;
}

FMRefiner::Rating FMRefiner::computeBestClique(Graph& graph, const NodeID u) {
	_rating.clear();
	const CliqueID u_c = graph.clique(u);
	_rating[u_c] = 0;

	for ( const Neighbor& n : graph.neighbors(u) ) {
		const CliqueID v_c = graph.clique(n.target);
		_rating[v_c] += graph.edgeWeight(n.id);
	}

	const EdgeWeight u_weighted_degree = graph.weightedDegree(u);
	const NodeWeight u_weight = graph.nodeWeight(u);
	const EdgeWeight current_rating = insertions(u_weight, _clique_weight[u_c] -
														   u_weight, _rating[u_c]) + (u_weighted_degree - _rating[u_c]) /* deletions */;

	// ignore the zero gain move of keeping u in its clique
	Rating best_rating = { INVALID_CLIQUE, std::numeric_limits<EdgeWeight>::max(), std::numeric_limits<EdgeWeight>::max() };

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

	// Check if it is beneficial to isolate the vertex again
	if ( !_empty_cliques.empty() && u_weighted_degree < best_rating.rating ) {
		best_rating.clique = ISOLATE_CLIQUE;
		best_rating.rating = u_weighted_degree;
		best_rating.delta = u_weighted_degree - current_rating;
	}

	return best_rating;
}

}