#pragma once

#include "cluster_editing/refinement/i_refiner.h"
#include "cluster_editing/context/context.h"
#include "cluster_editing/datastructures/sparse_map.h"

#include "cluster_editing/datastructures/pq.h"

namespace cluster_editing {

class FMRefiner final : public IRefiner {
private:
	static constexpr bool debug = false;

	struct Rating {
		CliqueID clique;
		EdgeWeight rating;
		EdgeWeight delta;
	};

	struct NodeData {
		CliqueID desired_target = INVALID_CLIQUE;
		NodeID ttv_index = INVALID_NODE;
		uint32_t num_skips = 0,  last_move_round = 0;
	};

	struct Move {
		NodeID node;
		CliqueID from, to;
	};

	CliqueID ISOLATE_CLIQUE = INVALID_CLIQUE - 1;

public:
	explicit FMRefiner(const Graph& graph, const Context& context) :
			_context(context),
			_nodes(),
			_clique_weight(graph.numNodes()),
			_empty_cliques(),
			_rating(graph.numNodes()),
			n(graph.numNodes()),
			target_clique_to_nodes(graph.numNodes()),
			pq(graph.numNodes())
			{ }

	size_t movedVertices() const {
		return _moved_vertices;
	}

private:

	void initializeImpl(Graph& graph) final;

	bool refineImpl(Graph& graph) final ;

	void moveVertex(Graph& graph, const NodeID u, const CliqueID to);

	Rating computeBestClique(Graph& graph, const NodeID u);

	void updateLookupDatastructure(NodeID u, CliqueID target) {
		removeFromLookupDataStructure(u);
		insertIntoLookupDatastructure(u, target);
	}

	void removeFromLookupDataStructure(NodeID u) {
		CliqueID old_target = n[u].desired_target;
		if (old_target != ISOLATE_CLIQUE) {
			target_clique_to_nodes[old_target][n[u].ttv_index] = target_clique_to_nodes[old_target].back();
			n[target_clique_to_nodes[old_target].back()].ttv_index = n[u].ttv_index;
			target_clique_to_nodes[old_target].pop_back();
			n[u].ttv_index = std::numeric_limits<uint32_t>::max();
		}
	}

	void insertIntoLookupDatastructure(NodeID u, CliqueID target) {
		n[u].desired_target = target;
		if (target != ISOLATE_CLIQUE) {
			n[u].ttv_index = target_clique_to_nodes[target].size();
			target_clique_to_nodes[target].push_back(u);
		}
	}

	const Context& _context;
	size_t _moved_vertices;
	std::vector<NodeID> _nodes;
	std::vector<NodeWeight> _clique_weight;
	std::vector<CliqueID> _empty_cliques;
	ds::SparseMap<CliqueID, EdgeWeight> _rating;

	uint32_t move_round;
	std::vector<NodeData> n;
	std::vector<std::vector<NodeID>> target_clique_to_nodes;
	mt_kahypar::ds::MinHeap<EdgeWeight, NodeID> pq;
	std::vector<Move> moves;
};
}  // namespace cluster_editing
