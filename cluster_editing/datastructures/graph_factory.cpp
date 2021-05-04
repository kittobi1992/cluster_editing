#include "graph_factory.h"

#include "cluster_editing/macros.h"

namespace cluster_editing {
namespace ds {

  using AdjacencyList = std::vector<std::vector<NodeID>>;

  Graph GraphFactory::construct(const AdjacencyList& adj_list) {
    Graph graph;
    graph._num_nodes = adj_list.size();
    graph._total_weight = graph._num_nodes; // graph is unweighted

    graph._nodes.resize(graph._num_nodes);
    graph._best_cliques.resize(graph._num_nodes);
    size_t current_idx = 0;
    for ( NodeID source = 0; source < graph._num_nodes; ++source ) {
      graph._nodes[source].setFirstEntry(current_idx);
      graph._nodes[source].setClique(source);
      graph._best_cliques[source] = source;
      graph._max_degree = std::max(
        graph._max_degree, static_cast<NodeID>(adj_list[source].size()));
      for ( const NodeID& target : adj_list[source] ) {
        graph._edges.emplace_back(target);
        ++current_idx;
      }
    }
    graph._num_edges = static_cast<EdgeID>(graph._edges.size());
    graph._best_quality = graph._num_edges;

    // Sentinel
    graph._nodes.emplace_back();
    graph._nodes.back().setFirstEntry(current_idx);

    return graph;
  }
}
} // namespace cluster_editing::ds