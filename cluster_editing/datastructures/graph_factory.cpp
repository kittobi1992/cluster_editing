#include "graph_factory.h"

#include "cluster_editing/macros.h"

namespace cluster_editing::ds {

  Graph GraphFactory::construct(const std::vector<size_t>& index_array,
                                const std::vector<NodeID> edges) {
    ASSERT(static_cast<size_t>(index_array.back()) == edges.size());
    Graph graph;
    graph._num_nodes = index_array.size() - 1;
    graph._num_edges = edges.size();
    graph._total_weight = graph._num_nodes; // graph is unweighted

    // Allocate nodes and edges
    graph._nodes.resize(graph._num_nodes);
    graph._nodes.emplace_back(); // Sentinel
    graph._edges.resize(graph._num_edges);

    for ( NodeID i = 0; i < index_array.size(); ++i ) {
      const size_t first_entry = index_array[i];
      graph._nodes[i].setFirstEntry(first_entry);
      graph._nodes[i].setClique(i);
      if ( i < index_array.size() - 1) {
        const NodeID source = static_cast<NodeID>(i);
        const size_t last_entry = index_array[i + 1];
        // Set source of each edge
        for ( size_t j = first_entry; j < last_entry; ++j ) {
          graph._edges[j].setSource(source);
        }
      }
    }

    // Set target of each edge
    for ( size_t i = 0; i < edges.size(); ++i ) {
      graph._edges[i].setTarget(edges[i]);
    }

    return graph;
  }

} // namespace cluster_editing::ds