#include "graph_io.h"

#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "cluster_editing/macros.h"
#include "cluster_editing/datastructures/graph_factory.h"

namespace cluster_editing {
namespace io {

using AdjacencyList = std::vector<std::vector<NodeID>>;

void readHeader(NodeID& num_nodes,
                EdgeID& num_edges) {
  std::string skip;
  std::cin >> skip >> skip >> num_nodes >> num_edges;
}

void readEdges(const EdgeID num_edges,
               AdjacencyList& adj_list) {
  for ( EdgeID i = 0; i < num_edges; ++i ) {
    NodeID u, v;
    std::cin >> u >> v;
    --u; --v;
    adj_list[u].push_back(v);
    adj_list[v].push_back(u);
  }
}

Graph readGraphFile() {
  AdjacencyList adj_list;
  NodeID num_nodes = 0;
  EdgeID num_edges = 0;
  readHeader(num_nodes, num_edges);

  // Read Egges
  adj_list.resize(num_nodes);
  readEdges(num_edges, adj_list);

  return ds::GraphFactory::construct(adj_list);
}
}
} // namespace cluster_editing::io