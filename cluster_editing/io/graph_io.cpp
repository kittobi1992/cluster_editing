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

void readHeader(std::ifstream& file,
                NodeID& num_nodes,
                EdgeID& num_edges) {
  std::string line;
  std::getline(file, line);

  // Skip comments
  while ( line[0] == 'c' )  {
    std::getline(file, line);
  }

  ASSERT(line[0] == 'p');
  std::istringstream sstream(line);
  std::string skip;
  sstream >> skip >> skip >> num_nodes >> num_edges;
}

void readEdges(std::ifstream& file,
               const EdgeID num_edges,
               AdjacencyList& adj_list) {
  std::string line;
  for ( EdgeID i = 0; i < num_edges; ++i ) {
    std::getline(file, line);
    while (line[0] == 'c') {
      std::getline(file, line);
    }

    std::istringstream line_stream(line);
    if (line_stream.peek() == EOF) {
      ERROR("Line should contain an edge");
    }

    NodeID u, v;
    line_stream >> u >> v;
    --u; --v;
    adj_list[u].push_back(v);
    adj_list[v].push_back(u);
  }
}

Graph readGraphFile(const std::string& filename) {
  ASSERT(!filename.empty(), "No filename for graph file specified");
  std::ifstream file(filename);
  AdjacencyList adj_list;
  if ( file ) {
    NodeID num_nodes = 0;
    EdgeID num_edges = 0;
    readHeader(file, num_nodes, num_edges);

    // Read Egges
    adj_list.resize(num_nodes);
    readEdges(file, num_edges, adj_list);

    file.close();
  } else {
    ERROR("File" << filename << "does not exist!");
  }

  return ds::GraphFactory::construct(adj_list);
}
}
} // namespace cluster_editing::io