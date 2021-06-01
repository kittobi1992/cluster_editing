/*******************************************************************************
 * This file is part of KaPoCE.
 *
 * KaPoCE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaPoCE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaPoCE.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "graph_io.h"

#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "cluster_editing/macros.h"
#include "cluster_editing/datastructures/graph_factory.h"

namespace cluster_editing::io {

using AdjacencyList = std::vector<std::vector<NodeID>>;

namespace {

void readHeader(std::istream& file,
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

void readEdges(std::istream& file,
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
} // namespace

Graph readGraphFile() {
  AdjacencyList adj_list;
  NodeID num_nodes = 0;
  EdgeID num_edges = 0;
  readHeader(std::cin, num_nodes, num_edges);

  // Read Egges
  adj_list.resize(num_nodes);
  readEdges(std::cin, num_edges, adj_list);

  return ds::GraphFactory::construct(adj_list);
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

} // namespace cluster_editing::io