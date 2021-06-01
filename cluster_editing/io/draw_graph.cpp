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

#include "cluster_editing/io/draw_graph.h"
#include "cluster_editing/io/graph_io.h"
#include "cluster_editing/io/ipe.h"

#include <algorithm>
#include <fstream>

void ignore_line(std::istream& in) {
  in.ignore(std::numeric_limits<std::streamsize>::max(), in.widen('\n'));
}

void draw_graph(const std::filesystem::path& graph_file_in,
                const std::filesystem::path& coords_file_in,
                const std::filesystem::path& partition_file_in,
                const std::filesystem::path& ipe_file_out)
{
  Graph G = io::readGraphFile(graph_file_in);
  draw_graph(G, read_coords(coords_file_in, G.numNodes()),
             read_partition(partition_file_in, G.numNodes()),
             ipe_file_out);
}

void draw_graph(const Graph& G,
                const std::vector<std::vector<double>>& coords,
                const std::vector<unsigned>& partition,
                const std::filesystem::path& ipe_file_out)
{
  unsigned max_partition = *std::max_element(partition.begin(), partition.end());
  max_partition++;

  // write ipe file
  IpeFile ipe(ipe_file_out);
  for (auto u : G.nodes()) {
    for (const NodeID& v : G.neighbors(u)) {
      Color v_col = hsv(partition[v] * 360 / max_partition, 1, 0.8);

      if (partition[u] == partition[v]) {
        ipe.line(coords[v][0], coords[v][1], coords[u][0], coords[u][1], v_col,
                 1.2, true);
      } else {
        ipe.line(coords[v][0], coords[v][1], coords[u][0], coords[u][1],
                 rgb(0.5, 0.5, 0.5), 0.4, true);
      }

      ipe.point(coords[v][0], coords[v][1], v_col);
    }
    Color u_col = hsv(partition[u] * 360 / max_partition, 1, 0.8);
    ipe.point(coords[u][0], coords[u][1], u_col);
  }
}

std::vector<std::vector<double>> read_coords(
    const std::filesystem::path& coords_file_in, unsigned num_nodes)
{
  std::vector<std::vector<double>> coords(num_nodes, {0, 0});
  std::ifstream in(coords_file_in);
  int offset = 1;
  while (in.peek() == '%') {
    ignore_line(in);
  }
  int v;
  double x, y;
  while (in >> v && in >> x && in >> y) {
    unsigned int id = v - offset;
    if (id >= coords.size()) {
      coords.resize(id + 1, {0.0, 0.0});
    }
    coords[id][0] = x;
    coords[id][1] = y;
  }
  return coords;
}

std::vector<unsigned> read_partition(
    const std::filesystem::path& partition_file_in, unsigned num_nodes)
{
  std::vector<unsigned> partition(num_nodes, 0);
  std::ifstream in(partition_file_in);
  int offset = 1;
  int v, p;
  while (in >> v && in >> p) {
    partition[v - offset] = p;
  }
  return partition;
}
