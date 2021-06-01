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

#include "gmock/gmock.h"

#include "cluster_editing/datastructures/graph_factory.h"

using ::testing::Test;

namespace cluster_editing::ds {

class AGraph : public Test {
 public:
  AGraph() :
    graph(GraphFactory::construct(
      { { 1, 2 }, { 0, 2 }, { 0, 1, 3 },
        { 2, 4, 5 }, { 3, 5 }, { 3, 4 } })) { }

  void verifyNeighbors(const Graph& test_graph,
                       const NodeID u,
                       const std::vector<NodeID>& neighbors) {
    size_t num_neighbors = 0;
    for ( const NodeID& v : test_graph.neighbors(u) ) {
      const size_t i = num_neighbors;
      ASSERT_EQ(v, neighbors[i]);
      ++num_neighbors;
    }
    ASSERT_EQ(num_neighbors, neighbors.size());
  }

  void verifyNeighbors(const NodeID u,
                       const std::vector<NodeID>& neighbors) {
    verifyNeighbors(graph, u, neighbors);
  }

  Graph graph;
};

TEST_F(AGraph, CheckStats) {
  ASSERT_EQ(6, graph.numNodes());
  ASSERT_EQ(14, graph.numEdges());
  ASSERT_EQ(6, graph.totalWeight());
}

TEST_F(AGraph, VerifiesThatEachNodeIsInitialInItsOwnClique) {
  for ( const NodeID& u : graph.nodes() ) {
    ASSERT_EQ(u, graph.clique(u));
  }
}

TEST_F(AGraph, VerifyNodeDegrees) {
  ASSERT_EQ(2, graph.degree(0));
  ASSERT_EQ(2, graph.degree(1));
  ASSERT_EQ(3, graph.degree(2));
  ASSERT_EQ(3, graph.degree(3));
  ASSERT_EQ(2, graph.degree(4));
  ASSERT_EQ(2, graph.degree(5));
}

TEST_F(AGraph, CheckNeighborsOfNode0) {
  verifyNeighbors(0, { 1, 2 });
}

TEST_F(AGraph, CheckNeighborsOfNode1) {
  verifyNeighbors(1, { 0, 2 });
}

TEST_F(AGraph, CheckNeighborsOfNode2) {
  verifyNeighbors(2, { 0, 1, 3 });
}

TEST_F(AGraph, CheckNeighborsOfNode3) {
  verifyNeighbors(3, { 2, 4, 5 });
}

TEST_F(AGraph, CheckNeighborsOfNode4) {
  verifyNeighbors(4, { 3, 5 });
}

TEST_F(AGraph, CheckNeighborsOfNode5) {
  verifyNeighbors(5, { 3, 4 });
}

} // namespace cluster_editing::ds