#include "gmock/gmock.h"

#include "cluster_editing/datastructures/graph_factory.h"

using ::testing::Test;

namespace cluster_editing::ds {

struct WeightedNeigbor {
  EdgeID id;
  NodeID target;
  EdgeWeight weight;
};

class AGraph : public Test {
 public:
  AGraph() :
    graph(GraphFactory::construct(
      { { 1, 2 }, { 0, 2 }, { 0, 1, 3 },
        { 2, 4, 5 }, { 3, 5 }, { 3, 4 } })) { }

  void verifyNeighbors(const Graph& test_graph,
                       const NodeID u,
                       const std::vector<Neighbor>& neighbors) {
    size_t num_neighbors = 0;
    for ( const Neighbor& neighbor : test_graph.neighbors(u) ) {
      const size_t i = num_neighbors;
      ASSERT_EQ(neighbor.id, neighbors[i].id);
      ASSERT_EQ(neighbor.target, neighbors[i].target);
      ++num_neighbors;
    }
    ASSERT_EQ(num_neighbors, neighbors.size());
  }

  void verifyNeighbors(const Graph& test_graph,
                       const NodeID u,
                       const std::vector<WeightedNeigbor>& neighbors) {
    size_t num_neighbors = 0;
    for ( const Neighbor& neighbor : test_graph.neighbors(u) ) {
      const size_t i = num_neighbors;
      ASSERT_EQ(neighbor.id, neighbors[i].id);
      ASSERT_EQ(neighbor.target, neighbors[i].target);
      ASSERT_EQ(1, neighbors[i].weight);
      ++num_neighbors;
    }
    ASSERT_EQ(num_neighbors, neighbors.size());
  }

  void verifyNeighbors(const NodeID u,
                       const std::vector<Neighbor>& neighbors) {
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

TEST_F(AGraph, VerifyEdgeSourceAndTarget1) {
  ASSERT_EQ(0, graph.source(0));
  ASSERT_EQ(1, graph.target(0));
}

TEST_F(AGraph, VerifyEdgeSourceAndTarget2) {
  ASSERT_EQ(0, graph.source(1));
  ASSERT_EQ(2, graph.target(1));
}

TEST_F(AGraph, VerifyEdgeSourceAndTarget3) {
  ASSERT_EQ(3, graph.source(7));
  ASSERT_EQ(2, graph.target(7));
}

TEST_F(AGraph, VerifyEdgeSourceAndTarget4) {
  ASSERT_EQ(5, graph.source(13));
  ASSERT_EQ(4, graph.target(13));
}

TEST_F(AGraph, CheckNeighborsOfNode0) {
  verifyNeighbors(0, { Neighbor { 0, 1 }, Neighbor { 1, 2 } });
}

TEST_F(AGraph, CheckNeighborsOfNode1) {
  verifyNeighbors(1, { Neighbor { 2, 0 }, Neighbor { 3, 2 } });
}

TEST_F(AGraph, CheckNeighborsOfNode2) {
  verifyNeighbors(2, { Neighbor { 4, 0 }, Neighbor { 5, 1 }, Neighbor { 6, 3 } });
}

TEST_F(AGraph, CheckNeighborsOfNode3) {
  verifyNeighbors(3, { Neighbor { 7, 2 }, Neighbor { 8, 4 }, Neighbor { 9, 5 } });
}

TEST_F(AGraph, CheckNeighborsOfNode4) {
  verifyNeighbors(4, { Neighbor { 10, 3 }, Neighbor { 11, 5 } });
}

TEST_F(AGraph, CheckNeighborsOfNode5) {
  verifyNeighbors(5, { Neighbor { 12, 3 }, Neighbor { 13, 4 } });
}

} // namespace cluster_editing::ds