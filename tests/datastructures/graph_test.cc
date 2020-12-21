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
      ASSERT_EQ(test_graph.edgeWeight(neighbor.id), neighbors[i].weight);
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

TEST_F(AGraph, VerifiesThatAllNodeWeightsAreInitialOne) {
  for ( const NodeID& u : graph.nodes() ) {
    ASSERT_EQ(1, graph.nodeWeight(u));
  }
}

TEST_F(AGraph, VerifiesThatAllEdgeWeightsAreInitialOne) {
  for ( const EdgeID& e : graph.edges() ) {
    ASSERT_EQ(1, graph.edgeWeight(e));
  }
}


TEST_F(AGraph, VerifiesThatAllSelfloopWeightsAreInitialZero) {
  for ( const NodeID& u : graph.nodes() ) {
    ASSERT_EQ(0, graph.selfloopWeight(u));
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

TEST_F(AGraph, CheckStatsAfterContraction1) {
  graph.setClique(2, 1);
  graph.setClique(4, 3);
  auto contraction = graph.contract();
  const Graph& c_graph = contraction.first;

  ASSERT_EQ(4, c_graph.numNodes());
  ASSERT_EQ(6, c_graph.numEdges());
  ASSERT_EQ(6, c_graph.totalWeight());
}

TEST_F(AGraph, CheckSelfloopAndNodeWeightAfterContraction1) {
  graph.setClique(2, 1);
  graph.setClique(4, 3);
  auto contraction = graph.contract();
  const Graph& c_graph = contraction.first;

  ASSERT_EQ(0, c_graph.selfloopWeight(0));
  ASSERT_EQ(1, c_graph.nodeWeight(0));
  ASSERT_EQ(1, c_graph.selfloopWeight(1));
  ASSERT_EQ(2, c_graph.nodeWeight(1));
  ASSERT_EQ(1, c_graph.selfloopWeight(2));
  ASSERT_EQ(2, c_graph.nodeWeight(2));
  ASSERT_EQ(0, c_graph.selfloopWeight(3));
  ASSERT_EQ(1, c_graph.nodeWeight(3));
}

TEST_F(AGraph, CheckIncidentEdgesAfterContraction1a) {
  graph.setClique(2, 1);
  graph.setClique(4, 3);
  auto contraction = graph.contract();
  const Graph& c_graph = contraction.first;

  verifyNeighbors(c_graph, 0, { WeightedNeigbor { 0, 1, 2 } });
}

TEST_F(AGraph, CheckIncidentEdgesAfterContraction1b) {
  graph.setClique(2, 1);
  graph.setClique(4, 3);
  auto contraction = graph.contract();
  const Graph& c_graph = contraction.first;

  verifyNeighbors(c_graph, 1, { WeightedNeigbor { 1, 0, 2 }, WeightedNeigbor { 2, 2, 1 } });
}

TEST_F(AGraph, CheckIncidentEdgesAfterContraction1c) {
  graph.setClique(2, 1);
  graph.setClique(4, 3);
  auto contraction = graph.contract();
  const Graph& c_graph = contraction.first;

  verifyNeighbors(c_graph, 2, { WeightedNeigbor { 3, 1, 1 }, WeightedNeigbor { 4, 3, 2 } });
}

TEST_F(AGraph, CheckIncidentEdgesAfterContraction1d) {
  graph.setClique(2, 1);
  graph.setClique(4, 3);
  auto contraction = graph.contract();
  const Graph& c_graph = contraction.first;

  verifyNeighbors(c_graph, 3, { WeightedNeigbor { 5, 2, 2 } });
}

TEST_F(AGraph, CheckStatsAfterContraction2) {
  graph.setClique(1, 0);
  graph.setClique(2, 0);
  graph.setClique(4, 3);
  graph.setClique(5, 3);
  auto contraction = graph.contract();
  const Graph& c_graph = contraction.first;

  ASSERT_EQ(2, c_graph.numNodes());
  ASSERT_EQ(2, c_graph.numEdges());
  ASSERT_EQ(6, c_graph.totalWeight());
}

TEST_F(AGraph, CheckSelfloopAndNodeWeightAfterContraction2) {
  graph.setClique(1, 0);
  graph.setClique(2, 0);
  graph.setClique(4, 3);
  graph.setClique(5, 3);
  auto contraction = graph.contract();
  const Graph& c_graph = contraction.first;

  ASSERT_EQ(3, c_graph.selfloopWeight(0));
  ASSERT_EQ(3, c_graph.nodeWeight(0));
  ASSERT_EQ(3, c_graph.selfloopWeight(1));
  ASSERT_EQ(3, c_graph.nodeWeight(1));
}

TEST_F(AGraph, CheckIncidentEdgesAfterContraction2a) {
  graph.setClique(1, 0);
  graph.setClique(2, 0);
  graph.setClique(4, 3);
  graph.setClique(5, 3);
  auto contraction = graph.contract();
  const Graph& c_graph = contraction.first;

  verifyNeighbors(c_graph, 0, { WeightedNeigbor { 0, 1, 1 } });
}

TEST_F(AGraph, CheckIncidentEdgesAfterContraction2b) {
  graph.setClique(1, 0);
  graph.setClique(2, 0);
  graph.setClique(4, 3);
  graph.setClique(5, 3);
  auto contraction = graph.contract();
  const Graph& c_graph = contraction.first;

  verifyNeighbors(c_graph, 1, { WeightedNeigbor { 1, 0, 1 } });
}

TEST_F(AGraph, CheckStatsAfterContraction3) {
  Graph t_graph = GraphFactory::construct(
    { { 1, 2 }, { 0, 3 }, { 0, 4 }, { 1, 4 }, { 2, 3, 5 },
      { 4, 6, 7 }, { 5, 8 }, { 5, 9 }, { 6, 9 }, { 7, 8 } });
  t_graph.setClique(1, 0);
  t_graph.setClique(4, 2);
  t_graph.setClique(6, 5);
  t_graph.setClique(9, 8);
  auto contraction = t_graph.contract();
  const Graph& c_graph = contraction.first;

  ASSERT_EQ(6, c_graph.numNodes());
  ASSERT_EQ(14, c_graph.numEdges());
  ASSERT_EQ(10, c_graph.totalWeight());
}

TEST_F(AGraph, CheckSelfloopAndNodeWeightAfterContraction3) {
  Graph t_graph = GraphFactory::construct(
    { { 1, 2 }, { 0, 3 }, { 0, 4 }, { 1, 4 }, { 2, 3, 5 },
      { 4, 6, 7 }, { 5, 8 }, { 5, 9 }, { 6, 9 }, { 7, 8 } });
  t_graph.setClique(1, 0);
  t_graph.setClique(4, 2);
  t_graph.setClique(6, 5);
  t_graph.setClique(9, 8);
  auto contraction = t_graph.contract();
  const Graph& c_graph = contraction.first;

  ASSERT_EQ(1, c_graph.selfloopWeight(0));
  ASSERT_EQ(2, c_graph.nodeWeight(0));
  ASSERT_EQ(1, c_graph.selfloopWeight(1));
  ASSERT_EQ(2, c_graph.nodeWeight(1));
  ASSERT_EQ(0, c_graph.selfloopWeight(2));
  ASSERT_EQ(1, c_graph.nodeWeight(2));
  ASSERT_EQ(1, c_graph.selfloopWeight(3));
  ASSERT_EQ(2, c_graph.nodeWeight(3));
  ASSERT_EQ(0, c_graph.selfloopWeight(4));
  ASSERT_EQ(1, c_graph.nodeWeight(4));
  ASSERT_EQ(1, c_graph.selfloopWeight(5));
  ASSERT_EQ(2, c_graph.nodeWeight(5));
}

TEST_F(AGraph, CheckIncidentEdgesAfterContraction3a) {
  Graph t_graph = GraphFactory::construct(
    { { 1, 2 }, { 0, 3 }, { 0, 4 }, { 1, 4 }, { 2, 3, 5 },
      { 4, 6, 7 }, { 5, 8 }, { 5, 9 }, { 6, 9 }, { 7, 8 } });
  t_graph.setClique(1, 0);
  t_graph.setClique(4, 2);
  t_graph.setClique(6, 5);
  t_graph.setClique(9, 8);
  auto contraction = t_graph.contract();
  const Graph& c_graph = contraction.first;

  verifyNeighbors(c_graph, 0, { WeightedNeigbor { 0, 1, 1 }, WeightedNeigbor { 1, 2, 1 } });
}

TEST_F(AGraph, CheckIncidentEdgesAfterContraction3b) {
  Graph t_graph = GraphFactory::construct(
    { { 1, 2 }, { 0, 3 }, { 0, 4 }, { 1, 4 }, { 2, 3, 5 },
      { 4, 6, 7 }, { 5, 8 }, { 5, 9 }, { 6, 9 }, { 7, 8 } });
  t_graph.setClique(1, 0);
  t_graph.setClique(4, 2);
  t_graph.setClique(6, 5);
  t_graph.setClique(9, 8);
  auto contraction = t_graph.contract();
  const Graph& c_graph = contraction.first;

  verifyNeighbors(c_graph, 1, { WeightedNeigbor { 2, 0, 1 }, WeightedNeigbor { 3, 2, 1 }, WeightedNeigbor { 4, 3, 1 } });
}

TEST_F(AGraph, CheckIncidentEdgesAfterContraction3c) {
  Graph t_graph = GraphFactory::construct(
    { { 1, 2 }, { 0, 3 }, { 0, 4 }, { 1, 4 }, { 2, 3, 5 },
      { 4, 6, 7 }, { 5, 8 }, { 5, 9 }, { 6, 9 }, { 7, 8 } });
  t_graph.setClique(1, 0);
  t_graph.setClique(4, 2);
  t_graph.setClique(6, 5);
  t_graph.setClique(9, 8);
  auto contraction = t_graph.contract();
  const Graph& c_graph = contraction.first;

  verifyNeighbors(c_graph, 2, { WeightedNeigbor { 5, 0, 1 }, WeightedNeigbor { 6, 1, 1 } });
}

TEST_F(AGraph, CheckIncidentEdgesAfterContraction3d) {
  Graph t_graph = GraphFactory::construct(
    { { 1, 2 }, { 0, 3 }, { 0, 4 }, { 1, 4 }, { 2, 3, 5 },
      { 4, 6, 7 }, { 5, 8 }, { 5, 9 }, { 6, 9 }, { 7, 8 } });
  t_graph.setClique(1, 0);
  t_graph.setClique(4, 2);
  t_graph.setClique(6, 5);
  t_graph.setClique(9, 8);
  auto contraction = t_graph.contract();
  const Graph& c_graph = contraction.first;

  verifyNeighbors(c_graph, 3, { WeightedNeigbor { 7, 1, 1 }, WeightedNeigbor { 8, 4, 1 }, WeightedNeigbor { 9, 5, 1 } });
}

TEST_F(AGraph, CheckIncidentEdgesAfterContraction3e) {
  Graph t_graph = GraphFactory::construct(
    { { 1, 2 }, { 0, 3 }, { 0, 4 }, { 1, 4 }, { 2, 3, 5 },
      { 4, 6, 7 }, { 5, 8 }, { 5, 9 }, { 6, 9 }, { 7, 8 } });
  t_graph.setClique(1, 0);
  t_graph.setClique(4, 2);
  t_graph.setClique(6, 5);
  t_graph.setClique(9, 8);
  auto contraction = t_graph.contract();
  const Graph& c_graph = contraction.first;

  verifyNeighbors(c_graph, 4, { WeightedNeigbor { 10, 3, 1 }, WeightedNeigbor { 11, 5, 1 } });
}

TEST_F(AGraph, CheckIncidentEdgesAfterContraction3f) {
  Graph t_graph = GraphFactory::construct(
    { { 1, 2 }, { 0, 3 }, { 0, 4 }, { 1, 4 }, { 2, 3, 5 },
      { 4, 6, 7 }, { 5, 8 }, { 5, 9 }, { 6, 9 }, { 7, 8 } });
  t_graph.setClique(1, 0);
  t_graph.setClique(4, 2);
  t_graph.setClique(6, 5);
  t_graph.setClique(9, 8);
  auto contraction = t_graph.contract();
  const Graph& c_graph = contraction.first;

  verifyNeighbors(c_graph, 5, { WeightedNeigbor { 12, 3, 1 }, WeightedNeigbor { 13, 4, 1 } });
}

TEST_F(AGraph, CheckStatsAfterTwoContractions) {
  Graph t_graph = GraphFactory::construct(
    { { 1, 2 }, { 0, 3 }, { 0, 4 }, { 1, 4 }, { 2, 3, 5 },
      { 4, 6, 7 }, { 5, 8 }, { 5, 9 }, { 6, 9 }, { 7, 8 } });
  t_graph.setClique(1, 0);
  t_graph.setClique(4, 2);
  t_graph.setClique(6, 5);
  t_graph.setClique(9, 8);
  auto contraction = t_graph.contract();
  Graph c_graph = std::move(contraction.first);
  c_graph.setClique(2, 1);
  c_graph.setClique(4, 3);
  auto contraction_2 = c_graph.contract();
  c_graph = std::move(contraction_2.first);

  ASSERT_EQ(4, c_graph.numNodes());
  ASSERT_EQ(6, c_graph.numEdges());
  ASSERT_EQ(10, c_graph.totalWeight());
}

TEST_F(AGraph, CheckSelfloopAndNodeWeightAfterTwoContractions) {
  Graph t_graph = GraphFactory::construct(
    { { 1, 2 }, { 0, 3 }, { 0, 4 }, { 1, 4 }, { 2, 3, 5 },
      { 4, 6, 7 }, { 5, 8 }, { 5, 9 }, { 6, 9 }, { 7, 8 } });
  t_graph.setClique(1, 0);
  t_graph.setClique(4, 2);
  t_graph.setClique(6, 5);
  t_graph.setClique(9, 8);
  auto contraction = t_graph.contract();
  Graph c_graph = std::move(contraction.first);
  c_graph.setClique(2, 1);
  c_graph.setClique(4, 3);
  auto contraction_2 = c_graph.contract();
  c_graph = std::move(contraction_2.first);

  ASSERT_EQ(1, c_graph.selfloopWeight(0));
  ASSERT_EQ(2, c_graph.nodeWeight(0));
  ASSERT_EQ(2, c_graph.selfloopWeight(1));
  ASSERT_EQ(3, c_graph.nodeWeight(1));
  ASSERT_EQ(2, c_graph.selfloopWeight(2));
  ASSERT_EQ(3, c_graph.nodeWeight(2));
  ASSERT_EQ(1, c_graph.selfloopWeight(3));
  ASSERT_EQ(2, c_graph.nodeWeight(3));
}

TEST_F(AGraph, CheckIncidentEdgesAfterTwoContractionsA) {
  Graph t_graph = GraphFactory::construct(
    { { 1, 2 }, { 0, 3 }, { 0, 4 }, { 1, 4 }, { 2, 3, 5 },
      { 4, 6, 7 }, { 5, 8 }, { 5, 9 }, { 6, 9 }, { 7, 8 } });
  t_graph.setClique(1, 0);
  t_graph.setClique(4, 2);
  t_graph.setClique(6, 5);
  t_graph.setClique(9, 8);
  auto contraction = t_graph.contract();
  Graph c_graph = std::move(contraction.first);
  c_graph.setClique(2, 1);
  c_graph.setClique(4, 3);
  auto contraction_2 = c_graph.contract();
  c_graph = std::move(contraction_2.first);

  verifyNeighbors(c_graph, 0, { WeightedNeigbor { 0, 1, 2 } });
}

TEST_F(AGraph, CheckIncidentEdgesAfterTwoContractionsB) {
  Graph t_graph = GraphFactory::construct(
    { { 1, 2 }, { 0, 3 }, { 0, 4 }, { 1, 4 }, { 2, 3, 5 },
      { 4, 6, 7 }, { 5, 8 }, { 5, 9 }, { 6, 9 }, { 7, 8 } });
  t_graph.setClique(1, 0);
  t_graph.setClique(4, 2);
  t_graph.setClique(6, 5);
  t_graph.setClique(9, 8);
  auto contraction = t_graph.contract();
  Graph c_graph = std::move(contraction.first);
  c_graph.setClique(2, 1);
  c_graph.setClique(4, 3);
  auto contraction_2 = c_graph.contract();
  c_graph = std::move(contraction_2.first);

  verifyNeighbors(c_graph, 1, { WeightedNeigbor { 1, 0, 2 }, WeightedNeigbor { 2, 2, 1 } });
}

TEST_F(AGraph, CheckIncidentEdgesAfterTwoContractionsC) {
  Graph t_graph = GraphFactory::construct(
    { { 1, 2 }, { 0, 3 }, { 0, 4 }, { 1, 4 }, { 2, 3, 5 },
      { 4, 6, 7 }, { 5, 8 }, { 5, 9 }, { 6, 9 }, { 7, 8 } });
  t_graph.setClique(1, 0);
  t_graph.setClique(4, 2);
  t_graph.setClique(6, 5);
  t_graph.setClique(9, 8);
  auto contraction = t_graph.contract();
  Graph c_graph = std::move(contraction.first);
  c_graph.setClique(2, 1);
  c_graph.setClique(4, 3);
  auto contraction_2 = c_graph.contract();
  c_graph = std::move(contraction_2.first);

  verifyNeighbors(c_graph, 2, { WeightedNeigbor { 3, 1, 1 }, WeightedNeigbor { 4, 3, 2 } });
}

TEST_F(AGraph, CheckIncidentEdgesAfterTwoContractionsD) {
  Graph t_graph = GraphFactory::construct(
    { { 1, 2 }, { 0, 3 }, { 0, 4 }, { 1, 4 }, { 2, 3, 5 },
      { 4, 6, 7 }, { 5, 8 }, { 5, 9 }, { 6, 9 }, { 7, 8 } });
  t_graph.setClique(1, 0);
  t_graph.setClique(4, 2);
  t_graph.setClique(6, 5);
  t_graph.setClique(9, 8);
  auto contraction = t_graph.contract();
  Graph c_graph = std::move(contraction.first);
  c_graph.setClique(2, 1);
  c_graph.setClique(4, 3);
  auto contraction_2 = c_graph.contract();
  c_graph = std::move(contraction_2.first);

  verifyNeighbors(c_graph, 3, { WeightedNeigbor { 5, 2, 2 } });
}


} // namespace cluster_editing::ds