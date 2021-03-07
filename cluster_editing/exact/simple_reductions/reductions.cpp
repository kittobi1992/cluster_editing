#include <iostream>
#include <vector>
#include <limits>

using namespace std;

double INF = numeric_limits<double>::infinity();

/*
  Rule 1: heavy non edge rule
    Set an edge to FORBIDDEN (weight -INF) if the degree of the node is greadter \
    than the weight between it and the other side of the edge in question

  @params:  graph - reference to the adjacency matrix of the graph
            u - first vertex in the edge
            v - second vertex in the edge
  @return:  none
*/
void heavy_non_edge_rule(vector<vector<int> > &graph, int u, int v) {
  int sum = 0;
  for (int w = 0; w < graph[u].size(); w++) {
    sum += graph[u][w];
  }

  if (abs(graph[u][v]) >= sum) {
    cout << "FORBIDDEN: (" << u << ", " << v << ")" << endl;
    graph[u][v] = -INF;
  }
}


/*
  Rule 2: heavy non edge single end rule
    MERGE vertices u and v if the weight of the edge is >= to the sum of abs of the
    edges between u and every node in V\{u, v}

  @params:  graph - reference to the adjacency matrix of the graph
            u - first vertex in the edge
            v - second vertex in the edge
  @return:  none
*/
void heavy_edge_single_end_rule(vector<vector<int> > &graph, int u, int v) {
  int sum = 0;
  for (int w = 0; w < graph.size(); w++) {
    // Don't count either uu or uv
    if (w == u or w == v)
      continue;

    sum += abs(graph[u][w]);
  }

  if (graph[u][v] >= sum) {
    cout << "MERGE: (" << u << ", " << v << ")" << endl;
  }
}


/*
  Rule 3: heavy edge both end rule
    MERGE vertices u and v if the weight of edge (u, v) is >= to the sum of all
    edge weights between u and N(u)\{v} added to the sum of all edge weights
    between v and N(v)\{u}


  @params:  graph - reference to the adjacency matrix of the graph
            u - first vertex in the edge
            v - second vertex in the edge
  @return:  none
*/
void heavy_edge_both_end_rule(vector<vector<int> > &graph, int u, int v) {
  int sum_u = 0;
  int sum_v = 0
  for (int w = 0; w < graph.size(); w++) {
    if (w != v)
      sum_u += graph[u][w];
    if (w != u)
      sum_v += graph[v][w];
  }

  if (graph[u][v] >= sum_u + sum_v) {
    cout << "MERGE: (" << u << ", " << v << ")" << endl;
  }
}

int main(int argc,  char **argv) {
  std::cout << "Implementing data reduction rules..." << std::endl;
  return 0;
}
