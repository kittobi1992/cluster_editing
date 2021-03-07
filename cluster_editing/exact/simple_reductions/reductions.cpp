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

}

int main(int argc,  char **argv) {
  std::cout << "Implementing data reduction rules..." << std::endl;
  return 0;
}
