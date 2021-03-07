#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <limits>

using namespace std;

int n = 0;

double INF = numeric_limits<double>::infinity();

void print(int **Adj, size_t size) {
  // Traverse the Adj[][]
  for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {

          // Print the value at Adj[i][j]
          printf("%d ", Adj[i][j]);
      }
      printf("\n");
  }
}

/*
  Rule 1: heavy non edge rule
    Set an edge to FORBIDDEN (weight -INF) if the degree of the node is greadter \
    than the weight between it and the other side of the edge in question

  @params:  graph - reference to the adjacency matrix of the graph
            u - first vertex in the edge
            v - second vertex in the edge
  @return:  none
*/
void heavy_non_edge_rule(int **graph, size_t size, int u, int v) {
  int sum = 0;
  for (int w = 0; w < size; w++) {
    sum += graph[u][w];
  }

  if (abs(graph[u][v]) >= sum) {
    cout << "FORBIDDEN: (" << u+1 << ", " << v+1 << ")" << endl;
    // graph[u][v] = -INF;
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
void heavy_edge_single_end_rule(int **graph, size_t size, int u, int v) {
  int sum = 0;
  for (int w = 0; w < size; w++) {
    // Don't count either uu or uv
    if (w == u or w == v)
      continue;

    sum += abs(graph[u][w]);
  }

  if (graph[u][v] >= sum) {
    cout << "MERGE: (" << u+1 << ", " << v+1 << ")" << endl;
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
void heavy_edge_both_end_rule(int **graph, size_t size, int u, int v) {
  int sum_u = 0;
  int sum_v = 0;
  for (int w = 0; w < size; w++) {
    if (w != v)
      sum_u += graph[u][w];
    if (w != u)
      sum_v += graph[v][w];
  }

  if (graph[u][v] >= sum_u + sum_v) {
    cout << "MERGE: (" << u+1 << ", " << v+1 << ")" << endl;
  }
}



int** makeAdjacencyMatrix() {
  string line;
  ifstream inputfile("graphs/exact001.gr");
  int** Adj = 0;

  if (inputfile.is_open()) {
    getline(inputfile, line);

    string edges = line.substr(line.find_last_of(" ") + 1);
    string rest = line.substr(0,line.find_last_of(" "));
    string vertices = rest.substr(rest.find_last_of(" ") + 1);

    n = stoi(vertices);
    int e = stoi(edges);

    // cout << n << e << endl;

    // int Adj[n+1][n+1];
    Adj = new int*[n];

    for (int i = 0; i < n; i++) {
      Adj[i] = new int[n];
      for (int j = 0; j < n; j++) {
        Adj[i][j] = 0;
      }
    }

    while (getline(inputfile, line)) {
      int v1 = stoi(line.substr(0, line.find(" "))) - 1;
      int v2 = stoi(line.substr(line.find(" ") + 1)) - 1;
      // cout << v1 << " " << v2 << endl;
      Adj[v1][v2] = 1;
      Adj[v2][v1] = 1;

    }
    inputfile.close();
  }
  return Adj;
}


int main(int argc,  char **argv) {
  std::cout << "Implementing data reduction rules..." << std::endl;
  int** Adj = makeAdjacencyMatrix();

  print(Adj, n);

  cout << "\nRule 1:" << endl;
  heavy_non_edge_rule(Adj, n, 0, 1);

  cout << "\nRule 2:" << endl;
  heavy_edge_single_end_rule(Adj, n, 0, 1);

  cout << "\nRule 3:" << endl;
  heavy_edge_single_end_rule(Adj, n, 0, 1);

  return 0;
}
