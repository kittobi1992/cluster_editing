#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <limits>

using namespace std;

int n = 0;

double INF = numeric_limits<double>::infinity();

void print(vector<vector<int> > adj) {
  // Traverse the adj[][]
  for (int i = 0; i < adj.size(); i++) {
      for (int j = 0; j < adj[i].size(); j++) {

          // Print the value at adj[i][j]
          printf("%d ", adj[i][j]);
      }
      printf("\n");
  }
}

void merge_vertices(vector<vector<int> > &adj, int u, int v) {
  // Combines the neighbors so they all onnect to u
  for (int idx = 0; idx < adj.size(); idx++) {
    if (adj[v][idx]) {
      adj[u][idx] = 1;
      adj[idx][u] = 1;
    }
  }

  // Remove v from both
  adj.erase(adj.begin()+v);
  for (int idx = 0; idx < adj.size(); idx++)
    adj[idx].erase(adj[idx].begin()+v);
}

/*
  Rule 1: heavy non edge rule
    Set an edge to FORBIDDEN (weight -INF) if the degree of the node is greadter \
    than the weight between it and the other side of the edge in question

  @params:  adj - reference to the adjacency matrix of the graph
            u - first vertex in the edge
            v - second vertex in the edge
  @return:  none
*/
void heavy_non_edge_rule(vector<vector<int> > &adj, int u, int v) {
  int sum = 0;
  for (int w = 0; w < adj.size(); w++) {
    sum += adj[u][w];
  }

  if (abs(adj[u][v]) >= sum) {
    cout << "FORBIDDEN: (" << u+1 << ", " << v+1 << ")" << endl;
    adj[u][v] = -INF;
    adj[v][u] = -INF;
  }
}


/*
  Rule 2: heavy non edge single end rule
    MERGE vertices u and v if the weight of the edge is >= to the sum of abs of the
    edges between u and every node in V\{u, v}

  @params:  adj - reference to the adjacency matrix of the graph
            u - first vertex in the edge
            v - second vertex in the edge
  @return:  none
*/
void heavy_edge_single_end_rule(vector<vector<int> > &adj, int u, int v) {
  int sum = 0;
  for (int w = 0; w < adj.size(); w++) {
    // Don't count either uu or uv
    if (w == u or w == v)
      continue;

    sum += abs(adj[u][w]);
  }

  if (adj[u][v] >= sum) {
    cout << "MERGE: (" << u+1 << ", " << v+1 << ")" << endl;
  }
}


/*
  Rule 3: heavy edge both end rule
    MERGE vertices u and v if the weight of edge (u, v) is >= to the sum of all
    edge weights between u and N(u)\{v} added to the sum of all edge weights
    between v and N(v)\{u}


  @params:  adj - reference to the adjacency matrix of the graph
            u - first vertex in the edge
            v - second vertex in the edge
  @return:  none
*/
void heavy_edge_both_end_rule(vector<vector<int> > &adj, int u, int v) {
  int sum_u = 0;
  int sum_v = 0;
  for (int w = 0; w < adj.size(); w++) {
    if (w != v)
      sum_u += adj[u][w];
    if (w != u)
      sum_v += adj[v][w];
  }

  if (adj[u][v] >= sum_u + sum_v) {
    cout << "MERGE: (" << u+1 << ", " << v+1 << ")" << endl;
  }
}



vector<vector<int> > makeAdjacencyMatrix(string fin) {
  string line;
  ifstream inputfile(fin);

  vector<vector<int> > adj;

  if (inputfile.is_open()) {
    getline(inputfile, line);

    string edges = line.substr(line.find_last_of(" ") + 1);
    string rest = line.substr(0,line.find_last_of(" "));
    string vertices = rest.substr(rest.find_last_of(" ") + 1);

    n = stoi(vertices);
    int e = stoi(edges);

    vector<int> temp(n, 0);
    for (int i = 0; i < n; i++) {
      adj.push_back(temp);
    }

    // cout << n << e << endl;

    for (int i = 0; i < n; i++) {
      // adj[i] = new int[n];
      for (int j = 0; j < n; j++) {
        adj[i][j] = 0;
      }
    }

    while (getline(inputfile, line)) {
      int v1 = stoi(line.substr(0, line.find(" "))) - 1;
      int v2 = stoi(line.substr(line.find(" ") + 1)) - 1;
      // cout << v1 << " " << v2 << endl;
      adj[v1][v2] = 1;
      adj[v2][v1] = 1;

    }
    inputfile.close();
  }
  return adj;
}


int main(int argc,  char **argv) {
  string fin = argv[1];

  std::cout << "Implementing data reduction rules: " << fin << std::endl;
  vector<vector<int> > adj = makeAdjacencyMatrix(fin);

  print(adj);

  cout << "\nRule 1:" << endl;
  // heavy_non_edge_rule(adj, 0, 1);

  print(adj);

  cout << "\nRule 2:" << endl;
  heavy_edge_single_end_rule(adj, 0, 1);

  cout << "\nRule 3:" << endl;
  heavy_edge_single_end_rule(adj, 0, 1);

  print(adj);
  cout << endl;
  merge_vertices(adj, 8, 9);
  print(adj);

  return 0;
}
