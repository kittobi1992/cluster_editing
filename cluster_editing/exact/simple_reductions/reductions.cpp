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
          if (adj[i][j] <= -1000000)
            cout << ". ";
          else
            printf("%d ", adj[i][j]);
      }
      printf("\n");
  }
}


bool computeCriticalClique(vector<vector<int> > adj, vector<vector<int> > &newAdj, int u, int v) {
  // vector<int> u_neighbors = adj[u];
  // vector<int> v_neighbors = adj[v];
  // vector<int> emptyRow (adj[u].size(), 0);
  // vector<vector<int> > newAdj (adj.size(), emptyRow);

  if (u > v) return false;

  // cout << u << " " << v << endl << endl;

  // for (int i = 0; i < adj.size(); i++) {
  //   for (int j = 0; j < adj[u].size(); j++) {
  //     if (i == j)
  //       adj[i][j] = 1;
  //   }
  // }

  bool hasNeighbors = false;
  for (int idx = 0; idx < adj[u].size(); idx++) {
    if (adj[u][idx] == 1 || adj[v][idx] == 1)
      hasNeighbors = true;

    if (adj[u][idx] != adj[idx][v] && u != idx && v != idx) {
      // cout << "(" << u << ", " << v << ")" << endl;
      // cout << idx << ": " << adj[u][idx] << " " << adj[idx][v] << endl << endl;
      return false;
    }
  }

  if (!hasNeighbors) return false;

  for (int i = 0; i < adj.size(); i++) {
    for (int j = 0; j < adj[u].size(); j++) {
      if (i == u || i == v || j == u || j == v) {
        // cout << i << " " << j  << endl;
        newAdj[i][j] = adj[i][j];
      }
    }
  }
  // cout << u << " " << v << endl;
  // print(newAdj);
  // cout << endl;

  // cout << "(" << u <<  ", " << v << ")"<< endl;
  return true;
}


/*
  Merge
    MERGE two vertices together

  TODO:
    How to merge when they have weights? do we need to consider this?
*/
void merge_vertices(vector<vector<int> > &adj, int u, int v) {
  // Combines the neighbors so they all onnect to u
  for (int idx = 0; idx < adj.size(); idx++) {
    if (adj[v][idx] != 0) {
      adj[u][idx] += adj[v][idx];
      adj[idx][u] += adj[v][idx];

      // // The else statement is the normal merge, this is trying to deal with seg faults
      // // from we think adding negative infinity
      // if (!(adj[v][idx] == -INF || adj[u][idx] == -INF || adj[idx][u] == -INF)) {
      //   adj[u][idx] = -INF;
      //   adj[idx][u] = -INF;
      // } else {
      //   adj[u][idx] += adj[v][idx];
      //   adj[idx][u] += adj[v][idx];
      // }
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
bool heavy_non_edge_rule(vector<vector<int> > &adj, int u, int v) {
  if (adj.size() < u || adj.size() < v) return false;
  if (u == -INF || v == -INF) return false;
  if (u == v) return false;

  int sum = 0;
  for (int w = 0; w < adj.size(); w++) {
    sum += adj[u][w];
  }

  if (abs(adj[u][v]) >= sum) {
    // cout << "FORBIDDEN: (" << u+1 << ", " << v+1 << ")" << endl;
    adj[u][v] = -INF;
    adj[v][u] = -INF;
    return true;
  }

  return false;
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
bool heavy_edge_single_end_rule(vector<vector<int> > &adj, int u, int v) {
  if (adj.size() < u || adj.size() < v) return false;
  if (u == -INF || v == -INF) return false;
  if (u == v) return false;

  int sum = 0;
  for (int w = 0; w < adj.size(); w++) {
    // Don't count either uu or uv
    if (w == u or w == v)
      continue;

    sum += abs(adj[u][w]);
  }

  if (adj[u][v] >= sum) {
    // cout << "MERGE: (" << u << ", " << v << ")" << endl;
    // merge_vertices(adj, u, v);
    return true;
  }

  return false;
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
bool heavy_edge_both_end_rule(vector<vector<int> > &adj, int u, int v) {
  if (adj.size() < u || adj.size() < v) return false;
  if (u == -INF || v == -INF) return false;
  if (u == v) return false;

  int sum_u = 0;
  int sum_v = 0;
  for (int w = 0; w < adj.size(); w++) {
    if (w != v)
      sum_u += adj[u][w];
    if (w != u)
      sum_v += adj[v][w];
  }

  if (adj[u][v] >= sum_u + sum_v) {
    // cout << "MERGE: (" << u << ", " << v << ")" << endl;
    // merge_vertices(adj, u, v);
    return true;
  }

  return false;
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

void apply1(vector<vector<int> > &adj) {
  for (int u = 0; u < adj.size(); u++) {
    for (int v = 0; v < adj[u].size(); v++) {
      heavy_non_edge_rule(adj, u, v);
    }
  }
}

void apply2(vector<vector<int> > &adj) {
  for (int u = 0; u < adj.size(); u++) {
    for (int v = 0; v < adj[u].size(); v++) {
      while(heavy_edge_single_end_rule(adj, u, v))
        continue;
        // print(adj);
    }
  }
}

void apply3(vector<vector<int> > &adj) {
  for (int u = 0; u < adj.size(); u++) {
    for (int v = 0; v < adj[u].size(); v++) {
      while(heavy_edge_both_end_rule(adj, u, v))
        continue;
        // print(adj);
    }
  }
}

int main(int argc,  char **argv) {
  string fin = argv[1];


  // std::cout << "Implementing data reduction rules: " << fin << std::endl;
  vector<vector<int> > adj = makeAdjacencyMatrix(fin);
  int total = adj.size();

  // print(adj);

  vector<vector<int> > clique;
  for (int u = 0; u < adj.size(); u++) {
    for (int v = 0; v < adj.size(); v++) {
      if (u != v) {
        vector<int> emptyRow (adj[u].size(), 0);
        vector<vector<int> > clique (adj.size(), emptyRow);
        if (computeCriticalClique(adj, clique, u, v)) {
          // print(clique);
          // cout << endl;
          // apply1(clique);
          apply2(clique);
          apply3(clique);
          // print(clique);
          total -= (total - clique.size());
          // if (clique.size() != n)
          //   cout << fin << ", " << n << ", " << clique.size() << endl;
        }
      }
    }
  }

  // print(adj);

  // // Apply Rule 1 to every pair of vertices
  // for (int u = 0; u < adj.size(); u++) {
  //   for (int v = 0; v < adj[u].size(); v++) {
  //     heavy_non_edge_rule(adj, u, v);
  //   }
  // }

  // print(adj);
  // cout << endl << "Rule 2" << endl;

  // Apply Rule 2 to every pair of vertices
  // for (int u = 0; u < adj.size(); u++) {
  //   for (int v = 0; v < adj[u].size(); v++) {
  //     while(heavy_edge_single_end_rule(adj, u, v))
  //       continue;
  //       // print(adj);
  //   }
  // }

  // cout << endl << "Rule 3" << endl;

  // Apply Rule 3 to every pair of vertices
  // for (int u = 0; u < adj.size(); u++) {
  //   for (int v = 0; v < adj[u].size(); v++) {
  //     while(heavy_edge_both_end_rule(adj, u, v))
  //       continue;
  //       // print(adj);
  //   }
  // }

  cout << fin << ", " << n << ", " << total << endl;

  return 0;
}
