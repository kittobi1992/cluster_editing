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

void print(vector<int> vec) {
  // Traverse the adj[][]
  for (int i = 0; i < vec.size(); i++) {
      printf("%d ", vec[i]);
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
void merge_vertices(vector<vector<int> > &adj, int u, int v,
                    vector< pair<pair<int, vector<int> >, pair<int, vector<int> > > > &merged_edges) {
  pair<int, vector<int> > u_pair;
  vector<int> u_neighbors;
  pair<int, vector<int> > v_pair;
  vector<int> v_neighbors;

  // Combines the neighbors so they all onnect to u
  for (int idx = 0; idx < adj.size(); idx++) {
    if (adj[u][idx] == 1)
      u_neighbors.push_back(idx);
    if (adj[v][idx] == 1) {
      // adj[u][idx] += adj[v][idx];
      // adj[idx][u] += adj[v][idx];
      adj[u][idx] = 1;
      adj[idx][u] = 1;
      v_neighbors.push_back(idx);
    }


  }

  // Remove v from both
  adj.erase(adj.begin()+v);
  for (int idx = 0; idx < adj.size(); idx++)
    if (v < adj[idx].size())
      adj[idx].erase(adj[idx].begin()+v);

  u_pair.first  = u;
  u_pair.second  = u_neighbors;
  v_pair.first  = v;
  v_pair.second  = v_neighbors;

  pair<pair<int, vector<int> >, pair<int, vector<int> > > merged_vertices;
  merged_vertices.first = u_pair;
  merged_vertices.second = v_pair;
  merged_edges.push_back(merged_vertices);

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
bool heavy_edge_single_end_rule(vector<vector<int> > &adj, int u, int v,
                                vector< pair<pair<int, vector<int> >, pair<int, vector<int> > > > &merged_edges) {
  if (adj.size() <= u || adj.size() <= v) return false;
  if (u == -INF || v == -INF) return false;
  if (u == v) return false;

  int sum = 0;
  for (int w = 0; w < adj.size(); w++) {
    // Don't count either uu or uv
    if (w == u or w == v)
      continue;

    sum += abs(adj[u][w]);
  }

  // cout << "aa" << endl;

  if (adj[u][v] >= sum) {
    // cout << "bb: (" << u << ", " << v << ")" << endl;
    // cout << "MERGE: (" << u << ", " << v << ")" << endl;
    cout << "rule 2 " << u << v << endl;
    merge_vertices(adj, u, v, merged_edges);
    // cout << "cc" << endl;
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
bool heavy_edge_both_end_rule(vector<vector<int> > &adj, int u, int v,
                              vector< pair<pair<int, vector<int> >, pair<int, vector<int> > > > &merged_edges) {
  if (adj.size() <= u || adj.size() <= v) return false;
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
    cout << "rule 3" << endl;
    merge_vertices(adj, u, v, merged_edges);
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

vector<int> apply1(vector<vector<int> > &adj) {
  vector<int> applied_edges;
  for (int u = 0; u < adj.size(); u++) {
    for (int v = 0; v < adj[u].size(); v++) {
      if (heavy_non_edge_rule(adj, u, v)) {
        applied_edges.push_back(u);
        applied_edges.push_back(v);
      }
    }
  }
  return applied_edges;
}

void apply2(vector<vector<int> > &adj,
            vector< pair<pair<int, vector<int> >, pair<int, vector<int> > > > &merged_edges) {
  for (int u = 0; u < adj.size(); u++) {
    for (int v = 0; v < adj[u].size(); v++) {
      while(heavy_edge_single_end_rule(adj, u, v, merged_edges))
        continue;
    }
  }
}

void apply3(vector<vector<int> > &adj,
            vector< pair<pair<int, vector<int> >, pair<int, vector<int> > > > &merged_edges) {
  for (int u = 0; u < adj.size(); u++) {
    for (int v = 0; v < adj[u].size(); v++) {
      while(heavy_edge_both_end_rule(adj, u, v, merged_edges))
        continue;
        // print(adj);
    }
  }
}


/*
  expand_to_original_size
  expand the graph to add rows and cols to the original size

  @params
    - reduced: the reduced adj matrix
    - original_size: the original size of the adj matrix
*/
vector<vector<int> > expand_to_original_size(vector<vector<int > > reduced,
                                             int original_size) {
  // Get the original size to expand into
  vector<vector<int> > rebuilt_graph;
  for (int i = 0; i < original_size; i++) {

   // Make the empty row
   vector<int> temp(original_size, 0);

   // fill the row with what is in the reduced graph
   if (i < reduced.size()) {
     for (int j = 0; j < reduced.size(); j++) {
       temp[j] = reduced[i][j];
     }
     rebuilt_graph.push_back(temp);
   }
   else {
     rebuilt_graph.push_back(temp);
   }
 }

  // Unmerge edges
  return rebuilt_graph;

}

/*
  expand_graph
  shift everything in row/col of idx to the right and down to make room for a spot there

  @params
   - adj: reference to the adj matrix we are expanding
   - idx: the idx of the row/col we have to add in
*/
void expand_graph(vector<vector<int> > &adj, int idx) {
  for (int i = adj.size()-1; i > idx; i--) {
    for (int j = adj.size()-1; j > idx; j--) {
      adj[i][j] = adj[i-1][j-1];
    }
  }
  for (int i = 0; i < adj.size(); i++) {
    adj[i][idx] = 0;
    adj[idx][i] = 0;
  }
}

/*
  rebuild_graph
  Rebuild the graph from the reductions we have applied to the original state

  @params
   - reduced: the reduced graph
   - merged_edges: the vector of all the reductions we have applied from previous functions
   - original_size: the original size of the adj matrix
  @return
   - the rebuilt graph in the form of an adj matrix
*/
vector<vector<int> > rebuild_graph(vector<vector<int > > reduced,
                                   vector< pair<pair<int, vector<int> >, pair<int, vector<int> > > > merged_edges,
                                   int original_size) {

  cout << endl << "REBUILDING..." << endl;

  // Get the original size to expand into
  vector<vector<int> > rebuilt_graph = expand_to_original_size(reduced, original_size);
  pair<pair<int, vector<int> >, pair<int, vector<int> > > reduction;
  pair<int, vector<int> > u_edge;
  pair<int, vector<int> > v_edge;
  int u;
  int v;
  vector<int> u_neighbors;
  vector<int> v_neighbors;

  for (int i = merged_edges.size()-1; i > -1; i--) {
    reduction = merged_edges[i];

    // Grab the vertext and its neighbors from the reduction pair
    u_edge = reduction.first;
    v_edge = reduction.second;

    // Pull out the vertex value and its neighbors
    u =  u_edge.first;
    u_neighbors = u_edge.second;
    v = v_edge.first;
    v_neighbors = v_edge.second;

    // Expand the graph so v has its row and column back
    expand_graph(rebuilt_graph, v);

    // clear u's neighbors
    for (int idx = 0; idx < rebuilt_graph.size(); idx++) {
      rebuilt_graph[u][idx] = 0;
      rebuilt_graph[idx][u] = 0;
    }

    // Set u's neighbors
    for (int idx = 0; idx < u_neighbors.size(); idx++) {
      rebuilt_graph[u][u_neighbors[idx]] = 1;
      rebuilt_graph[u_neighbors[idx]][u] = 1;
    }

    // Set v's neighbors
    for (int idx = 0; idx < v_neighbors.size(); idx++) {
      rebuilt_graph[v][v_neighbors[idx]] = 1;
      rebuilt_graph[v_neighbors[idx]][v] = 1;
    }
  }

  return rebuilt_graph;
}

int main(int argc,  char **argv) {
  string fin = argv[1];

  // std::cout << "Implementing data reduction rules: " << fin << std::endl;
  vector<vector<int> > adj = makeAdjacencyMatrix(fin);
  int total = adj.size();
  int sum = 0;
  int total_cost = 0;
  for (int u = 0; u < adj.size(); u++) {
    for (int v = 0; v < adj.size(); v++) {
      sum += adj[u][v];
    }
  }
  print(adj);
  int og_size = adj.size();

  vector<vector<int> > clique;
  vector<int> applied_edges;

  // We will store the merges as a vector of pairs, so it will be like
  // a pair of these: pair<vertex, neighbors> all shoved in a vector.
  vector< pair<pair<int, vector<int> >, pair<int, vector<int> > > > merged_edges;
  for (int u = 0; u < adj.size(); u++) {
    for (int v = 0; v < adj.size(); v++) {
      if (u != v) {
        vector<int> emptyRow (adj[u].size(), 0);
        vector<vector<int> > clique (adj.size(), emptyRow);
        if (computeCriticalClique(adj, clique, u, v)) {
          // print(clique);
          // cout << endl;
          apply1(clique);
          // print(applied_edges);
          apply2(clique, merged_edges);
          apply3(clique, merged_edges);
          // print(clique);
          total -= (total - clique.size());

          for (int i = 0; i < clique.size(); i++) {
            for (int j = 0; j < clique.size(); j++) {
              if (clique[i][j] != -INF) {
                sum = clique[i][j];
              }
            }
          }
        }
      }
    }
  }

  // for (int i = 0; i < merged_edges.size(); i++) {
  //   cout << i << ": " << merged_edges[i].first.first << " merged with " << merged_edges[i].second.first << endl;
  //   cout << "\tu neighbors: ";
  //   for (int u = 0; u < merged_edges[i].first.second.size(); u++)
  //     cout << merged_edges[i].first.second[u] << ", ";
  //   cout << endl << "\tv neighbors: ";
  //   for (int v = 0; v < merged_edges[i].second.second.size(); v++) {
  //     cout << merged_edges[i].second.second[v] << ", ";
  //   }
  //   cout << endl;
  // }


  int s = total_cost-sum;
  if (total_cost-sum <= 0) s = -1;

  cout << fin << ", " << n << ", " << total << ",  " <<  s << endl;

  print(rebuild_graph(adj, merged_edges, og_size));

  return 0;
}
