
#include <algorithm>
#include <iostream>
#include <numeric>
#include <sstream>
#include <cassert>

#include <cluster_editing/exact/instance.h>
#include <cluster_editing/exact/reductions.h>


using namespace std;

// read graph into matrix
vector<vector<int>> readGraph() {

  std::istringstream sstream;
  auto getline = [&]() {
    std::string line;
    do {
      std::getline(cin, line);
    } while ( line[0] == 'c' );
    sstream = std::istringstream(line);
  };


  getline();
  std::string skip;
  int n, m;
  sstream >> skip >> skip >> n >> m;

  vector<vector<int>> res(n, vector<int>(n, -1));

  for ( int i = 0; i < m; ++i ) {
    getline();

    int u, v;
    sstream >> u >> v;
    --u; --v;
    res[u][v] = 1;
    res[v][u] = 1;
  }

  return res;
}

// isClique
bool isClique(const Edges& graph) {
  int n = size(graph);
  for (int u = 0; u < n; ++u)
    for (int v = u+1; v < n; ++v)
      if (graph[u][v] <= 0)
        return false;
  
  return true;
}


// BB recursive function
struct Solution {
  bool worked = false;
  int cost = 1e9;
  vector<vector<int>> cliques;
};

auto selectBranchingEdge(const Instance& graph) {
  // find conflict triple
  int n = size(graph.edges);
  for(int v=0; v<n; ++v) {
    for(int u=0; u<n; ++u) {
      if (graph.edges[u][v] <= 0) continue; // uv must be edge
      for(int w=0; w<n; ++w) {
        if(w==u || w==v) continue;
        if(graph.edges[u][w]<=0) continue; // uw must be edge
        if(graph.edges[v][w]>0) continue; // vw must be non-edge
        return pair(u,v);
      }
    }
  }
  return pair(-1,-1);
}


Solution solveMaybeUnconnected(Instance graph, int budget, bool highL = false); // FWD


Solution solve(Instance graph, int budget, bool highL = false) {
  if (budget < graph.spendCost)
    return {};

  if(isClique(graph.edges)) {
    Solution solution;
    solution.cost = graph.spendCost;
    solution.worked = true;
    vector<int> expandedNodeSet;
    for(int v=0; v<size(graph.edges); ++v)
      expandedNodeSet.insert(end(expandedNodeSet), begin(graph.idmap[v]), end(graph.idmap[v]));
    solution.cliques.push_back(expandedNodeSet);
    return solution;
  }

  // compute some lower bounds
  int n = size(graph.edges);
  vector<vector<int>> icf(n, vector<int>(n, 0));
  vector<vector<int>> icp(n, vector<int>(n, 0));
  for(int u=0; u<n; ++u) {
    for(int v=0; v<n; ++v) {
      if (u==v) continue;
      icf[u][v] += max(0, graph.edges[u][v]);
      icp[u][v] += max(0, -graph.edges[u][v]);

      for(int w=0; w<n; ++w) {
        if (w==u || w==v) continue;
        if(graph.edges[u][w]>0 && graph.edges[v][w]>0) // icf
          icf[u][v] += min(graph.edges[u][w], graph.edges[v][w]);

        if((graph.edges[u][w]>0) != (graph.edges[v][w]>0)) // icp
          icp[u][v] += min(abs(graph.edges[u][w]), abs(graph.edges[v][w]));
      }

      if(min(icf[u][v], icp[u][v]) > budget - graph.spendCost) 
        return {};
    }
  }

  for(int u=0; u<n; ++u) {
    for(int v=u+1; v<n; ++v) {

      // we must merge
      if(icf[u][v]>budget)
        return solveMaybeUnconnected(merge(graph,u,v), budget);

      // must be forbidden
      if(icp[u][v]+graph.spendCost>budget && graph.edges[u][v]>-1e7) { // TODO magic number?
        auto finstance = graph;
        finstance.spendCost += max(0,finstance.edges[u][v]); // cost for deletion
        finstance.edges[u][v] = -1e8;
        finstance.edges[v][u] = -1e8;
        return solveMaybeUnconnected(finstance, budget);
      }
    }
  }

  // here we choose the edge to branch on
  auto [u,v] = selectBranchingEdge(graph);
  int bestReduction = min(icp[u][v], icf[u][v]);
  for(int i=0; i<n; ++i) {
    for(int j=i+1; j<n; ++j) {
      if(graph.edges[i][j]>0 && min(icf[i][j], icp[i][j]) > bestReduction) {
        bestReduction =  min(icf[i][j], icp[i][j]);
        u = i;
        v = j;
      }
    }
  }
  assert(u!=v && graph.edges[u][v]>0);


  auto minstance = merge(graph,u,v); // merged instance
  auto finstance = graph; // forbidden instance
  finstance.spendCost += max(0,graph.edges[u][v]); // cost for deletion
  finstance.edges[u][v] = -1e8;
  finstance.edges[v][u] = -1e8;

  for(int k=graph.spendCost; k<=budget; ++k) {

    if(highL) 
      cout << k << endl; // debug stuff

    // try merging
    auto mergedSolution = solveMaybeUnconnected(minstance, k);
    assert(!mergedSolution.worked || mergedSolution.cost == k);
    if(mergedSolution.worked)
      return mergedSolution;

    // try permanent deletion
    auto forbiddenSolution = solveMaybeUnconnected(finstance, k);
    assert(!forbiddenSolution.worked || forbiddenSolution.cost == k);
    if(forbiddenSolution.worked) 
      return forbiddenSolution;
  }

  return {};
}

Solution solveMaybeUnconnected(Instance graph, int budget, bool highL) {
  if(budget<graph.spendCost) return {}; // that should probably never happen
  Solution solution;
  solution.worked = true;
  solution.cost = graph.spendCost;
  for(auto comp : constructConnectedComponents(graph)) {
    auto subsolution = solve(comp, budget - solution.cost, highL);
    if(!subsolution.worked || solution.cost + subsolution.cost > budget) {
      solution.worked = false;
      break;
    };
    solution.cost += subsolution.cost;
    // add cliques of component to solution for complete graph
    for(auto clique : subsolution.cliques) 
      solution.cliques.push_back(clique);
  }
  return solution;
}


int main(int argc, char* argv[]) {

  /*
  we can solve
  - 001 (k=3)
  - 065 (k=128) longest with 16s
  - 077 (k=78)
  - 137 (k=16)
  - 153 (k=6)
  - 155 (k=63)
  - 173 (k=100)
  */
  string file = "../../../instances/exact/exact137.gr";
  if(argc>1)
    file = argv[1];
  freopen(file.c_str(), "r", stdin);
  cout << "reading file " << file << endl;
  auto edges = readGraph();
  cout << "n = " << size(edges) << endl;

  Instance graph(edges.size());
  graph.edges = edges;

  auto solution = solveMaybeUnconnected(graph, 1e9, true);
  cout << solution.worked << endl;
  cout << "k=" << solution.cost << endl;
  return 0;
}
