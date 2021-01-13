#include <algorithm>
#include <iostream>
#include <numeric>

#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>


using namespace std;

struct DSF {
    vector<int> parent;
    DSF(int n) : parent(n) {
        iota(begin(parent),end(parent),0);
    }
    int find(int v) {
        if(parent[v] == v) return v;
        return parent[v] = find(parent[v]);
    }
    void join(int u, int v) {
        parent[find(u)] = find(v);
    }
};


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

using Edges = vector<vector<int>>;
using IDMap = vector<vector<int>>;
struct Graph {
  Edges edges;
  IDMap idmap;

  Graph(const Edges& _edges, const IDMap& _idmap) : edges(_edges), idmap(_idmap) {};

  Graph(int n) : edges(n, vector<int>(n,-1)), idmap(n) {
    for(int i=0; i<n; ++i) idmap[i] = {i};
  }
};

// CC
vector<int> connectedComponents(Edges graph) {

  int n = size(graph);
  vector comp(n,-1);
  int numComps = 0;
  for(int i=0; i<n; ++i) {
    if(comp[i]!=-1) continue;
    // bfs
    vector q{i};
    comp[i] = numComps;
    while(!empty(q)) {
      auto v = q.back(); q.pop_back();
      for(int u=0; u<n; ++u) {
        if(comp[u]!=-1 || graph[v][u]<=0) continue;
        comp[u] = numComps;
        q.push_back(u);
      }
    }

    numComps++;
  }
  
  return comp;
}

// isClique
bool isClique(Edges graph) {
  int n = size(graph);
  for (int u = 0; u < n; ++u)
    for (int v = u+1; v < n; ++v)
      if (graph[u][v] <= 0)
        return false;
  
  return true;
}

// merge
Graph merge(Graph graph, int u, int v, int& k) {
  // u is representant
  // all below v has same index
  // all above v has index -1
  int n = size(graph.edges);
  if(u>v) swap(u,v); // u<v
  Edges merged(n-1, vector<int>(n-1, -1));
  auto newid = [n,v](int node) { return node - (node>v); };

  for(int i=0; i<n; ++i) {
    for(int j=0; j<n; ++j) {
      if(i==u || i==v) continue;
      if(j==u || j==v) continue;
      merged[newid(i)][newid(j)] = graph.edges[i][j];
    }
  }

  for(int w=0; w<n; ++w) {
      if(w==u || w==v) continue;
      // merge (u,w) and (v,w) into (u',w)
      auto uw = graph.edges[u][w];
      auto vw = graph.edges[v][w];
      merged[newid(w)][u] = uw + vw;
      merged[u][newid(w)] = uw + vw;
      if((uw<0 && vw>0) || (uw>0 && vw<0))
        k -= min(abs(uw), abs(vw));
  }

  IDMap mergedmap(n-1);
  for(int i=0; i<n; ++i) {
    if(i==v) continue;
    mergedmap[newid(i)] = graph.idmap[i];
  }
  mergedmap[u].insert(end(mergedmap[u]), begin(graph.idmap[v]), end(graph.idmap[v]));

  assert(merged[u][u]==-1);

  return {merged, mergedmap};
}

// BB recursive function
struct Solution {
  bool worked = false;
  int k = 1e9;
  vector<vector<int>> cliques;
};

auto selectBranchingEdge(const Graph& graph) {
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


vector<Graph> constructComps(const Graph& graph) {
  int n = graph.edges.size();
  auto compNum = connectedComponents(graph.edges);
  auto numComps = *max_element(begin(compNum), end(compNum)) + 1;
  vector<vector<int>> nodesPerComp(numComps);
  vector<Graph> comps;
  for(int i=0; i<n; ++i) 
    nodesPerComp[compNum[i]].push_back(i);
  for(int c=0; c<numComps; ++c) {
    Graph component(nodesPerComp[c].size());
    for (int u = 0; u < nodesPerComp[c].size(); ++u) {
      component.idmap[u] = graph.idmap[nodesPerComp[c][u]];
      for (int v = 0; v < nodesPerComp[c].size(); ++v)
        component.edges[u][v] = graph.edges[nodesPerComp[c][u]][nodesPerComp[c][v]];
    }

    comps.push_back(component);
  }

  return comps;
}

Solution solveMaybeUnconnected(Graph graph, int k, bool highL = false); // FWD


Solution solve(Graph graph, int budget, bool highL = false) {
  if (budget < 0)
    return {};

  if(isClique(graph.edges)) {
    Solution solution;
    solution.k = 0;
    solution.worked = true;
    vector<int> expandedNodeSet;
    for(int v=0; v<size(graph.edges); ++v)
      expandedNodeSet.insert(end(expandedNodeSet), begin(graph.idmap[v]), end(graph.idmap[v]));
    solution.cliques.push_back(expandedNodeSet);
    return solution;
  }

  for(int k=0; k<=budget; ++k) {

    if(highL) cout << k << endl;
    // compute lower bounds
    auto [u,v] = selectBranchingEdge(graph);
    assert(u!=v && graph.edges[u][v]>0);

    // forbidden
    auto finstance = graph;
    int delcost = finstance.edges[u][v];
    finstance.edges[u][v] = -1e9;
    finstance.edges[v][u] = -1e9;
    auto forbiddenSolution = solveMaybeUnconnected(finstance, k-delcost);
    forbiddenSolution.k += delcost;
    assert(!forbiddenSolution.worked || forbiddenSolution.k == k);
    if(forbiddenSolution.worked) 
      return forbiddenSolution;

    // merge
    int newk = k;
    auto minstance = merge(graph, u, v, newk);
    auto mergedSolution = solveMaybeUnconnected(minstance, newk);
    mergedSolution.k += (k-newk);
    assert(!mergedSolution.worked || mergedSolution.k == k);
    if(mergedSolution.worked)
      return mergedSolution;
  }

  return {};
}

Solution solveMaybeUnconnected(Graph graph, int k, bool highL) {
  Solution solution;
  solution.worked = true;
  solution.k = 0;
  for(auto comp : constructComps(graph)) {
    auto subsolution = solve(comp, k - solution.k, highL);
    if(!subsolution.worked || solution.k + subsolution.k > k) {
      solution.worked = false;
      break;
    };
    solution.k += subsolution.k;
    // add cliques of component to solution for complete graph
    for(auto clique : subsolution.cliques) {
      solution.cliques.push_back(clique);
    } 
  }
  return solution;
}


int main(int argc, char* argv[]) {
  auto edges = readGraph();

  Graph graph(edges.size());
  graph.edges = edges;

  auto solution = solveMaybeUnconnected(graph, 1e9, true);
  cout << solution.worked << endl;
  cout << "k=" << solution.k << endl;
  return 0;
}
