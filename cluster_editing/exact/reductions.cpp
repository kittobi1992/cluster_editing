
#include "reductions.h"

#include <algorithm>
#include <cassert>

using namespace std;


vector<int> connectedComponents(const Edges& graph) {

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

vector<Instance> constructConnectedComponents(const Instance& graph) {
  int n = graph.edges.size();
  auto compNum = connectedComponents(graph.edges);
  auto numComps = *max_element(begin(compNum), end(compNum)) + 1;
  vector<vector<int>> nodesPerComp(numComps);
  vector<Instance> comps;
  for(int i=0; i<n; ++i) 
    nodesPerComp[compNum[i]].push_back(i);
  for(int c=0; c<numComps; ++c) {
    Instance component(nodesPerComp[c].size());
    for (int u = 0; u < nodesPerComp[c].size(); ++u) {
      component.idmap[u] = graph.idmap[nodesPerComp[c][u]];
      for (int v = 0; v < nodesPerComp[c].size(); ++v)
        component.edges[u][v] = graph.edges[nodesPerComp[c][u]][nodesPerComp[c][v]];
    }

    comps.push_back(component);
  }

  return comps;
}


Instance merge(const Instance& inst, int u, int v) {
  int n = size(inst.edges);

  int mergeCost = 0;
  if(inst.edges[u][v]<0) // we must insert the edge if it did not exist yet
    mergeCost += abs(inst.edges[u][v]);

  // u is representant
  // all below v has same index
  // all above v has index -1
  if(u>v) swap(u,v); // u<v
  Edges merged(n-1, vector<int>(n-1, -1));
  auto newid = [n,v](int node) { return node - (node>v); };

  // keep all edges not not involving u or v
  for(int i=0; i<n; ++i) {
    for(int j=0; j<n; ++j) {
      if(i==u || i==v) continue;
      if(j==u || j==v) continue;
      merged[newid(i)][newid(j)] = inst.edges[i][j];
    }
  }

  // hanlde edges to u,v
  for(int w=0; w<n; ++w) {
      if(w==u || w==v) continue;
      // merge (u,w) and (v,w) into (u',w)
      auto uw = inst.edges[u][w];
      auto vw = inst.edges[v][w];
      merged[newid(w)][u] = uw + vw;
      merged[u][newid(w)] = uw + vw;
      if((uw<0 && vw>0) || (uw>0 && vw<0))
        mergeCost += min(abs(uw), abs(vw));
  }
  assert(merged[u][u]==-1);

  IDMap mergedmap(n-1);
  // copy all old ones to the right place
  for(int i=0; i<n; ++i)
    if(i!=v) 
      mergedmap[newid(i)] = inst.idmap[i];
  // insert all from v into u
  mergedmap[u].insert(end(mergedmap[u]), begin(inst.idmap[v]), end(inst.idmap[v]));

  Instance result{merged, mergedmap};
  result.spendCost = inst.spendCost + mergeCost;
  return result;
}