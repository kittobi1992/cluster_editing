
#include "reductions.h"

#include <algorithm>
#include <cassert>
#include <queue>

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
  auto newid = [v](int node) { return node - (node>v); };

  // keep all edges not not involving u or v
  for(int i=0; i<n; ++i) {
    for(int j=0; j<n; ++j) {
      if(i==u || i==v) continue;
      if(j==u || j==v) continue;
      merged[newid(i)][newid(j)] = inst.edges[i][j];
    }
  }

  // handle edges to u,v
  for(int w=0; w<n; ++w) {
      if(w==u || w==v) continue;
      // merge (u,w) and (v,w) into (u',w)
      auto uw = inst.edges[u][w];
      auto vw = inst.edges[v][w];
      assert( ! (min(uw,vw)==-INF && max(uw,vw)==INF));
      auto newval = uw+vw;
      if(min(uw,vw)==-INF) newval = -INF;
      if(max(uw,vw)==INF) newval = INF;
      merged[newid(w)][u] = newval;
      merged[u][newid(w)] = newval;
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

std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>> computeICFandICP(const Edges& edges) {
  int n = size(edges);
  vector<vector<int>> icf(n, vector<int>(n, 0));
  vector<vector<int>> icp(n, vector<int>(n, 0));
  for(int u=0; u<n; ++u) {
    for(int v=0; v<n; ++v) {
      if (u==v) continue;
      icf[u][v] += max(0, edges[u][v]);
      icp[u][v] += max(0, -edges[u][v]);

      for(int w=0; w<n; ++w) {
        if (w==u || w==v) continue;
        if(edges[u][w]>0 && edges[v][w]>0) // icf
          icf[u][v] += min(edges[u][w], edges[v][w]);

        if((edges[u][w]>0) != (edges[v][w]>0)) // icp
          icp[u][v] += min(abs(edges[u][w]), abs(edges[v][w]));
      }
    }
  }

  return {icf,icp};
}

std::optional<Instance> icxReductions(const Instance& inst, int budget) {
  auto [icf, icp] = computeICFandICP(inst.edges);
  int n = size(inst.edges);

  for(int u=0; u<n; ++u) {
    for(int v=u+1; v<n; ++v) {  

      // unsolvable instance
      if(min(icf[u][v], icp[u][v]) > budget - inst.spendCost) 
        return Instance{};

      // we must merge
      if(icf[u][v]>budget-inst.spendCost)
        return merge(inst,u,v);

      // must be forbidden
      if(icp[u][v]>budget-inst.spendCost && inst.edges[u][v]>-INF/2) { // TODO magic number?
        auto finstance = inst;
        finstance.spendCost += max(0,finstance.edges[u][v]); // cost for deletion
        finstance.edges[u][v] = -INF;
        finstance.edges[v][u] = -INF;
        return finstance;
      }
    }
  }
  return {};
}

vector<int> distances(const Instance& inst, int source) {
    int n = size(inst.edges);
    vector<int> dists(n, INF);
    queue<int> q;

    q.emplace(source);
    dists[source] = 0;

    while (!q.empty()) {
        auto node = q.front(); q.pop();
        for (int neigh=0; neigh<n; ++neigh) {
            if (inst.edges[node][neigh] > 0 && dists[neigh]==INF) {
                dists[neigh] = dists[node] + 1;
                q.emplace(neigh);
            }
        }
    }
    return dists;
}

std::optional<Instance> distance4Reduction(const Instance &inst) {
    // APSP
    int n = size(inst.edges);
    vector<vector<int>> apsp(n);
    for(int v=0; v<n; ++v)
        apsp[v] = distances(inst, v);

    bool applicable = false;
    for(auto& row : apsp)
        for(auto val : row)
            applicable |= val >= 3;
    if(!applicable)
        return {};

    auto reduced = inst;
    for(int v=0; v<n; ++v)
        for(int u=v+1; u<n; ++u)
            if(apsp[u][v]>=3)
                reduced.edges[u][v] = reduced.edges[v][u] = -INF;

    return reduced;
}
