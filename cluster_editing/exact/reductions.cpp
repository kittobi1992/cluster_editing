/*******************************************************************************
 * This file is part of KaPoCE.
 *
 * KaPoCE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaPoCE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaPoCE.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "reductions.h"

#include <algorithm>
#include <cassert>
#include <queue>
#include <iostream>

#include <cluster_editing/exact/lower_bounds.h>
#include <cluster_editing/exact/solver.h>

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
      for (int v = 0; v < nodesPerComp[c].size(); ++v) {
          component.edges[u][v] = graph.edges[nodesPerComp[c][u]][nodesPerComp[c][v]];
          component.orig[u][v] = graph.orig[nodesPerComp[c][u]][nodesPerComp[c][v]];
      }
    }

    comps.push_back(component);
  }

  return comps;
}


Instance merge(const Instance& inst, int u, int v) {
  int n = size(inst.edges);

  for(int i=0; i<n; ++i)
      for(int j=0; j<n; ++j)
          assert(abs(inst.edges[i][j])==INF || inst.edges[i][j]==inst.orig[i][j]);

  int mergeCost = 0;
  if(inst.edges[u][v]<0) // we must insert the edge if it did not exist yet
    mergeCost += abs(inst.edges[u][v]);

  // u is representant
  // all below v has same index
  // all above v has index -1
  if(u>v) swap(u,v); // u<v
  Edges merged(n-1, vector<int>(n-1, -1));
    Edges merged_orig(n-1, vector<int>(n-1, -1));
  auto newid = [v](int node) { return node - (node>v); };

  // keep all edges not not involving u or v
  for(int i=0; i<n; ++i) {
    for(int j=0; j<n; ++j) {
      if(i==u || i==v) continue;
      if(j==u || j==v) continue;
      merged[newid(i)][newid(j)] = inst.edges[i][j];
        merged_orig[newid(i)][newid(j)] = inst.orig[i][j];
    }
  }

  // handle edges to u,v
  for(int w=0; w<n; ++w) {
      if(w==u || w==v) continue;
      // merge (u,w) and (v,w) into (u',w)
      auto uw = inst.edges[u][w];
      auto vw = inst.edges[v][w];
      merged_orig[newid(w)][u] = merged_orig[u][newid(w)] = inst.orig[u][w] + inst.orig[v][w];
      if(min(uw,vw)==-INF && max(uw,vw)==INF)
          return {}; // in this case the instance is not solvable
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
  result.orig = merged_orig;
  result.done_clusters = inst.done_clusters;

  // debug checks
  for(int i=0; i<n-1; ++i)
      for(int j=0; j<n-1; ++j)
          assert(abs(result.edges[i][j])==INF || result.edges[i][j] == result.orig[i][j]);

  return result;
}

void mergeAllINF(Instance& inst) {

    for(int v=0; v<size(inst.edges); ++v) {
        auto it = find(begin(inst.edges[v]), end(inst.edges[v]), INF);
        while(it!=end(inst.edges[v])) {
            assert(v!=it-begin(inst.edges[v]));
            inst = merge(inst, v, it-begin(inst.edges[v]));
            if(inst.spendCost==INF) return; // merge had cost infinity
            it = find(begin(inst.edges[v]), end(inst.edges[v]), INF);
        }
    }
    for(auto& row : inst.edges)
        for(auto v : row)
            assert(v!=INF);
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

std::optional<Instance> forcedChoices(const Instance& inst, int upper_bound, bool verbose) {
    auto packing = getAPacking(inst);
    auto potential = inst.edges;
    auto lower = 0;
    for(auto t : packing) lower += t.apply(potential);

    int n=size(inst.edges);
    // contains for each edge u,v the triple in the packing corresponding to that edges (if any)
    // -1 if no triple
    // -2 if more than one triple
    vector tripleMap(n, vector(n, -1));
    for(int i=0; i<size(packing); ++i) {
        auto [u,v,w,cost,valid] = packing[i];
        for(auto [a,b] : {pair(u,v), pair(u,w), pair(v,w)})
            tripleMap[a][b] = tripleMap[b][a] = (tripleMap[a][b]==-1 ? i : -2);
    }

    vector<tuple<int,int,int>> forced;
    for(int u=0; u<n; ++u) {
        for(int v=u+1; v<n; ++v) {
            if(inst.edges[u][v]==-INF) continue;

            // remove triple
            auto old = potential[u][v];
            //if(tripleMap[u][v]==-2) continue; // TODO for now only try edges that have 0 or 1 assinged triples
            if(tripleMap[u][v] >= 0) {
                auto t = packing[tripleMap[u][v]];
                lower += t.apply(potential,true);
                old = potential[u][v];
            }

            // see what happens if i set uv to permanent or forbidden
            for(auto choice : {INF, -INF}) {
                int increase = tripleMap[u][v]==-2 ? -abs(inst.edges[u][v]-old) : 0;
                auto uv_cost = inst.edges[u][v];
                if(choice==INF && uv_cost<0) increase += abs(uv_cost);
                if(choice==-INF && uv_cost>0) increase += uv_cost;
                potential[u][v] = potential[v][u] = choice;
                for(int x=0; x<n; ++x) {
                    auto t = Triple(u,v,x,potential);
                    if(t.valid) increase += t.cost;
                }
                if(lower+increase+inst.spendCost > upper_bound && -choice != inst.edges[u][v])
                    forced.emplace_back(u,v,-choice);
            }

            // reapply triple
            potential[u][v] = potential[v][u] = old;
            if(tripleMap[u][v] >= 0) {
                auto t = packing[tripleMap[u][v]];
                lower += t.apply(potential);
            }
        }
    }

    int perms = 0;
    for(auto [u,v,c] : forced) perms += (c==INF);
    if(verbose) cout << "Forced: " << size(forced) << "\tPerm: " << perms << "\tForb: " << size(forced)-perms << endl;
    if(empty(forced)) return {};

    // apply all forced choices
    auto res = inst;
    for(auto [u,v,c] : forced) {
        if(c==INF && res.edges[u][v]<0) res.spendCost -= res.edges[u][v];
        if(c==-INF&& res.edges[u][v]>0) res.spendCost += res.edges[u][v];
        res.edges[u][v] = res.edges[v][u] = c;
        if(res.spendCost>upper_bound)
            return Instance{}; // not solvable
    }

    mergeAllINF(res);

    return res;
}

std::optional<Instance> simpleTwin(const Instance &inst) {

    int n = size(inst.edges);
    if(n<3) return {};

    Instance res = inst;
    bool found = false;
    for(int v=0; v<n; ++v) {
        for(int u=v+1; u<n; ++u) {
            if(inst.edges[u][v]<0) continue;

            bool haveFactor = false;
            int a=1, b=1;
            int w=0;
            for(; w<n; ++w) {
                if(w==u || w==v) continue;
                auto uw = inst.edges[u][w];
                auto vw = inst.edges[v][w];
                // handle all 0 and -INF cases, break for conflicts
                if((uw==0) != (vw==0)) break;
                if((uw==-INF) != (vw==-INF)) break;
                if((uw==0) && (vw==0)) continue;
                if((uw==-INF) && (vw==-INF)) continue;
                // handle conflicts from differing sign
                if((uw>0) != (vw>0)) break;
                // (possibly set and) check if factor is the same for all
                if(!haveFactor)
                  a = vw, b=uw, haveFactor=true;
                if(a*uw != b*vw) break;
            }
            // if no conflict until now we have a twin
            if(w==n) { // if loop went until the end
                res.edges[u][v] = res.edges[v][u] = INF;
                found = true;
            }
        }
    }

    if(!found) return {};

    mergeAllINF(res);

    return res;
}


int max_over_subsetsDP(vector<tuple<int,int,int>>& B, long long delta_u, long long delta_v) {

    int range_x = 0;
    for(auto [x,y,f] : B) {
        assert(min(x,y)!=-INF);
        range_x += abs(x);
    }

    vector last(2*range_x + 1, -INF);
    last[range_x] = 0; // at index range_x is there point point 0 of the interval [-range, +range]
    for(auto [x,y,f] : B) {

        auto next = last; // OPTION 0: don't use pair
        for(int i=0; i<size(last); ++i) {
            if(last[i]==-INF) continue;
            // valid values should never be out of x_range
            assert(0<=i-x && i+x<size(last));
            if(f!=1) next[i+x] = max(next[i+x], last[i] - y); // OPTION 1: put in bucket 1
            if(f!=2) next[i-x] = max(next[i-x], last[i] + y); // OPTION 2: put in bucket 2
        }

        swap(next,last);
    }

    auto res = min(delta_u, delta_v);
    for(int i=0; i<size(last); ++i) {
        auto x = i-range_x;
        auto y = last[i];
        res = max(res, min(x+delta_u,y+delta_v));
    }

    return res;
}

std::optional<Instance> complexTwin(const Instance &inst, bool calc_dp) {
    if(size(inst.edges)<2) return {};
    int n = size(inst.edges);
    auto& g = inst.edges;
    auto& o = inst.orig;

    vector<pair<int,int>> perms;

    bool usedEq = false;
    for(int u=0; u<n; ++u) {
        for(int v=u+1; v<n; ++v) {
            if(g[u][v]<=0) continue;

            vector W(n,false);
            long long delta_u=0, delta_v=0; // 64 bit to sum up all those -INFs
            for(int w=0; w<n; ++w) {
                if(w==u || w==v) continue;
                if(g[u][w]==-INF && g[v][w]==-INF) continue;
                if(min(g[u][w],g[v][w])>-INF && (g[u][w]>0) != (g[v][w]>0)) { // no -INFs in play and w is in some exclusive neighborhood
                    delta_u += abs(g[u][w]);
                    delta_v += abs(g[v][w]);
                    continue;
                }
                W[w] = true; // either -INFs in play or (+,+) or (-,-)
            }

            assert(delta_u<INF && delta_v<INF);
            assert(delta_u>=0 && delta_v>=0);

            // try early breaks
            if(g[u][v]<min(delta_u, delta_v)) continue;
            auto upper_bound = 2*delta_u + 2*delta_v;
            for(int w=0; w<n; ++w)
                if(W[w]) upper_bound += abs(o[u][w] - o[v][w]);
            if(2*g[u][v]>upper_bound) {
                perms.emplace_back(u,v);
                continue;
            }

            if(!calc_dp) continue;
            vector<tuple<int,int,int>> B; // uw, vw, forbiddenBag
            for(int w=0; w<n; ++w) if(W[w]) {
                int forbidden = -1;
                if(g[u][w]==-INF) forbidden = 1;
                if(g[v][w]==-INF) forbidden = 2;
                assert(max(g[u][w],g[v][w])>-INF);
                B.emplace_back(o[u][w], o[v][w], forbidden); // use original here!
            }
            auto val = max_over_subsetsDP(B, delta_u, delta_v);
            assert(2*val <= upper_bound);
            if(g[u][v]>val)
                perms.emplace_back(u,v);
            else if(g[u][v]==val && !usedEq) {
                perms.emplace_back(u,v);
                usedEq = true;
            }
        }
    }

    if(empty(perms)) return {};
    auto res = inst;
    for(auto [u,v] : perms) res.edges[u][v] = res.edges[v][u] = INF;
    mergeAllINF(res);
    return res;
}

std::optional<Instance> removeCliques(const Instance &inst) {
    int n = size(inst.edges);
    auto compId = connectedComponents(inst.edges);
    int num = 1+ *max_element(begin(compId), end(compId));
    vector isClique(num, true);
    vector comps(num, vector<int>());
    for(int v=0; v<n; ++v) comps[compId[v]].push_back(v);
    for(int v=0; v<n; ++v)
        for(int u=v+1; u<n; ++u)
            if(compId[v]==compId[u] && inst.edges[v][u]<0)
                isClique[compId[v]] = false;

    if(count(begin(isClique), end(isClique), true) == 0)
        return {};

    Instance res = inst;
    vector<int> toDelete;
    for(int i=0; i<num; ++i) {
        if(!isClique[i]) continue;
        toDelete.insert(end(toDelete),begin(comps[i]),end(comps[i]));
        vector<int> expandedNodeSet;
        for(auto v : comps[i])
            expandedNodeSet.insert(end(expandedNodeSet), begin(inst.idmap[v]), end(inst.idmap[v]));
        res.done_clusters.push_back(expandedNodeSet);
    }
    res = remove_nodes(res, toDelete);

    return res;
}

std::optional<Instance> weightedKernel(const Instance &inst) {
    int n = size(inst.edges);
    const auto& g = inst.edges;

    for(int v=0; v<n; ++v) {

        vector<int> Nv{v};
        for(int x=0; x<n; ++x)
            if(x!=v && g[v][x]>0)
                Nv.push_back(x);

        if(size(Nv) == 1) continue;

        long long deficiency = 0;
        int cut = 0;
        bool hasZero = false;
        for(int u : Nv) {
            // edges to outside
            for(int x=0; x<n; ++x)
                if(x!=v && g[v][x]<=0) // x is not in Nv
                    if(g[u][x]>0) // x connected to u
                        cut += g[u][x];

            // anti-edges inside
            for(auto w : Nv) {
                if(u<=w) continue; // we only want to look at pairs u>w
                if(g[u][w]==0) hasZero = true;
                if(g[u][w]<0) deficiency += abs(g[u][w]);
            }
        }

        if(hasZero) continue;
        if(deficiency>=INF) continue;
        int pv = 2*deficiency + cut;
        assert(pv<INF);
        if(pv >= size(Nv)) continue;
        // we are reducible

        // merge Nv and make the smallest index in Nv the representative
        auto res = inst;
        sort(begin(Nv), end(Nv));
        int rep = Nv.front();
        while(size(Nv)>1) res = merge(res, rep, Nv.back()), Nv.pop_back();

        int single_connected = -1;
        for(int x=0; x<size(res.edges); ++x) {
            if(x==rep) continue;
            if(res.edges[x][rep]>0) {
                assert(single_connected==-1);
                single_connected = x;
            } else {
                res.edges[rep][x] = res.edges[x][rep] = -INF;
            }
        }

        return res;
    }
    return {};
}

std::optional<Instance> force_small_components(const Instance &inst) {
    int n = inst.edges.size();
    auto compNum = connectedComponents(inst.edges);
    auto numComps = *max_element(begin(compNum), end(compNum)) + 1;
    if(numComps==1) return {};

    vector<vector<int>> nodesPerComp(numComps);
    for(int i=0; i<n; ++i) nodesPerComp[compNum[i]].push_back(i);

    // solve all components that are smaller than 21 nodes
    Instance res = inst;
    vector<int> to_delete;
    int rem_comps = numComps;
    for(int c=0; c<numComps && rem_comps>1; ++c) {
        if(nodesPerComp[c].size()>20) continue;
        rem_comps--;

        // build sub-instance
        Instance comp(nodesPerComp[c].size());
        for (int u = 0; u < size(nodesPerComp[c]); ++u) {
            comp.idmap[u] = inst.idmap[nodesPerComp[c][u]];
            for (int v = 0; v < nodesPerComp[c].size(); ++v) {
                comp.edges[u][v] = inst.edges[nodesPerComp[c][u]][nodesPerComp[c][v]];
                comp.orig[u][v] = inst.orig[nodesPerComp[c][u]][nodesPerComp[c][v]];
            }
        }

        // solve it
        auto s = solve_exact(comp);

        // add clusters to global solution and remove component from instance
        res.spendCost += s.cost;
        for(auto& cluster : s.cliques) res.done_clusters.push_back(cluster);
        for(auto v : nodesPerComp[c]) to_delete.push_back(v);
    }

    if(rem_comps == numComps) return {};

    res = remove_nodes(res, to_delete);
    return res;
}
