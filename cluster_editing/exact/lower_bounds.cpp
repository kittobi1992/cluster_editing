
#include <cluster_editing/exact/lower_bounds.h>

#include <array>
#include <cassert>
#include <vector>
#include <optional>
#include <numeric>
#include <algorithm>

using namespace std;

int packing_local_search_bound(const Instance& inst, int limit) {
    auto edges = inst.edges;
    using Triple = array<int,4>; // array (u,v,w,cost) represents triple v-u-w with cost being min(uv,uw,-vw)

    int n = size(edges);
    long long cost = 0; // maintains sum of packing[i].cost
    auto apply = [&](const Triple& t, bool undo=false) {
        auto [u,v,w,uvw_cost] = t;
        if(undo) uvw_cost = -uvw_cost;
        cost += uvw_cost;
        edges[v][u] = edges[u][v] = edges[u][v] - uvw_cost;
        edges[w][u] = edges[u][w] = edges[u][w] - uvw_cost;
        edges[w][v] = edges[v][w] = edges[v][w] + uvw_cost;

        if(!undo) {
            assert(edges[u][v]>=0);
            assert(edges[u][w]>=0);
            assert(edges[v][w]<=0);
        }
    };

    vector<Triple> packing;
    for(int u=0; u<n; ++u) {
        for(int v=0; v<n; ++v) {
            auto& uv = edges[u][v];
            if(uv<=0) continue;
            for(int w=v+1; w<n && uv>0; ++w) {
                auto& uw = edges[u][w];
                auto& vw = edges[v][w];
                if(uw<=0 || vw>=0) continue;
                packing.push_back({u,v,w, min(uv,min(uw,-vw))});
                assert(edges[u][v]>0 && edges[u][w]>0 && edges[v][w]<0);
                assert(packing.back()[3]>0);
                apply(packing.back());
                if(cost>limit) return cost;
            }
        }
    }

    // local search part
    int num_no_improvement = 0;
    for(int iter=0; iter<INF && cost<=limit && num_no_improvement<5; ++iter) {

        bool has_improved = false;

        for(int i=0; i<size(packing); ++i) {
            // try to replace this triple in packing
            auto [u,v,w,uvw_cost] = packing[i];
            // undo triple
            apply(packing[i], true);
            assert(edges[u][v]>0 && edges[u][w]>0 && edges[v][w]<0);

            optional<Triple> rep_uv;
            for(int x=0; x<n; ++x) { // find other triple using the edge uv
                if(x==w || x==v || x==u) continue;
                if(edges[v][x]==0 || edges[u][x]==0) continue;
                if((edges[v][x]>0)==(edges[u][x]>0)) continue;
                int c = min(edges[u][v], min(abs(edges[u][x]), abs(edges[v][x])));
                Triple t{u,v,x,c};
                if(edges[u][x]<0) swap(t[0],t[1]);
                assert(edges[t[0]][t[1]]>0 && edges[t[0]][t[2]]>0 && edges[t[1]][t[2]]<0);
                if(!rep_uv || rep_uv->at(3)<c) rep_uv = t;
            }

            optional<Triple> rep_uw;
            for(int x=0; x<n; ++x) { // find other triple using the edge uw
                if(x==w || x==v || x==u || (rep_uv && x==rep_uv->at(2))) continue;
                if(edges[w][x]==0 || edges[u][x]==0) continue;
                if((edges[w][x]>0)==(edges[u][x]>0)) continue;
                int c = min(edges[u][w], min(abs(edges[u][x]), abs(edges[w][x])));
                Triple t{u,w,x,c};
                if(edges[u][x]<0) swap(t[0],t[1]);
                assert(edges[t[0]][t[1]]>0 && edges[t[0]][t[2]]>0 && edges[t[1]][t[2]]<0);
                if(!rep_uw || rep_uw->at(3)<c) rep_uw = t;
            }

            optional<Triple> rep_vw;
            for(int x=0; x<n; ++x) { // find other triple using the edge vw
                if(x==w || x==v || x==u) continue;
                if(rep_uv && x==rep_uv->at(2)) continue;
                if(rep_uw && x==rep_uw->at(2)) continue;
                if(edges[w][x]<=0 || edges[v][x]<=0) continue; // both must be present
                int c = min(-edges[v][w], min(edges[v][x], edges[w][x]));
                Triple t{x,v,w,c};
                assert(edges[t[0]][t[1]]>0 && edges[t[0]][t[2]]>0 && edges[t[1]][t[2]]<0);
                if(!rep_vw || rep_vw->at(3)<c) rep_vw = t;
            }

            vector<Triple> reps;
            if(rep_uv) reps.push_back(*rep_uv);
            if(rep_uw) reps.push_back(*rep_uw);
            if(rep_vw) reps.push_back(*rep_vw);
            auto sum = accumulate(begin(reps), end(reps), 0, [](int acc, auto& t){ return acc+t[3]; });
            if(uvw_cost>=sum) { // uvw was good
                apply(packing[i]);
                continue;
            }

            // replace uvw
            has_improved = true;
            sort(begin(reps), end(reps), [](auto& t1, auto& t2){ return t1[3]>t2[3]; });
            packing[i] = reps[0];
            apply(reps[0]);
            for(int k=1; k<size(reps); ++k) {
                packing.push_back(reps[k]);
                apply(reps[k]);
            }
        }

        num_no_improvement = has_improved ? 0 : num_no_improvement+1;

    }

    return cost;
}
