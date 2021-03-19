
#include <cluster_editing/exact/kab_bounds.h>

#include <vector>
#include <cassert>
#include <algorithm>

using namespace std;

struct KAB {
    vector<int> A,B;
    bool valid = false;

    KAB(const vector<int>& left, const vector<int>& right, const Edges& potential) : A(left), B(right) {
        for(auto a : A) for(auto b : A) if(a!=b && potential[a][b]>=0) return;
        for(auto a : B) for(auto b : B) if(a!=b && potential[a][b]>=0) return;
        for(auto a : A) for(auto b : B) if(potential[a][b]<=0) return;

        valid = true;
    }

    int apply(Edges& potential, bool undo = false) {
        int mod = undo ? -1 : 1;
        for(auto a : A) for(auto b : A) if(a!=b) potential[a][b] += mod, assert(potential[a][b]<=0);
        for(auto a : B) for(auto b : B) if(a!=b) potential[a][b] += mod, assert(potential[a][b]<=0);
        for(auto a : A) for(auto b : B) potential[b][a] = potential[a][b] -= mod, assert(potential[a][b]>=0);
        int a = size(A), b = size(B);
        if(a>b) swap(a,b);
        return mod * a*(b-1);
    }

    bool canAdd(int v, bool left, Edges& potential) {
        if(find(begin(A),end(A),v) != end(A)) return false;
        if(find(begin(B),end(B),v) != end(B)) return false;
        for(auto a : (left ? A : B)) if(potential[a][v]>=0) return false;
        for(auto b : (left ? B : A)) if(potential[b][v]<=0) return false;
        return true;
    }
};


int kab_bound(const Instance &inst, int limit) {

    // build initial packing
    auto potential = inst.edges;
    int n = size(potential);

    int cost = 0;
    vector<KAB> packing;

    for(int u=0; u<n; ++u) {
        for(int v=0; v<n; ++v) {
            for(int w=v+1; w<n; ++w) {
                if(u==v || u==w) continue;
                auto& uv = potential[u][v];
                auto& uw = potential[u][w];
                auto& vw = potential[v][w];
                if(uv<=0 || uw<=0 || vw>=0) continue;

                KAB kab({u}, {v,w}, potential);
                assert(kab.valid);

                // extent kab as far as possible
                for(int x=0; x<n; ++x) {
                    bool left = size(kab.A) < size(kab.B);
                    if(kab.canAdd(x,left,potential)) (left ? kab.A : kab.B).push_back(x);
                    if(kab.canAdd(x,!left,potential)) (!left ? kab.A : kab.B).push_back(x);
                }

                cost += kab.apply(potential);
                packing.push_back(kab);
                if(cost>limit) return cost;
            }
        }
    }

    return cost;
}
