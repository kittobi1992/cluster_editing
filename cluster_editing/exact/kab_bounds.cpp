
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

    [[nodiscard]] int apply(Edges& potential, bool undo = false) const {
        int mod = undo ? -1 : 1;
        for(auto a : A) for(auto b : A) if(a!=b) potential[a][b] += mod, assert(potential[a][b]<=0);
        for(auto a : B) for(auto b : B) if(a!=b) potential[a][b] += mod, assert(potential[a][b]<=0);
        for(auto a : A) for(auto b : B) potential[b][a] = potential[a][b] -= mod, assert(potential[a][b]>=0);
        int a = size(A), b = size(B);
        if(a>b) swap(a,b);
        return mod * a*(b-1);
    }

    // check if we can grow
    bool canAdd(int v, bool left, const Edges& potential) const {
        if(find(begin(A),end(A),v) != end(A)) return false;
        if(find(begin(B),end(B),v) != end(B)) return false;
        for(auto a : (left ? A : B)) if(potential[a][v]>=0) return false;
        for(auto b : (left ? B : A)) if(potential[b][v]<=0) return false;
        return true;
    }
};

// checks if we can add t to packing; assumes that t is already a valid triple for an empty packing
bool canAdd(const Edges& potential, array<int,3> t) {
    auto [u,v,w] = t;
    assert((potential[u][v]<0) + (potential[u][w]<0) + (potential[v][w]<0) <= 1);
    assert((potential[u][v]>0) + (potential[u][w]>0) + (potential[v][w]>0) <= 2);
    return potential[u][v] && potential[u][w] && potential[v][w];
}

// returns nodes x that form a P3 with u and v
vector<int> reps4edge(const Edges& potential, int u, int v) {
    vector<int> res;
    int n = size(potential);
    int pos_edge = potential[u][v]>0;
    int neg_edge = potential[u][v]<0;
    for(int x=0; x<n; ++x) {
        if(x==u || x==v) continue;
        if(pos_edge + (potential[x][v]>0) + (potential[x][u]>0) != 2) continue;
        if(neg_edge + (potential[x][v]<0) + (potential[x][u]<0) != 1) continue;
        res.push_back(x);
    }
    return res;
}

auto toKAB(array<int,3> t, const Edges& potential) {
    auto [u,v,w] = t;
    vector<int> A{u}, B;
    (potential[u][v]>0 ? B : A).push_back(v);
    (potential[u][w]>0 ? B : A).push_back(w);
    assert(KAB(A,B,potential).valid);
    return KAB(A,B,potential);
}

vector<KAB> mutate_P3(const KAB& old, Edges& potential, int& cost) {

    cost += old.apply(potential,true);

    // find replacements
    int u,v,w;
    if(size(old.A)==1) u = old.A[0], v = old.B[0], w = old.B[1];
    else u = old.B[0], v = old.A[0], w = old.A[1];
    vector<array<int,3>> reps;
    for(auto [x,y] : {pair(u,v), pair(u,w), pair(v,w)})
        for(auto z : reps4edge(potential, x,y))
            reps.push_back({x,y,z});

    // try to pick 2
    for(auto r1 : reps) {
        auto k1 = toKAB(r1,potential);
        cost += k1.apply(potential);
        for(auto r2 : reps) {
            if(!canAdd(potential, r2)) continue;
            auto k2 = toKAB(r2,potential);
            cost += k2.apply(potential);
            for(auto r3 : reps) {
                if(!canAdd(potential,r3)) continue;
                auto k3 = toKAB(r3,potential);
                cost += k3.apply(potential);
                return {k1,k2,k3};
            }
            return {k1,k2};
        }
        cost += k1.apply(potential,true);
    }

    auto score = [&](array<int,3> t) {
        auto overlap = abs(potential[t[0]][t[1]])
                       + abs(potential[t[0]][t[2]])
                       + abs(potential[t[1]][t[2]]);
        auto rnd = rand();
        return pair(overlap, rnd);
    };

    array<int,3> best{u,v,w};
    auto val = score(best);
    for(auto t : reps) if(score(t)>val) best = t;

    auto k = toKAB(best,potential);
    cost += k.apply(potential);
    return {k};
}

vector<KAB> mutate(const KAB& old, Edges& potential, int& cost) {

    // try to grow
    for(int x=0; x<size(potential); ++x) {
        for(bool left : {0,1}) {
            if(!old.canAdd(x,left,potential)) continue;
            auto bigger = old;
            (left ? bigger.A : bigger.B).push_back(x);
            cost += old.apply(potential, true);
            cost += bigger.apply(potential);
            return {bigger};
        }
    }

    // handle p3 separately
    if(size(old.A) + size(old.B)==3) return mutate_P3(old, potential, cost);

    auto cost_before = cost;
    cost += old.apply(potential, true);
    auto cost_before_wo_old = cost;

    // try to take as many p3s in as possible after we remove one node x from old to yield kab
    auto test = [&](const KAB& kab, int x) {
        auto free_edges = kab.A;
        copy(begin(kab.B), end(kab.B), back_inserter(free_edges));
        vector<KAB> reps;
        for(auto u : free_edges) {
            int best = -1, val = -1;
            for(auto rep : reps4edge(potential, u, x)) {
                auto score = rand();
                if(score>val) best = rep, val = score;

            }
            if(best != -1) {
                auto k = toKAB({x,u,best}, potential);
                cost += k.apply(potential);
                reps.push_back(k);
            }
        }
        return reps;
    };

    // take away from A
    if(size(old.A)>1) {
        for(int i=0; i<size(old.A); ++i) {
            auto A_cpy = old.A;
            swap(A_cpy.back(), A_cpy[i]);
            A_cpy.pop_back();
            auto kab = KAB(A_cpy, old.B, potential);
            cost += kab.apply(potential);
            auto reps = test(kab, old.A[i]);
            reps.push_back(kab);
            if(cost > cost_before || (cost_before == cost && rand() > RAND_MAX/2))
                return reps;
            for(auto& rep : reps) cost += rep.apply(potential,true);
            assert(cost==cost_before_wo_old);
        }
    }

    // take away from B
    if(size(old.B)>1) {
        for(int i=0; i<size(old.B); ++i) {
            auto B_cpy = old.B;
            swap(B_cpy.back(), B_cpy[i]);
            B_cpy.pop_back();
            auto kab = KAB(B_cpy, old.A, potential);
            cost += kab.apply(potential);
            auto reps = test(kab, old.B[i]);
            reps.push_back(kab);
            if(cost > cost_before || (cost_before == cost && rand() > RAND_MAX/2))
                return reps;
            for(auto& rep : reps) cost += rep.apply(potential,true);
            assert(cost==cost_before_wo_old);
        }
    }

    cost += old.apply(potential);
    return {old};
}


int kab_bound(const Instance &inst, int limit) {

    // build initial packing
    auto potential = inst.edges;
    int n = size(potential);

    int cost = 0;
    vector<KAB> packing;

    // start with simple P3 packing
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

                cost += kab.apply(potential);
                packing.push_back(kab);
                if(cost>limit) return cost;
            }
        }
    }

    int num_stall=0;
    while(num_stall<5) {
        bool improved = false;

        vector<KAB> nextPacking;
        for(auto& kab : packing) {
            int old_cost = cost;
            auto reps = mutate(kab, potential, cost);
            if(cost>old_cost) improved = true;
            copy(begin(reps), end(reps), back_inserter(nextPacking));
        }
        packing = nextPacking;

        num_stall = improved ? 0 : num_stall+1;
    }

    return cost;
}
