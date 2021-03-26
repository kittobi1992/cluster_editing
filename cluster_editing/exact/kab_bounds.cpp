
#include <cluster_editing/exact/kab_bounds.h>

#include <vector>
#include <cassert>
#include <algorithm>
#include <random>
#include <set>

using namespace std;

struct KAB {
    vector<int> A,B;
    mutable bool is_in_packing = false;

    // construct KAB from p3
    KAB(const vector<int>& a, const vector<int>& b) : A(a), B(b) {};
    KAB(array<int,3> t, const Edges& potential) {
        auto [u,v,w] = t;
        A.push_back(u);
        (potential[u][v]>0 ? B : A).push_back(v);
        (potential[u][w]>0 ? B : A).push_back(w);
        for(auto a : A) for(auto b : A) assert(a==b || potential[a][b]<0);
        for(auto a : B) for(auto b : B) assert(a==b || potential[a][b]<0);
        for(auto a : A) for(auto b : B) assert(potential[a][b]>0);
    }

    static int improvement(int a, int b) {
        if(a>b) swap(a,b);
        return a*(b-1);
    }

    // check if we can grow
    bool canGrow(int v, bool addToA, const Edges& potential) const {
        if(find(begin(A),end(A),v) != end(A)) return false;
        if(find(begin(B),end(B),v) != end(B)) return false;
        for(auto a : (addToA ? A : B)) if(potential[a][v]>=0) return false;
        for(auto b : (addToA ? B : A)) if(potential[b][v]<=0) return false;
        return true;
    }
};

struct Candidate {
    bool isP3; // p3 or kab growth
    int improvement;

    array<int,3> p3;
    int growth_vertex;
    int growing_kab_index;
    bool insert_into_A;

    bool operator<(const Candidate& rhs) const { return improvement < rhs.improvement; }
};

struct Kab_bound {

    const Instance& inst;
    int limit;
    Edges potential;
    int n;
    vector<KAB> packing;
    int cost = 0;
    vector<vector<int>> kabs_per_node;
    mt19937 gen;

    Kab_bound (const Instance& _inst, int _limit) :
        inst(_inst),
        limit(_limit),
        potential(_inst.edges),
        n(size(_inst.edges)),
        kabs_per_node(n)
    {}

    int packing_size = 0;
    void apply(const KAB& kab, bool undo = false) {
        int mod = undo ? -1 : 1;
        assert(undo == kab.is_in_packing);
        kab.is_in_packing = !undo;
        packing_size += mod;
        for(auto a : kab.A) for(auto b : kab.A) if(a!=b) potential[a][b] += mod, assert(potential[a][b]<=0);
        for(auto a : kab.B) for(auto b : kab.B) if(a!=b) potential[a][b] += mod, assert(potential[a][b]<=0);
        for(auto a : kab.A) for(auto b : kab.B) potential[b][a] = potential[a][b] -= mod, assert(potential[a][b]>=0);
        cost += mod * KAB::improvement(size(kab.A), size(kab.B));
    }

    void remove(const KAB& kab) { apply(kab, true); }

    // checks if we can add t to packing; assumes that t is already a valid triple for an empty packing
    bool canAdd(array<int,3> t) {
        auto [u,v,w] = t;
        assert((inst.edges[u][v]<0) + (inst.edges[u][w]<0) + (inst.edges[v][w]<0) == 1);
        assert((inst.edges[u][v]>0) + (inst.edges[u][w]>0) + (inst.edges[v][w]>0) == 2);
        return potential[u][v] && potential[u][w] && potential[v][w];
    }

    bool canAdd(Candidate c) {
        if(c.isP3) return canAdd(c.p3);
        const KAB& kab = packing[c.growing_kab_index];
        return kab.is_in_packing && kab.canGrow(c.growth_vertex, c.insert_into_A, potential);
    }

    KAB applyCandidate(Candidate c) {
        assert(canAdd(c));
        if(c.isP3) {
            KAB kab = KAB(c.p3,potential);
            apply(kab);
            return kab;
        }

        // we are an extension
        auto& old_kab = packing[c.growing_kab_index];
        KAB new_kab(old_kab.A, old_kab.B);
        remove(old_kab);
        assert(new_kab.canGrow(c.growth_vertex, c.insert_into_A, potential));
        (c.insert_into_A ? new_kab.A : new_kab.B).push_back(c.growth_vertex);
        apply(new_kab);
        return new_kab;
    }

    void unapplyCandidate(Candidate c, const KAB& created_kab) {
        remove(created_kab);
        if(c.isP3) {
            assert(canAdd(c));
            return;
        }
        auto& old_kab = packing[c.growing_kab_index];
        apply(old_kab);
        assert(canAdd(c));
    }

    // returns nodes x that form a P3 with u and v
    vector<Candidate> candidatesForFreedPair(int u, int v) {
        vector<Candidate> res;

        // find p3s
        for(int x=0; x<n; ++x) {
            if(x==u || x==v) continue;
            if((potential[u][v]>0) + (potential[x][v]>0) + (potential[x][u]>0) != 2) continue;
            if((potential[u][v]<0) + (potential[x][v]<0) + (potential[x][u]<0) != 1) continue;
            res.push_back({.isP3=true, .improvement=1, .p3={u,v,x}});
        }

        // find extensions
        for(auto [growing, other] : {pair(u,v), pair(v,u)}) {
            for(auto kab_id : kabs_per_node[other]) {
                const KAB& kab = packing[kab_id];
                if(!kab.is_in_packing) continue;
                auto canLeft = potential[kab.A[0]][growing] < 0;
                if(kab.canGrow(growing,canLeft,potential)) {
                    Candidate c{};
                    c.isP3 = false;
                    c.improvement = KAB::improvement(size(kab.A) + canLeft, size(kab.B) + !canLeft);
                    c.growth_vertex = growing;
                    c.growing_kab_index = kab_id;
                    c.insert_into_A = canLeft;
                    res.push_back(c);
                }
            }
        }

        return res;
    }

    vector<KAB> mutate_P3(const KAB& old) {

        remove(old);

        // find replacements
        int u,v,w;
        if(size(old.A)==1) u = old.A[0], v = old.B[0], w = old.B[1];
        else u = old.B[0], v = old.A[0], w = old.A[1];
        vector<Candidate> cands;
        for(auto [x,y] : {pair(u,v), pair(u,w), pair(v,w)})
            for(auto c : candidatesForFreedPair(x,y))
                cands.push_back(c);

        // if we have one candidate with improvement>1 take it
        assert(size(cands)>1);
        auto mx_cand = *max_element(begin(cands), end(cands));
        if(mx_cand.improvement>1)
            return {applyCandidate(mx_cand)};

        // try to pick 2
        for(auto c1 : cands) {
            auto k1 = applyCandidate(c1);
            for(auto c2 : cands) {
                if(!canAdd(c2)) continue;
                // found two candidates that work together
                auto k2 = applyCandidate(c2);

                // fill up with everything else
                vector<KAB> result{k1,k2};
                for(auto c : cands)
                    if(canAdd(c))
                        result.push_back(applyCandidate(c));
                return result;
            }
            unapplyCandidate(c1, k1);
        }

        auto score = [&](const Candidate& c) {
            return rand();
        };

        assert(!empty(cands));
        auto best = pair(score(cands.front()), cands.front());
        for(auto c : cands)
            best = max(best, pair(score(c), c));

        return {applyCandidate(best.second)};
    }

    vector<KAB> mutate(const KAB& old) {
        assert(old.is_in_packing);

        // try to grow
        for(int x=0; x<n; ++x) {
            for(bool left : {0,1}) {
                if(!old.canGrow(x,left,potential)) continue;
                KAB bigger(old.A, old.B);
                (left ? bigger.A : bigger.B).push_back(x);
                remove(old);
                apply(bigger);
                return {bigger};
            }
        }

        // handle p3 separately
        if(size(old.A) + size(old.B)==3) return mutate_P3(old);

        auto cost_before = cost;
        remove(old);
        auto cost_before_wo_old = cost;

        // take away from A or B of old kab
        for(auto remove_from_A : {0,1}) {
            auto sideRemove = (remove_from_A ? old.A : old.B);
            auto sideKeep = (!remove_from_A ? old.A : old.B);
            if(size(sideRemove)==1) continue;
            for(auto x : sideRemove) {

                // insert old KAB without the node x
                auto side_removed = sideRemove;
                auto it = find(begin(side_removed), end(side_removed), x);
                assert(it!=end(side_removed));
                side_removed.erase(it);
                KAB kab(side_removed, sideKeep);
                apply(kab);

                // fill in edges from x to other nodes in old kab that have become free
                vector<pair<Candidate,KAB>> cands;
                auto free_edges = kab.A; // kab.A ∪ kab.B
                copy(begin(kab.B), end(kab.B), back_inserter(free_edges));
                shuffle(begin(free_edges), end(free_edges), gen);
                for(auto u : free_edges) {
                    Candidate best;
                    int val = -1;
                    for(auto c : candidatesForFreedPair(u, x)) {
                        assert(canAdd(c));
                        auto score = rand();
                        if(score>val) best = c, val = score;

                    }
                    if(val != -1) {
                        auto new_kab = applyCandidate(best);
                        cands.push_back(pair(best,new_kab));
                    }
                }


                //
                if(cost > cost_before || (cost_before == cost && rand() > RAND_MAX/2)) {
                    assert(kab.is_in_packing);
                    vector<KAB> result{kab}; // the old kab without x
                    for(auto [c,new_kab] : cands)
                        assert(new_kab.is_in_packing), result.push_back(new_kab);
                    return result;
                }

                // undo everything
                remove(kab);

                for(auto& [c,new_kab] : cands)
                    unapplyCandidate(c,new_kab);
                assert(cost==cost_before_wo_old);
            }
        }

        KAB same_old(old.A, old.B);
        apply(same_old);
        assert(cost==cost_before);
        return {same_old};
    }

    void constructKABperNodeStructure() {
        for(auto& vec : kabs_per_node)
            vec.clear();
        for(int i=0; i<size(packing); ++i) {
            const auto& kab = packing[i];
            assert(kab.is_in_packing);
            for(auto v : kab.A) kabs_per_node[v].push_back(i);
            for(auto v : kab.B) kabs_per_node[v].push_back(i);
        }
    }

    void verify() {
        int real_cost = 0;
        for(auto& kab : packing) {
            real_cost += KAB::improvement(size(kab.A), size(kab.B));
            set<int> nodes(begin(kab.A), end(kab.A));
            nodes.insert(begin(kab.B), end(kab.B));
            assert(size(nodes)==size(kab.A) + size(kab.B));
            assert(kab.is_in_packing);
        }
        assert(packing_size == size(packing));
        assert(real_cost == cost);
    }

    int compute_bound() {

        // start with simple P3 packing
        for(int u=0; u<n; ++u) {
            for(int v=0; v<n; ++v) {
                for(int w=v+1; w<n; ++w) {
                    if(u==v || u==w) continue;
                    auto& uv = potential[u][v];
                    auto& uw = potential[u][w];
                    auto& vw = potential[v][w];
                    if(uv<=0 || uw<=0 || vw>=0) continue;

                    packing.emplace_back(vector{u}, vector{v,w});
                    apply(packing.back());
                    if(cost>limit) return cost;
                }
            }
        }

        int num_stall=0;
        while(num_stall<5 && cost<=limit) {
            bool improved = false;

            verify();
            constructKABperNodeStructure();
            vector<KAB> nextPacking;
            for(auto& kab : packing) {
                if(!kab.is_in_packing) continue;
                int old_cost = cost;
                auto reps = mutate(kab);
                assert(!kab.is_in_packing);
                for(auto& e : reps) assert(e.is_in_packing);
                if(cost>old_cost) improved = true;
                copy(begin(reps), end(reps), back_inserter(nextPacking));
            }
            assert(size(nextPacking) == packing_size);
            packing = nextPacking;
            for(auto& kab : packing)
                assert(kab.is_in_packing);

            num_stall = improved ? 0 : num_stall+1;
        }
        verify();

        return cost;
    }

};

int kab_bound(const Instance &inst, int limit) {
    Kab_bound bound(inst,limit);
    return bound.compute_bound();
}
