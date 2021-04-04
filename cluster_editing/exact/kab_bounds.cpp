
#include <cluster_editing/exact/kab_bounds.h>

#include <vector>
#include <cassert>
#include <algorithm>
#include <random>
#include <set>
#include <chrono>
#include <map>
#include <iostream>

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
    int improvement() const { return improvement(size(A), size(B)); }

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
    uniform_real_distribution<> dist;
    Edges p3_count;

    Kab_bound (const Instance& _inst, int _limit) :
        inst(_inst),
        limit(_limit),
        potential(_inst.edges),
        n(size(_inst.edges)),
        kabs_per_node(n),
        p3_count(n, vector(n,0))
    {
        for(int u=0; u<n; ++u) {
            for(int v=0; v<n; ++v) {
                for(int w=v+1; w<n; ++w) {
                    if(u==v || u==w) continue;
                    auto& uv = potential[u][v];
                    auto& uw = potential[u][w];
                    auto& vw = potential[v][w];
                    if(uv<=0 || uw<=0 || vw>=0) continue;
                    p3_count[u][w]++; p3_count[w][u]++;
                    p3_count[u][v]++; p3_count[v][u]++;
                    p3_count[w][v]++; p3_count[v][w]++;
                }
            }
        }
    }

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

    // rate a set of candidates that can be applied together
    auto score_set(const vector<Candidate>& vec) {
        int gain = 0;
        for(auto& c : vec) gain += c.improvement;

        int covered_edges = 0;
        for(auto& c : vec) {
            if(c.isP3) covered_edges += 2;
            else {
                auto& kab = packing[c.growing_kab_index];
                covered_edges += size(c.insert_into_A ? kab.B : kab.A);
            }
        }

        double overlap = 0;
        for(auto& c : vec) {
            if(c.isP3) overlap += p3_count[c.p3[0]][c.p3[1]] + p3_count[c.p3[0]][c.p3[2]] + p3_count[c.p3[2]][c.p3[1]];
            else {
                auto& kab = packing[c.growing_kab_index];
                for(auto a : kab.A) for(auto b : kab.A) if(a<b) overlap += p3_count[a][b];
                for(auto a : kab.B) for(auto b : kab.B) if(a<b) overlap += p3_count[a][b];
                for(auto a : kab.A) for(auto b : kab.B) overlap += p3_count[a][b];
            }
        }
        overlap *= 0.5 + dist(gen);
        return tuple(gain, -overlap, rand());
    }

    // rate a set of candidates and an associated kab that is inserted with them
    auto score_with_reinsertion(const vector<Candidate>& vec, const KAB& kab) {
        auto score = score_set(vec);
        get<0>(score) += kab.improvement();
        for(auto a : kab.A) for(auto b : kab.A) if(a<b) get<1>(score) -= dist(gen) * p3_count[a][b];
        for(auto a : kab.B) for(auto b : kab.B) if(a<b) get<1>(score) -= dist(gen) * p3_count[a][b];
        for(auto a : kab.A) for(auto b : kab.B) get<1>(score) -= dist(gen) * p3_count[a][b];
        //get<1>(score) -= size(kab.A) * size(kab.B); // covered_edges
        return score;
    }

    // returns best set of candidates that can be applied together
    vector<Candidate> findReplacementSet(const vector<Candidate>& cands) {
        //assert(!empty(cands));
        auto old_cost = cost;

        // keep track of the best set of replacements so far
        vector<Candidate> best_set;
        decltype(score_set(best_set)) best_score;
        auto update_best = [&](const vector<Candidate>& contender) {
            auto new_score = score_set(contender);
            if(new_score>best_score)
                best_set = contender, best_score = new_score;
        };

        // try to pick one
        for(auto c : cands)
            update_best({c});

        // try to pick 2
        for(auto c1 : cands) {
            auto k1 = applyCandidate(c1);
            for(auto c2 : cands) {
                if(!canAdd(c2)) continue;
                // found two candidates that work together
                auto k2 = applyCandidate(c2);

                // fill up with everything else
                vector<Candidate> filled_cands;
                vector<KAB> filled_kabs;
                for(auto c : cands) {
                    if(!canAdd(c)) continue;
                    filled_kabs.push_back(applyCandidate(c));
                    filled_cands.push_back(c);
                }

                // update best
                auto replacement_set = filled_cands;
                replacement_set.push_back(c1);
                replacement_set.push_back(c2);
                update_best(replacement_set);

                // revert the fill-up
                for(int i=0; i<size(filled_cands); ++i)
                    unapplyCandidate(filled_cands[i], filled_kabs[i]);

                unapplyCandidate(c2, k2);
            }
            unapplyCandidate(c1, k1);
        }

        assert(cost == old_cost);
        return best_set;
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

        auto best_cands = findReplacementSet(cands);
        vector<KAB> res;
        for(auto c : best_cands)
            res.push_back(applyCandidate(c));
        return res;
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

        // remove old and remember some costs
        auto cost_before = cost;
        remove(old);
        auto cost_before_wo_old = cost;

        // keep track of best found set (as in findReplacementSet)
        vector<Candidate> best_set;
        KAB reduced_reinsertion = old;
        auto best_score = score_with_reinsertion(best_set, old);
        auto update_best = [&](const KAB& old_wo_x, const vector<Candidate>& best_replacements_for_x) {
            auto new_score = score_with_reinsertion(best_replacements_for_x, old_wo_x);
            if(new_score<=best_score) return;
            best_set = best_replacements_for_x;
            reduced_reinsertion = old_wo_x;
            best_score = new_score;
        };
        update_best(old, {});

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

                // use edges from x to other nodes in old kab that have become free
                vector<Candidate> cands;
                for(auto u : kab.A) {
                    auto cands_for_pair = candidatesForFreedPair(u, x);
                    cands.insert(end(cands), begin(cands_for_pair), end(cands_for_pair));
                }
                for(auto u : kab.B) {
                    auto cands_for_pair = candidatesForFreedPair(u, x);
                    cands.insert(end(cands), begin(cands_for_pair), end(cands_for_pair));
                }

                remove(kab);
                assert(cost==cost_before_wo_old);

                update_best(kab, findReplacementSet(cands));
            }
        }

        // apply the best found replacement set along with the associated reinsertion of the reduced old cap
        // if nothing was found the set is empty and the reinsertion is just the old kab
        assert(!empty(best_set) || reduced_reinsertion.A==old.A);
        apply(reduced_reinsertion);
        vector<KAB> res{reduced_reinsertion};
        for(auto c : best_set)
            res.push_back(applyCandidate(c));
        return res;
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

    int compute_bound(int time_limit, bool verbose) {
        auto t1 = chrono::steady_clock::now();
        auto end_time = chrono::steady_clock::now() + chrono::seconds(time_limit);

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

        int rounds = 0;
        int improvements = 0;
        int num_stall=0;
        while(cost<limit) {
            if(5*improvements<rounds && chrono::steady_clock::now()>end_time)
                break;
            bool improved = false;
            rounds++;

            verify();
            shuffle(begin(packing), end(packing), gen);
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
            verify();

            num_stall = improved ? 0 : num_stall+1;
            improvements += improved;
        }

        auto t2 = chrono::steady_clock::now();

        if(verbose) {
            int fpos = 0;
            int fneg = 0;
            for(int i=0; i<n; ++i)
                for(int j=i+1; j<n; ++j) {
                    fpos += potential[i][j]>0;
                    fneg += potential[i][j]<0;
                }
            map<pair<int,int>, int> sizes;
            for(const auto& kab : packing) {
                int a = size(kab.A), b = size(kab.B);
                sizes[minmax(a,b)]++;
            }
            auto duration = chrono::duration_cast<chrono::duration<double>>(t2-t1).count();
            cout << endl;
            cout << "free edges " << fpos << endl;
            cout << "free non-edges " << fneg << endl;
            cout << "rounds " << rounds << endl;
            cout << "improvements " << improvements << endl;
            cout << "time needed " << duration << endl;
            cout << "time/round " << duration/rounds << endl;
            cout << "packing sizes " << endl;
            for(auto [k,v] : sizes) cout << "\t( " << k.first << " / " << k.second << " ) " << v << endl;
        }

        return cost;
    }

};

int kab_bound(const Instance &inst, int limit, bool verbose, int time_limit) {
    Kab_bound bound(inst,limit);
    return bound.compute_bound(time_limit,verbose);
}
