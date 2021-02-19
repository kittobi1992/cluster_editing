
#include <cluster_editing/exact/lower_bounds.h>

#include <cassert>
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>

using namespace std;

struct Triple {
    int u=-1,v=-1,w=-1;
    int cost=0;
    bool valid = false;

    Triple() = default;
    Triple(int a, int b, int c, const Edges& potential) {
        if (a==b || b==c || a==c)
            return;

        if (potential[a][b] < 0)
            u = c, v = a, w = b;
        else if (potential[b][c] < 0)
            u = a, v = b, w = c;
        else if (potential[a][c] < 0)
            u = b, v = a, w = c;

        if (u==-1)
            return;
        if (v > w)
            swap(v, w);
        auto uv = potential[u][v];
        auto vw = potential[v][w];
        auto uw = potential[u][w];
        if (!(uv > 0 && uw > 0))
            return;
        assert(uv > 0);
        assert(uw > 0);
        assert(vw < 0);

        cost = min(uv, min(uw, -vw));
        assert(cost > 0);
        valid = true;
    }

    // returns cost modification of total packing
    int apply(Edges& potential, bool undo=false) {
        int mod = undo ? -cost : cost;
        potential[v][u] = potential[u][v] = potential[u][v] - mod;
        potential[w][u] = potential[u][w] = potential[u][w] - mod;
        potential[w][v] = potential[v][w] = potential[v][w] + mod;
        assert(potential[u][v] >= 0);
        assert(potential[u][w] >= 0);
        assert(potential[v][w] <= 0);
        if(undo) {
            assert(potential[u][v]);
            assert(potential[u][w]);
            assert(potential[v][w]);
        }
        return mod;
    }

    bool operator<(const Triple& rhs) const {
        return !rhs.valid || this->cost < rhs.cost;
    }
};

int packingCost(const vector<Triple> & packing) {
    int cost = 0;
    for (auto & triple : packing)
        cost += triple.cost;

    return cost;
}

void maximizePacking(vector<Triple>& packing, Edges& potential, int limit) {
    int n = size(potential);
    // build initial packing
    auto cost = packingCost(packing);
    for(int u=0; u<n; ++u) {
        for(int v=0; v<n; ++v) {
            auto& uv = potential[u][v];
            for(int w=v+1; w<n && uv>0; ++w) {
                auto& uw = potential[u][w];
                auto& vw = potential[v][w];
                if(uw<=0 || vw>=0) continue;
                Triple t(u,v,w,potential);
                assert(t.valid);
                cost += t.apply(potential);
                packing.push_back(t);
                if(cost>limit) return;
            }
        }
    }
}

void swapLocal(vector<Triple>& packing, Edges& potential, int limit) {
    int n = size(potential);
    mt19937 gen(1337);
    uniform_real_distribution<> dist; // [0,1)


    // local search part
    int num_no_improvement = 0;
    auto cost = packingCost(packing);
    for (int iter = 0; iter < INF && num_no_improvement < 5; ++iter) {

        bool has_improved = false;

        for (int i = 0; i < size(packing); ++i) {
            // try to replace this triple in packing
            auto[u, v, w, uvw_cost, valid] = packing[i];
            // undo triple
            cost += packing[i].apply(potential, true);

            vector<pair<Triple, int>> alternatives; // 0-3 elemns
            for (auto[a, b] : {pair(u, v), pair(u, w), pair(v, w)}) {

                tuple<Triple, double, int> mx;
                auto&[best, rand, cand] = mx;
                for (int x = 0; x < n; ++x) {
                    if (x == u || x == v || x == w) continue;
                    Triple t(a, b, x, potential);
                    if (!t.valid) continue;
                    mx = max(mx, tuple(t, dist(gen), x));
                }
                for (auto&[t, x] : alternatives)
                    best.valid &= cand != x;
                if (best.valid)
                    alternatives.push_back(pair(best, cand));
            }

            int sum = 0;
            for (auto&[t, x] : alternatives) sum += t.cost;

            if (uvw_cost > sum || (uvw_cost == sum && dist(gen) < 0.7)) { // uvw was good
                cost += packing[i].apply(potential); // do not replace uvw
                continue;
            }

            // replace uvw
            has_improved = uvw_cost < sum;
            for (auto&[t, x] : alternatives) cost += t.apply(potential);
            packing[i] = alternatives.back().first;
            alternatives.pop_back();
            for (auto&[t, x] : alternatives) packing.push_back(t);

            if(cost>limit) return;
        }

        num_no_improvement = has_improved ? 0 : num_no_improvement + 1;

    }
}

int packing_local_search_bound(const Instance& inst, int limit) {
    // TODO use limit
    mt19937 gen(random_device{}());
    uniform_real_distribution<> dist; // [0,1)

    // build initial packing
    vector<Triple> packing;
    auto potential = inst.edges;
    maximizePacking(packing, potential,limit);
    swapLocal(packing, potential,limit);
    auto cost = packingCost(packing);

    for(int throws = 0; throws<0; ++throws) {

        vector<Triple> nextPacking;
        nextPacking.reserve(packing.size());
        auto next_potential = inst.edges;
        for(auto& t : packing) if(dist(gen)<0.8) nextPacking.push_back(t), t.apply(next_potential);
        maximizePacking(nextPacking, next_potential, limit);
        swapLocal(nextPacking, next_potential, limit);
        auto next_cost = packingCost(nextPacking);

        if(next_cost > cost) {
            packing = nextPacking;
            potential = next_potential;
            cost = next_cost;
        }
    }

    return cost;
}
