
#include <cluster_editing/exact/lower_bounds.h>

#include <cassert>
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>

using namespace std;

Triple::Triple(int a, int b, int c, const Edges& potential) {
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
int Triple::apply(Edges& potential, bool undo) {
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

    // Return early if the bound is large enough already
    if(cost>limit) return;

    for (int iter = 0; iter < INF && num_no_improvement < 5; ++iter) {

        const int old_cost = cost;

        for (size_t i = 0; i < size(packing); ++i) {
            // try to replace this triple in packing
            auto[u, v, w, uvw_cost, valid] = packing[i];
            // undo triple
            cost += packing[i].apply(potential, true);

            struct Alternative {
                int node;
                int cost;
                double random;
            };
            std::array<vector<Alternative>, 3> alternatives;
            std::array<pair<int, int>, 3> node_pairs = {pair(u, v), pair(u, w), pair(v, w)};
            for (int pi = 0; pi < 3; ++pi) {
                auto[a, b]  = node_pairs[pi];
                for (int x = 0; x < n; ++x) {
                    if (x == u || x == v || x == w) continue;
                    Triple t(a, b, x, potential);
                    if (!t.valid) continue;
                    alternatives[pi].push_back(
                        Alternative{x, t.cost, dist(gen)});
                }

                // Sort dcreasingly by cost
                std::sort(alternatives[pi].begin(), alternatives[pi].end(),
                          [](const Alternative &a, const Alternative &b) {
                              return std::make_tuple(a.cost, a.random) > std::make_tuple(b.cost, b.random);
                });
            }

            int best_cost = -1;
            double best_random = -1;
            Triple best_single_replacement;
            bool two_found = false;
            for (int pi = 0; pi < 3; ++pi) {
                for (auto [x1, cost1, r1] : alternatives[pi]) {
                    Triple t1(node_pairs[pi].first, node_pairs[pi].second, x1, potential);
                    assert(t1.valid);
                    cost += t1.apply(potential);
                    if (std::make_tuple(t1.cost, r1) > std::make_tuple(best_cost, best_random)) {
                        best_single_replacement = t1;
                        best_cost = t1.cost;
                        best_random = r1;
                    }
                    for (int j = pi + 1; j < 3; ++j) {
                        for (auto [x2, cost2, r2] : alternatives[j]) {
                            Triple t2(node_pairs[j].first, node_pairs[j].second, x2, potential);
                            if (!t2.valid) continue;
                            cost += t2.apply(potential);
                            packing.push_back(t1);
                            packing.push_back(t2);
                            two_found = true;
                            break;
                            cost += t2.apply(potential, true);
                        }
                        if (two_found) break;
                    }
                    if (two_found) break;
                    cost += t1.apply(potential, true);
                }
                if (two_found) break;
            }

            if (two_found) {
                for (int pi = 0; pi < 3; ++pi) {
                    for (auto [x1, cost1, r1] : alternatives[pi]) {
                        Triple t1(node_pairs[pi].first, node_pairs[pi].second, x1,
                                  potential);
                        if (!t1.valid)
                            continue;
                        cost += t1.apply(potential);
                        packing.push_back(t1);
                    }
                }
            } else if (best_cost > uvw_cost || (best_cost > 0 && dist(gen) < 0.7)) {
                    cost += best_single_replacement.apply(potential);
                    packing.push_back(best_single_replacement);
            }

            Triple orig(u, v, w, potential);
            if (orig.valid) {
                cost += orig.apply(potential);
                packing[i] = orig;
            } else {
                packing[i] = packing.back();
                packing.pop_back();
            }

            assert(cost == packingCost(packing));

            if(cost>limit) return;
        }

        num_no_improvement = cost > old_cost ? 0 : num_no_improvement + 1;

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

vector<Triple> getAPacking(const Instance &inst) {
    vector<Triple> packing;
    auto potential = inst.edges;
    maximizePacking(packing, potential,INF);
    swapLocal(packing, potential,INF);
    return packing;
}
