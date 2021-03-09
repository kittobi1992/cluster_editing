
#include <cluster_editing/exact/star_bound.h>

#include <cassert>
#include <utility>
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>
#include <map>
#include <set>
#include <unordered_set>
#include <iostream>

using namespace std;

template<class T>
static inline void hash_combine(std::size_t &seed, T const &v) {
    // https://www.boost.org/doc/libs/1_75_0/doc/html/hash/reference.html#boost.hash_combine
    seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

struct pair_hash {
    template<class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2> &pair) const {
        std::size_t seed = 0;
        hash_combine(seed, pair.first);
        hash_combine(seed, pair.second);
        return seed;
    }
};

// TODO Put this back into StarBound namespace with hash function?
struct Star {
    std::vector<int> nodes;

    explicit Star(std::vector<int> p_nodes) : nodes(std::move(p_nodes)) {
        assert(nodes.size() >= 3);
        make_canonical();
    }

    void make_canonical() {
        sort(nodes.begin() + 1, nodes.end());
    }

    bool is_canonical() {
        return nodes.size() >= 3 && is_sorted(nodes.begin() + 1, nodes.end());
    }

    [[nodiscard]] int center() const {
        return nodes[0];
    }

    bool operator<(const Star &rhs) const {
        return nodes < rhs.nodes;
    }

    bool operator==(const Star &rhs) const {
        return nodes == rhs.nodes;
    }
};

namespace std {
    template<>
    struct hash<Star> {
        size_t operator()(const Star &star) const {
            std::size_t seed = 0;
            for (auto node : star.nodes)
                hash_combine(seed, node);
            return seed;
        }
    };
}

class StarBound {
protected:
    enum CandidateType {
        P3,
        STAR_EXTENSION
    };


    struct Candidate {
        Candidate(CandidateType p_type, std::vector<int> p_nodes, int p_shared_p3) :
                type(p_type), nodes(std::move(p_nodes)), shared_p3(p_shared_p3) {
            assert(!nodes.empty());
        }

        CandidateType type;
        std::vector<int> nodes;
        int shared_p3;
    };
    //Candidate = collections.namedtuple('Candidate', ['candidate_type',
    //        'nodes', 'shared_p3'])

    Edges g;
    std::vector<std::vector<bool>> used_by_bound;
    std::vector<std::unordered_set<Star>> stars;
    std::vector<std::set<int>> free_edges_g;
    std::unordered_map<pair<int, int>, Star, pair_hash> stars_for_edge;
    std::vector<std::vector<int>> p3_count;

    std::mt19937_64 gen;

public:
    int bound = 0;

    explicit StarBound(Edges g)
            : g(g), used_by_bound(g.size(), vector<bool>(g.size(), false)), stars(g.size()), free_edges_g(g.size()),
              p3_count(g.size(), vector<int>(g.size(), 0)) {
        //TODO Make this weighted somehow?
        for (size_t u = 0; u < g.size(); ++u) {
            for (size_t v = 0; v < g.size(); ++v)
                if (v != u && g[u][v] > 0)
                    free_edges_g[u].insert(v);
        }

        //TODO This is cubic and unweighted
        for (size_t u = 0; u < g.size(); ++u) {
            for (size_t v = 0; v < g.size(); ++v) {
                if (v == u || g[u][v] < 0)
                    continue;
                for (size_t w = 0; w < g.size(); ++w) {
                    if (w >= v || w == u || g[u][w] < 0 || g[v][w] > 0)
                        continue;
                    p3_count[u][v]++;
                    p3_count[v][u]++;
                    p3_count[v][w]++;
                    p3_count[w][v]++;
                    p3_count[u][w]++;
                    p3_count[w][u]++;
                }
            }
        }

    }

    std::set<int> free_neighbors(int u) const {
        assert(0 <= u && u < (int) size(g));
        return free_edges_g[u];
    };

    bool pair_used(int u, int v) const {
        assert(0 <= u && u < (int) size(g) && 0 <= v && v < (int) size(g));
        return used_by_bound[u][v];
    };

    bool has_star(const Star &star) const {
        assert(std::all_of(star.nodes.begin(), star.nodes.end(), [&](int u) { return 0 <= u && u < (int) size(g); }));
        const auto &entry = stars[star.center()];
        return entry.find(star) != entry.end();
    };

    std::vector<Star> stars_in_random_order() {
        std::vector<Star> ans;
        for (auto &u_stars : stars) {
            std::copy(u_stars.begin(), u_stars.end(), std::back_inserter(ans));
        }
        std::shuffle(ans.begin(), ans.end(), gen);
        return ans;

    };

    bool can_add(const Star &star) const {
        for (auto v1 : star.nodes)
            for (auto v2 : star.nodes)
                if (v1 != v2 && pair_used(v1, v2))
                    return false;
        return true;
    };

    bool can_add_candidate(const Candidate &candidate) const {
        assert(std::all_of(candidate.nodes.begin(), candidate.nodes.end(),
                           [&](int u) { return 0 <= u && u < (int) size(g); }));
        if (candidate.type == P3)
            return can_add(Star(candidate.nodes));

        auto star = Star(vector<int>(candidate.nodes.begin(), candidate.nodes.end() - 1));
        auto ext = candidate.nodes.back();
        if (!has_star(star))
            return false;
        for (auto node : star.nodes) {
            if (pair_used(node, ext))
                return false;
        }
        return true;
    };

    void add_star(const Star &star) {
        assert(std::all_of(star.nodes.begin(), star.nodes.end(), [&](int u) { return 0 <= u && u < (int) size(g); }));
        assert(star.nodes.size() > 2);
        for (size_t i = 1; i < star.nodes.size(); ++i) {
            int u = star.nodes[i];
            assert(g[u][star.center()] > 0);
            for (size_t j = i + 1; j < star.nodes.size(); ++j) {
                int v = star.nodes[j];
                assert(g[u][v] < 0);
            }
        }

        stars[star.center()].insert(star);
        bound += static_cast<int>(star.nodes.size() - 2);


        for (size_t i = 1; i < star.nodes.size(); ++i) {
            int u = star.nodes[i];
            stars_for_edge.insert({{star.center(), u}, star});
        }
        for (auto u : star.nodes) {
            for (auto v : star.nodes) {
                if (u == v)
                    continue;
                assert(!pair_used(u, v));
                used_by_bound[u][v] = true;
                if (g[u][v] > 0) {
                    free_edges_g[u].erase(v);
                }
            }
        }
    };

    void remove_star(const Star &star) {
#ifndef NDEBUG
        assert(!star.nodes.empty());
        assert(std::all_of(star.nodes.begin(), star.nodes.end(), [&](int u) { return 0 <= u && u < (int) size(g); }));
        auto star_copy = star;
#endif
        stars[star.center()].erase(star);
        assert(star ==
               star_copy); // Make sure that the input parameter star is not a reference into stars[star.center()].
        bound -= static_cast<int>(star.nodes.size() - 2);
        for (size_t i = 1; i < star.nodes.size(); ++i) {
            int u = star.nodes[i];
            stars_for_edge.erase({star.center(), u});
        }
        assert(star ==
               star_copy); // Make sure that the input parameter star is not a reference into some stars_for_edge[...].
        assert(std::all_of(star.nodes.begin(), star.nodes.end(), [&](int u) { return 0 <= u && u < (int) size(g); }));
        for (auto u : star.nodes) {
            for (auto v : star.nodes) {
                if (u == v)
                    continue;
                used_by_bound[u][v] = false;
                if (g[u][v] > 0)
                    free_edges_g[u].insert(v);
            }
        }
    };

    void add_candidate(const Candidate &candidate) {
        assert(std::all_of(candidate.nodes.begin(), candidate.nodes.end(),
                           [&](int u) { return 0 <= u && u < (int) size(g); }));
        if (candidate.type == CandidateType::P3) {
            add_star(Star(candidate.nodes));
        } else {
            auto existing_star = candidate.nodes;
            existing_star.pop_back();
            remove_star(Star(existing_star));
            add_star(Star(candidate.nodes));
        }
    }

    void remove_candidate(const Candidate &candidate) {
        assert(std::all_of(candidate.nodes.begin(), candidate.nodes.end(),
                           [&](int u) { return 0 <= u && u < (int) size(g); }));
        if (candidate.type == CandidateType::P3) {
            remove_star(Star(candidate.nodes));
        } else {
            auto existing_star = candidate.nodes;
            existing_star.pop_back();
            remove_star(Star(candidate.nodes));
            add_star(Star(existing_star));
        }
    }

    std::vector<Candidate> get_candidates(int u, int v) const {
        assert(0 <= u && u < (int) size(g) && 0 <= v && v < (int) size(g));
        vector<Candidate> star_extensions;

        vector<array<int, 3>> p3s;
        if (g[u][v] > 0) {
            for (int y : free_edges_g[v])
                if (y != u && g[u][y] < 0 && !pair_used(u, y))
                    p3s.push_back({v, u, y});
            for (int y : free_edges_g[u])
                if (y != v && g[v][y] < 0 && !pair_used(v, y))
                    p3s.push_back({u, v, y});

            for (auto[center, ext] : {pair{u, v}, {v, u}}) {
                for (const auto &star : stars[center]) {
                    bool valid = true;
                    for (size_t i = 1; i < star.nodes.size(); ++i) {
                        int x = star.nodes[i];
                        if (pair_used(ext, x) || g[ext][x] >= 0) {
                            valid = false;
                            break;
                        }
                    }
                    if (valid) {
                        auto new_nodes = star.nodes;
                        new_nodes.push_back(ext);
                        int count = 0;
                        for (auto x : star.nodes)
                            count += p3_count[ext][x];
                        star_extensions.emplace_back(STAR_EXTENSION, std::move(new_nodes), count);
                    }
                }
            }
        } else if (g[u][v] < 0) { //TODO or <= 0?
            for (auto y : free_edges_g[u])
                if (free_edges_g[v].find(y) != free_edges_g[v].end())
                    p3s.push_back({y, v, u});


            for (auto[ext, existing] : {pair{u, v}, {v, u}}) {
                for (auto center : free_edges_g[ext]) {
                    if (auto star_it = stars_for_edge.find({center, existing}); star_it != stars_for_edge.end()) {
                        const auto&[_, star] = *star_it;

                        bool valid = true;
                        for (size_t i = 1; i < star.nodes.size(); ++i) {
                            int x = star.nodes[i];
                            if (pair_used(ext, x) || g[ext][x] >= 0) {
                                valid = false;
                                break;
                            }
                        }
                        if (valid) {
                            auto new_nodes = star.nodes;
                            new_nodes.push_back(ext);
                            int count = 0;
                            for (auto x : star.nodes)
                                count += p3_count[ext][x];
                            assert(std::all_of(new_nodes.begin(), new_nodes.end(),
                                               [&](int u) { return 0 <= u && u < (int) size(g); }));
                            star_extensions.emplace_back(STAR_EXTENSION, std::move(new_nodes), count);
                        }
                    }
                }
            }
        }

        for (auto nodes : p3s) {
            int count = 0;
            auto[a, b, c] = nodes;
            for (auto[x, y] : {pair{a, b}, {a, c}, {b, c}})
                count += p3_count[x][y];
            assert(0 <= a && a < (int) size(g));
            assert(0 <= b && b < (int) size(g));
            assert(0 <= c && c < (int) size(g));
            star_extensions.emplace_back(P3, std::vector{a, b, c}, count);
        }
        return star_extensions;
    };

    void try_improve(const Star &star) {

        // First attempt: merge with another star
        auto star2s = stars[star.center()]; // Must be a copy so that remove_star(star2) does not invalidate star2
        for (const auto &star2 : star2s) {
            if (star == star2) {
                continue;
            }
            bool all_unused_non_edges = true;
            for (auto u : star.nodes) {
                if (u == star.center())
                    continue;
                for (auto v : star2.nodes) {
                    if (v == star2.center())
                        continue;
                    // if bound.pair_used(u, v) or v in bound.g[u]:
                    if (pair_used(u, v) || (g[u][v] >= 0)) { // NOTE: > 0?
                        all_unused_non_edges = false;
                        break;
                    }
                }
                if (!all_unused_non_edges)
                    break;
            }
            if (all_unused_non_edges) {
                remove_star(star);
                remove_star(star2);
                Star new_star(star);
                assert(!star2.nodes.empty());
                new_star.nodes.insert(new_star.nodes.end(), star2.nodes.begin() + 1, star2.nodes.end());
                new_star.make_canonical();
                add_star(new_star);
                return;
            }
        }

        Star candidate_star = star;
        std::vector<int> nodes(candidate_star.nodes.begin() + 1, candidate_star.nodes.end());
        for (auto v : nodes) {
            remove_star(candidate_star);

            bool is_p3 = candidate_star.nodes.size() == 3;

            if (!is_p3) {
                assert(candidate_star.nodes.size() > 3);
                candidate_star.nodes.erase(std::remove(candidate_star.nodes.begin(), candidate_star.nodes.end(), v),
                                           candidate_star.nodes.end());
                add_star(candidate_star);
            }

            vector<pair<pair<int, int>, vector<Candidate>>> candidates_per_pair;
            if (is_p3) {
                int a = candidate_star.nodes[0];
                int b = candidate_star.nodes[1];
                int c = candidate_star.nodes[2];
                for (auto[x, y] : {pair{a, b}, {a, c}, {b, c}})
                    candidates_per_pair.emplace_back(pair{x, y}, get_candidates(x, y));
            } else {
                for (auto x : candidate_star.nodes) {
                    assert(g[x][v] < 0 || x == star.center());
                    candidates_per_pair.emplace_back(pair{v, x}, get_candidates(v, x));
                }
            }

            bool two_found = false;
            for (const auto &[pair, candidates] : candidates_per_pair) {
                assert(!two_found);
                for (const auto &candidate : candidates) {
                    assert(!two_found);
                    add_candidate(candidate);

                    for (const auto &[pair2, candidates2] : candidates_per_pair) {
                        if (pair == pair2)
                            continue;
                        for (const auto &candidate2 : candidates2) {
                            if (can_add_candidate(candidate2)) {
                                two_found = true;
                                add_candidate(candidate2);
                            }
                        }
                    }
                    if (two_found)
                        break;
                    remove_candidate(candidate);
                }
                if (two_found)
                    break;
            }

            // Re-insert v
            if (!two_found) {
                // Replace by random candidate
                std::vector<Candidate> all_candidates;
                for (const auto &[pair, candidates] : candidates_per_pair)
                    all_candidates.insert(all_candidates.end(), candidates.begin(), candidates.end());
                if (!all_candidates.empty()) {
                    auto uni_dist = std::uniform_real_distribution<float>(0, 1);
                    auto idx_dist = std::uniform_int_distribution<size_t>(0, all_candidates.size() - 1);

                    const auto cmp = [](const Candidate &a, const Candidate &b) {
                        return a.shared_p3 < b.shared_p3;
                    };
                    Candidate replacement = (uni_dist(gen) < 0.8)
                                            ? *std::min_element(all_candidates.begin(), all_candidates.end(), cmp)
                                            : all_candidates[idx_dist(gen)];

                    auto r_nodes = replacement.nodes;
                    auto s_nodes = candidate_star.nodes;
                    s_nodes.push_back(v);
                    std::sort(r_nodes.begin(), r_nodes.end());
                    std::sort(s_nodes.begin(), s_nodes.end());
                    if (r_nodes == s_nodes) {
                        candidate_star.nodes.push_back(v);
                        candidate_star.make_canonical();
                    }
                    add_candidate(replacement);
                } else {
                    if (!is_p3) {
                        remove_star(candidate_star);
                        candidate_star.nodes.push_back(v);
                        candidate_star.make_canonical();
                    }
                    add_star(candidate_star);
                }
            }
            if (is_p3) {
                break;
            }
        }
    }
};

vector<int> degeneracyOrdering(const map<int, vector<int>> &g) {
    vector<int> answer;
    map<int, size_t> deg;
    vector<set<int>> nodesByDeg;

    for (const auto &[u, neighbors] : g) {
        auto d = deg[u] = size(neighbors);
        nodesByDeg.resize(max<size_t>(nodesByDeg.size(), d + 1));
        nodesByDeg[d].insert(u);
    }
    for (size_t d = 0; d < nodesByDeg.size(); ++d) {
        auto curNodes = vector<int>(nodesByDeg[d].begin(), nodesByDeg[d].end());
        for (size_t i = 0; i < curNodes.size(); ++i) {
            auto u = curNodes[i];
            if (auto g_u = g.find(u); g_u != g.end()) {
                for (auto v : g_u->second) {
                    if (auto dv_it = deg.find(v); dv_it != deg.end()) {
                        auto[_, dv] = *dv_it;
                        if (dv > d) {
                            nodesByDeg[dv].erase(v);
                            if (dv > d + 1)
                                nodesByDeg[dv - 1].insert(v);
                            else
                                curNodes.push_back(v);
                        }
                        deg[v]--;
                    }
                }
            }
            deg.erase(u);
            answer.push_back(u);
        }
    }

    //sort(begin(nodes), end(nodes), [&g](auto& a, auto& b){ return size(g[a])>size(g[b]); });
    reverse(answer.begin(), answer.end());

    return answer;
}

vector<vector<int>> coloring(const map<int, vector<int>> &g) {
    int maxColor = -1;
    auto color = map<int, int>();
    auto ans = vector<vector<int>>();

    for (auto u : degeneracyOrdering(g)) {
        set<int> neighborColors;
        if (auto g_u = g.find(u); g_u != g.end()) {
            for (auto v : g_u->second) {
                if (auto color_it = color.find(v); color_it != color.end())
                    neighborColors.insert(color_it->second);
            }
        }
        int c = 0;
        while (neighborColors.find(c) != neighborColors.end())
            c++;
        if (c > maxColor) {
            maxColor = c;
            ans.emplace_back();
        }
        ans[c].push_back(u);
        color[u] = c;
    }
    return ans;
}

int star_bound(const Instance &inst, int limit) {
    auto potential = inst.edges;
    auto bound = StarBound(potential);

    int n = size(potential);

    // Calculate degrees
    auto node_sizes = vector<pair<int, int>>();
    for (int i = 0; i < n; ++i) {
        int neighbors = 0;
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            neighbors += potential[i][j] > 0;
        }
        node_sizes.emplace_back(neighbors, i);
    }

    sort(node_sizes.begin(), node_sizes.end(), greater<>());

    for (auto[_, u] : node_sizes) {
        // Calc neighbor candidates
        auto candidates = bound.free_neighbors(u);
        // Build neighborhood graph
        auto neighborgraph = map<int, vector<int>>();
        for (auto cand1 : candidates) {
            for (auto cand2 : candidates)
                if (cand1 != cand2)
                    if (bound.pair_used(cand1, cand2) || potential[cand1][cand2] >= 0) // > or >= ?
                        neighborgraph[cand1].push_back(cand2);
        }
        // Build stars based on neighborhood coloring
        auto stars = coloring(neighborgraph);
        for (const auto &star_leaves : stars) {
            if (size(star_leaves) < 2)
                continue;
            auto new_nodes = vector<int>();
            new_nodes.push_back(u);
            std::copy(star_leaves.begin(), star_leaves.end(), std::back_inserter(new_nodes));
            bound.add_star(Star(new_nodes));
            if (bound.bound > limit)
                return bound.bound;
        }
    }
    int old_bound = 0;
    int num_unchanged = 0;
    while (old_bound < bound.bound || num_unchanged < 5) {
        if (old_bound < bound.bound)
            num_unchanged = 0;
        else
            num_unchanged++;

        old_bound = bound.bound;
        //cout << bound.bound << endl;
        for (const auto &star : bound.stars_in_random_order()) {
            // Check if for some reason we have removed this star
            if (!bound.has_star(star))
                continue;

            bound.try_improve(star);
            if (bound.bound > limit)
                return bound.bound;
        }
    }
    return bound.bound;
    /*
    auto potential = inst.edges;
    auto n = size(potential);

    // Calculate degrees
    auto node_sizes = vector<pair<int, int>>();
    for (int i = 0; i < n; ++i) {
        int neighbors = 0;
        for (int j = 0; j < n; ++j) {
            if (i==j) continue;
            neighbors += potential[i][j] > 0;
        }
        node_sizes.push_back({neighbors, i});
    }

    sort(node_sizes.begin(), node_sizes.end(), greater<>());
    
    int bound = 0;
    for (auto [_, u] : node_sizes) {
        // Calc neighbor candidates
        auto candidates = vector<int>();
        for (int v = 0; v < n; ++v)
            if (v != u && potential[u][v] > 0)
                candidates.push_back(v);
        // Build neighborhood graph
        auto neighborgraph = map<int, vector<int>>();
        for (auto cand1 : candidates) {
            for (auto cand2 : candidates)
                if (cand1 != cand2 && potential[cand1][cand2] >= 0) // why >= 0?
                    neighborgraph[cand1].push_back(cand2);
        }
        // Build stars based on neighborhood coloring
        auto stars = coloring(neighborgraph);
        for (auto & star : stars) {
            if (size(star) < 2)
                continue;

            auto minEdge = INF;
            for (auto v1 : star) {
                minEdge = min(minEdge, potential[u][v1]);

                for (auto v2 : star)
                    if (v1 != v2)
                        minEdge = min(minEdge, -potential[v1][v2]);
            }

            bound += (size(star) - 1) * minEdge;
            if (bound > limit)
                return bound;

            for (auto v1 : star) {
                potential[u][v1] -= minEdge;
                potential[v1][u] -= minEdge;

                for (auto v2 : star)
                    if (v1 != v2) {
                        potential[v1][v2] += minEdge;
                    }
            }

        }
    }

    // TODO Local search

    return bound;*/
}
