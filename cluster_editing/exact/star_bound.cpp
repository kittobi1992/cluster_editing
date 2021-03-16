
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

    template<class It>
    class Range {
        It m_begin, m_end;
    public:
        Range(It begin, It end) : m_begin(begin), m_end(end) {}

        [[nodiscard]] auto begin() const { return m_begin; }

        [[nodiscard]] auto end() const { return m_end; }
    };

    [[nodiscard]] auto leaves() const { return Range(nodes.begin() + 1, nodes.end()); }

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


class Potential {
    int n;
    Edges potential;
    std::vector<std::vector<bool>> potential_is_zero;
    std::vector<std::set<int>> free_edges_g;
    int bound = 0;

public:
    explicit Potential(Edges edges) : n(edges.size()), potential(std::move(edges)), potential_is_zero(potential.size(), vector<bool>(potential.size(), false)), free_edges_g(potential.size()) {
        for (int u = 0; u < n; ++u) {
            for (int v = 0; v < n; ++v) {
                assert(-INF <= potential[u][v] && potential[u][v] <= INF);
                if (u == v)
                    continue;
                if (potential[u][v] > 0) {
                    free_edges_g[u].insert(v);
                } else if (potential[u][v] == 0) {
                    potential_is_zero[u][v] = true;
                }
            }
        }
    }

    [[nodiscard]] bool is_consistent(const Edges &original_edges) const {
        bool valid = true;
        for (int u = 0; u < n; ++u) {
            for (int v = 0; v < n; ++v) {
                if (u == v)
                    continue;
                if (potential[u][v] == 0) {
                    if (!potential_is_zero[u][v]) {
                        valid = false;
                    }
                    if (free_edges_g[u].count(v) == 1) {
                        valid = false;
                    }
                } else if (potential[u][v] > 0) {
                    if (potential_is_zero[u][v]) {
                        valid = false;
                    }
                    if (free_edges_g[u].count(v) == 0) {
                        valid = false;
                    }
                    if (!(0 <= potential[u][v] && potential[u][v] <= original_edges[u][v])) {
                        valid = false;
                    }
                } else {
                    if (potential_is_zero[u][v]) {
                        valid = false;
                    }
                    if (free_edges_g[u].count(v) == 1) {
                        valid = false;
                    }
                    if (!(original_edges[u][v] <= potential[u][v] && potential[u][v] <= 0)) {
                        valid = false;
                    }
                }
            }
        }
        return valid;
    }

    bool is_consistent(const Edges &original_edges, const std::unordered_map<Star, int> &stars_in_bound) {
        bool valid = is_consistent(original_edges);
        Potential other(original_edges);
        for (auto [star, weight] : stars_in_bound) {
            other.apply_star_to_potential(star, weight);
        }
        if (potential != other.potential) {
            valid = false;
        }
        if (free_edges_g != other.free_edges_g) {
            valid = false;
        }
        if (potential_is_zero != other.potential_is_zero) {
            valid = false;
        }
        if (bound != other.bound) {
            valid = false;
        }
        return valid;
    }

    const auto &operator[](int u) const {
        return potential[u];
    }

    [[nodiscard]] bool has_positive_potential(int u, int v) const {
        assert(0 <= u && u < n && 0 <= v && v < n && u != v);
        return potential[u][v] > 0;
    };

    [[nodiscard]] bool has_negative_potential(int u, int v) const {
        assert(0 <= u && u < n && 0 <= v && v < n && u != v);
        return potential[u][v] < 0;
    };

    [[nodiscard]] bool pair_used(int u, int v) const {
        assert(0 <= u && u < n && 0 <= v && v < n && u != v);
        return potential_is_zero[u][v];
    };

    [[nodiscard]] const auto &free_edges(int u) const {
        assert(0 <= u && u < n);
        return free_edges_g[u];
    }

    [[nodiscard]] int calculate_costs(const Star &star) const {
        int weight = INF;
        if (star.nodes.size() <= 2)
            return 0;

        for (size_t i = 1; i < star.nodes.size(); ++i) {
            int u = star.nodes[i];
            if (potential[star.center()][u] <= 0)
                return 0;
            weight = min(weight, potential[star.center()][u]);

            for (size_t j = i+1; j < star.nodes.size(); ++j) {
                int v = star.nodes[j];
                if (potential[u][v] >= 0)
                    return 0;
                weight = min(weight, -potential[u][v]);
            }
        }

        assert(0 < weight && weight < INF);
        return weight;
    }

    int apply_star_to_potential(Star star, int weight) {
#ifndef NDEBUG
        assert(0 < weight && weight < INF);
        assert(std::all_of(star.nodes.begin(), star.nodes.end(), [&](int u) { return 0 <= u && u < n; }));
        assert(star.nodes.size() > 2);
        for (size_t i = 1; i < star.nodes.size(); ++i) {
            int u = star.nodes[i];
            assert(potential[u][star.center()] > 0);
            for (size_t j = i + 1; j < star.nodes.size(); ++j) {
                int v = star.nodes[j];
                assert(u != v);
                assert(potential[u][v] < 0);
            }
            assert(u != star.center());
        }
#endif

        bound += ((int)star.nodes.size() - 2) * weight;


        /*for (int i = 1; i < star.nodes.size(); ++i) {
            int u = star.nodes[i];
            stars_for_edge[make_pair(star.center(), u)] = star;
        }*/
        for (auto u : star.nodes) {
            for (auto v : star.nodes) {
                if (u==v)
                    continue;
                if (u == star.center() || v == star.center()) {
                    assert(has_positive_potential(u, v));
                    potential[u][v] -= weight;
                    assert(potential[u][v] >= 0);
                    if (potential[u][v] == 0) {
                        potential_is_zero[u][v] = true;
                        free_edges_g[u].erase(v);
                    }
                }
                else {
                    assert(has_negative_potential(u, v));
                    potential[u][v] += weight;
                    assert(potential[u][v] <= 0);
                    if (potential[u][v] == 0) {
                        potential_is_zero[u][v] = true;
                    }
                }
            }
        }

        return ((int)star.nodes.size() - 2) * weight;
    }

    int unapply_star_to_potential(Star star, int weight) {
        assert(0 < weight && weight < INF);
        bound -= ((int)star.nodes.size() - 2) * weight;
        for (auto u : star.nodes) {
            for (auto v : star.nodes) {
                if (u == v)
                    continue;
                if (u == star.center() || v == star.center()) {
                    if (potential[u][v] == 0) {
                        potential_is_zero[u][v] = false;
                        free_edges_g[u].insert(v);
                    }
                    assert(potential[u][v] >= 0);
                    potential[u][v] += weight;
                    assert(potential[u][v] > 0);

                }
                else {
                    assert(potential[u][v] <= 0);
                    if (potential[u][v] == 0) {
                        potential_is_zero[u][v] = false;
                    }
                    potential[u][v] -= weight;
                    assert(potential[u][v] < 0);
                }
            }
        }
        return -((int)star.nodes.size() - 2) * weight;
    };

    [[nodiscard]] int get_bound() const {
        return bound;
    }
};


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

    int n;
    Edges original_edges; // only for debugging purposes
public:
    Potential potential;
private:
    std::vector<std::unordered_set<Star>> stars;
    std::unordered_map<pair<int, int>, std::set<Star>, pair_hash> stars_for_edge;
    std::unordered_map<Star, int> stars_in_bound;
    std::vector<std::vector<int>> p3_count;

    std::mt19937_64 gen;

public:
    explicit StarBound(Edges g)
            : n(g.size()), original_edges(g), potential(g), stars(g.size()),
              p3_count(g.size(), vector<int>(g.size(), 0)) {
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

    bool has_star(const Star &star) const {
        assert(std::all_of(star.nodes.begin(), star.nodes.end(), [&](int u) { return 0 <= u && u < n; }));
        return stars[star.center()].find(star) != stars[star.center()].end();
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
                if (v1 != v2 && potential.pair_used(v1, v2))
                    return false;
        return true;
    };

    bool can_add_candidate(const Candidate &candidate) const {
        assert(std::all_of(candidate.nodes.begin(), candidate.nodes.end(),
                           [&](int u) { return 0 <= u && u < n; }));
        if (candidate.type == P3)
            return can_add(Star(candidate.nodes));

        auto star = Star(vector<int>(candidate.nodes.begin(), candidate.nodes.end() - 1));
        auto ext = candidate.nodes.back();
        if (!has_star(star))
            return false;
        return std::all_of(star.nodes.begin(), star.nodes.end(), [&](int node) { return !potential.pair_used(node, ext); });
    };

    [[nodiscard]] int add_star(Star star, int weight = -1) {
        // `insert` might invalidates iterators to `stars[...]` and `stars_for_edge[...]`
        // `erase` might invalidate iterators to `free_edges_g[...]`
        // `stars_in_bound[star]` might invalidate iterators to `stars_in_bound`
#ifndef NDEBUG
        assert(std::all_of(star.nodes.begin(), star.nodes.end(), [&](int u) { return 0 <= u && u < n; }));
        assert(star.nodes.size() > 2);
        for (size_t i = 1; i < star.nodes.size(); ++i) {
            int u = star.nodes[i];
            assert(potential[u][star.center()] > 0);
            for (size_t j = i + 1; j < star.nodes.size(); ++j) {
                int v = star.nodes[j];
                assert(u != v);
                assert(potential[u][v] < 0);
            }
            assert(u != star.center());
        }
#endif

        if (weight == -1)
            weight = potential.calculate_costs(star);
        assert(0 < weight && weight < INF);

        potential.apply_star_to_potential(star, weight);
        stars_in_bound[star] += weight;

        // Skip other data structures if star is already in bound
        if (stars_in_bound[star] > weight) {
            //assert(potential.is_consistent(original_edges, stars_in_bound));
            return weight;
        }

        stars[star.center()].insert(star);

        for (auto u : star.leaves()) {
            stars_for_edge[{star.center(), u}].insert(star);
        }

        //assert(potential.is_consistent(original_edges, stars_in_bound));
        return weight;
    };

    void remove_star(Star star, int weight) {
        // `erase` might invalidates iterators to `stars[...]` and `stars_for_edge[...]`
        // `insert` might invalidate iterators to `free_edges_g[...]`
#ifndef NDEBUG
        assert(0 < weight && weight < INF);
        for (size_t i = 1; i < star.nodes.size(); ++i) {
            int u = star.nodes[i];
            for (size_t j = i + 1; j < star.nodes.size(); ++j) {
                int v = star.nodes[j];
                assert(u != v);
            }
            assert(u != star.center());
        }
#endif

        potential.unapply_star_to_potential(star, weight);
        stars_in_bound[star] -= weight;
        assert(stars_in_bound[star] >= 0);

        // Only touch other data structures if star is now fully removed
        if (stars_in_bound[star] > 0) {
            //assert(potential.is_consistent(original_edges, stars_in_bound));
            return;
        }

        stars_in_bound.erase(star);
        stars[star.center()].erase(star);

        //bound -= static_cast<int>(star.nodes.size() - 2);
        for (auto u : star.leaves()) {
            stars_for_edge[{star.center(), u}].erase(star);
        }
        //assert(potential.is_consistent(original_edges, stars_in_bound));
    };

    [[nodiscard]] int add_candidate(const Candidate &candidate) {
        assert(std::all_of(candidate.nodes.begin(), candidate.nodes.end(),
                           [&](int u) { return 0 <= u && u < n; }));
        if (candidate.type == CandidateType::P3) {
            return add_star(Star(candidate.nodes));
        } else {
            auto existing_star = candidate.nodes;
            existing_star.pop_back();
            remove_star(Star(existing_star), 1);
            return add_star(Star(candidate.nodes));
        }
    }

    void remove_candidate(const Candidate &candidate, int weight) {
        assert(std::all_of(candidate.nodes.begin(), candidate.nodes.end(),
                           [&](int u) { return 0 <= u && u < n; }));
        if (candidate.type == CandidateType::P3) {
            remove_star(Star(candidate.nodes), weight);
        } else {
            auto existing_star = candidate.nodes;
            existing_star.pop_back();
            remove_star(Star(candidate.nodes), weight);
            add_star(Star(existing_star), 1);
        }
    }

    std::vector<Candidate> get_candidates(int u, int v) const {
        assert(0 <= u && u < n && 0 <= v && v < n);
        vector<Candidate> star_extensions;

        vector<array<int, 3>> p3s;
        if (potential[u][v] > 0) {
            for (int y : potential.free_edges(v))
                if (y != u && potential[u][y] < 0 && !potential.pair_used(u, y))
                    p3s.push_back({v, u, y});
            for (int y : potential.free_edges(u))
                if (y != v && potential[v][y] < 0 && !potential.pair_used(v, y))
                    p3s.push_back({u, v, y});

            for (auto p : {pair{u, v}, {v, u}}) {
                auto center = p.first;
                auto ext = p.second;
                for (const auto &star : stars[center]) {
                    auto valid = std::all_of(star.leaves().begin(), star.leaves().end(), [&](auto x) {
                        return ext != x && potential[ext][x] < 0;
                    });
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
        } else if (potential[u][v] < 0) { //TODO or <= 0?
            for (auto y : potential.free_edges(u))
                if (potential.free_edges(v).find(y) != potential.free_edges(v).end())
                    p3s.push_back({y, v, u});


            for (auto p : {pair{u, v}, {v, u}}) {
                auto ext = p.first;
                auto existing = p.second;
                for (auto center : potential.free_edges(ext)) {
                    if (auto stars_it = stars_for_edge.find({center, existing}); stars_it != stars_for_edge.end()) {
                        const auto&[_, stars] = *stars_it;

                        for (const auto &star : stars) { // TODO: Check if that's the correct generalization for multiple stars per edge
                            bool valid = std::all_of(star.leaves().begin(), star.leaves().end(), [&](int x) {
                                return ext != x && potential[ext][x] < 0;
                            });
                            if (valid) {
                                auto new_nodes = star.nodes;
                                new_nodes.push_back(ext);
                                int count = 0;
                                for (auto x : star.nodes)
                                    count += p3_count[ext][x];
                                assert(std::all_of(new_nodes.begin(), new_nodes.end(),
                                                   [&](int u) { return 0 <= u && u < n; }));
                                star_extensions.emplace_back(STAR_EXTENSION, std::move(new_nodes), count);
                            }
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
            assert(0 <= a && a < n);
            assert(0 <= b && b < n);
            assert(0 <= c && c < n);
            assert(potential[a][b] > 0);
            assert(potential[a][c] > 0);
            assert(potential[b][c] < 0);
            star_extensions.emplace_back(P3, std::vector{a, b, c}, count);
            assert(potential.calculate_costs(Star(star_extensions.back().nodes)) > 0);
        }
        return star_extensions;
    };

    void try_improve(const Star &star) {
        //assert(potential.is_consistent(original_edges, stars_in_bound));
        // First attempt: merge with another star
        // Must be a copy because remove_star(star2) might invalidate iterators to stars[star.center()]

        auto star2s = stars[star.center()];
        for (const auto &star2 : star2s) {
            if (star == star2) {
                continue;
            }
            bool all_unused_non_edges = true;
            for (auto u : star.leaves()) {
                for (auto v : star2.leaves()) {
                    // if bound.pair_used(u, v) or v in bound.g[u]:
                    if (u == v || !(potential[u][v] < 0)) { // NOTE: > 0?
                        all_unused_non_edges = false;
                        break;
                    }
                }
                if (!all_unused_non_edges)
                    break;
            }
            if (all_unused_non_edges) {
                remove_star(star, 1);
                remove_star(star2, 1);
                Star new_star(star);
                assert(!star2.nodes.empty());
                new_star.nodes.insert(new_star.nodes.end(), star2.leaves().begin(), star2.leaves().end());
                new_star.make_canonical();
                add_star(new_star);
                return;
            }
        }


        Star candidate_star = star;
        std::vector<int> nodes(candidate_star.leaves().begin(), candidate_star.leaves().end());
        for (auto v : nodes) {
            remove_star(candidate_star, 1);

            bool is_p3 = candidate_star.nodes.size() == 3;

            if (!is_p3) {
                assert(candidate_star.nodes.size() > 3);
                candidate_star.nodes.erase(std::remove(candidate_star.nodes.begin(), candidate_star.nodes.end(), v),
                                           candidate_star.nodes.end());
                assert(potential.calculate_costs(candidate_star) > 0);
                add_star(candidate_star, 1);
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
                    assert(potential[x][v] < 0 || x == star.center());
                    candidates_per_pair.emplace_back(pair{v, x}, get_candidates(v, x));
                }
            }

            bool two_found = false;
            for (const auto &[pair, candidates] : candidates_per_pair) {
                assert(!two_found);
                for (const auto &candidate : candidates) {
                    assert(!two_found);
                    int weight = add_candidate(candidate);

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
                    remove_candidate(candidate, weight);
                }
                if (two_found)
                    break;
            }

            // Re-insert v
            if (!two_found) {
                // Replace by random candidate
                std::vector<Candidate> all_candidates;
                for (auto &[pair, candidates] : candidates_per_pair) {
                    std::move(candidates.begin(), candidates.end(), back_inserter(all_candidates));
                    candidates.clear();
                }
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
                        remove_star(candidate_star, 1);
                        candidate_star.nodes.push_back(v);
                        candidate_star.make_canonical();
                    }
                    assert(potential.calculate_costs(candidate_star) > 0);
                    add_star(candidate_star, 1);
                }
            }
            if (two_found || is_p3) {
                break;
            }
        }
    }

    bool is_consistent() {
        return potential.is_consistent(original_edges, stars_in_bound);
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
    auto bound = StarBound(inst.edges);

    int n = size(inst.edges);

    // Calculate degrees
    auto node_sizes = vector<pair<int, int>>();
    for (int i = 0; i < n; ++i) {
        int neighbors = 0;
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            neighbors += inst.edges[i][j] > 0;
        }
        node_sizes.emplace_back(neighbors, i);
    }

    sort(node_sizes.begin(), node_sizes.end(), greater<>());

    for (auto[_, u] : node_sizes) {
        // Calc neighbor candidates
        auto candidates = bound.potential.free_edges(u);
        // Build neighborhood graph
        auto neighborgraph = map<int, vector<int>>();
        for (auto cand1 : candidates) {
            for (auto cand2 : candidates)
                if (cand1 != cand2)
                    if (bound.potential.pair_used(cand1, cand2) || bound.potential[cand1][cand2] >= 0) // > or >= ?
                        neighborgraph[cand1].push_back(cand2);
        }
        // Build stars based on neighborhood coloring
        auto stars = coloring(neighborgraph);
        for (const auto &star_leaves : stars) {
            if (size(star_leaves) < 2)
                continue;
            vector new_nodes = {u};
            std::copy(star_leaves.begin(), star_leaves.end(), std::back_inserter(new_nodes));
            bound.add_star(Star(new_nodes));
            if (bound.potential.get_bound() > limit) {
                assert(bound.is_consistent());
                return bound.potential.get_bound();
            }
        }
    }
    int old_bound = 0;
    int num_unchanged = 0;
    while (old_bound < bound.potential.get_bound() || num_unchanged < 5) {
        if (old_bound < bound.potential.get_bound())
            num_unchanged = 0;
        else
            num_unchanged++;

        old_bound = bound.potential.get_bound();
        //cout << bound.bound << endl;
        auto stars = bound.stars_in_random_order();
        for (const auto &star : stars) {
            // Check if for some reason we have removed this star
            if (!bound.has_star(star))
                continue;

            bound.try_improve(star);
            if (bound.potential.get_bound() > limit) {
                assert(bound.is_consistent());
                return bound.potential.get_bound();
            }
        }

        assert(bound.is_consistent());
    }
    return bound.potential.get_bound();
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
