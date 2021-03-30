
#include <cluster_editing/exact/star_bound.h>

#include <cluster_editing/exact/reductions.h>
#include <cluster_editing/exact/lower_bounds.h>

#include <cluster_editing/datastructures/adjacency_row.h>
#include <cluster_editing/datastructures/robin_hood.h>

#include <cassert>
#include <utility>
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>
#include <map>
#include <set>
#include <iostream>
#include <optional>

using namespace std;
using namespace cluster_editing;

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

    bool includes_edge(int u, int v) {
        auto has_u = nodes[0]==u || binary_search(begin(nodes)+1, end(nodes), u);
        auto has_v = nodes[0]==v || binary_search(begin(nodes)+1, end(nodes), v);
        return has_u && has_v;
    }
};

namespace std {
    template<>
    struct hash<Star> {
        size_t operator()(const Star &star) const {
            const void *ptr = star.nodes.data();
            size_t len = star.nodes.size() * sizeof(decltype(star.nodes[0]));
            return robin_hood::hash_bytes(ptr, len);
        }
    };
}


class Potential {
    int n;
    Edges potential;
    std::vector<ds::AdjacencyRow> m_non_zero_potential;
    std::vector<ds::AdjacencyRow> m_adj;
    int bound = 0;

public:
    explicit Potential(Edges edges) :
            n(edges.size()), potential(std::move(edges)),
            m_non_zero_potential(potential.size(), ds::AdjacencyRow(potential.size())),
            m_adj(potential.size(), ds::AdjacencyRow(potential.size())) {
        for (int u = 0; u < n; ++u) {
            for (int v = 0; v < n; ++v) {
                assert(-INF <= potential[u][v] && potential[u][v] <= INF);
                if (u == v) {
                    assert(!m_non_zero_potential[u].test(v));
                    continue;
                }
                if (potential[u][v] > 0) {
                    m_non_zero_potential[u].set(v);
                    m_adj[u].set(v);
                } else if (potential[u][v] == 0) {
                    assert(!m_non_zero_potential[u].test(v));
                } else {
                    m_non_zero_potential[u].set(v);
                    assert(!m_adj[u].test(v));
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
                    if (m_non_zero_potential[u][v]) {
                        valid = false;
                    }
                } else if (potential[u][v] > 0) {
                    if (!m_non_zero_potential[u][v]) {
                        valid = false;
                    }
                    if (!m_adj[u].test(v)) {
                        valid = false;
                    }
                    if (!(0 <= potential[u][v] && potential[u][v] <= original_edges[u][v])) {
                        valid = false;
                    }
                } else {
                    if (!m_non_zero_potential[u][v]) {
                        valid = false;
                    }
                    if (m_adj[u].test(v)) {
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

    bool is_consistent(const Edges &original_edges,
                       const robin_hood::unordered_map<Star, int> &stars_in_bound) {
        bool valid = is_consistent(original_edges);
        Potential other(original_edges);
        for (auto[star, weight] : stars_in_bound) {
            other.apply_star_to_potential(star, weight);
        }
        if (potential != other.potential) {
            valid = false;
        }
        if (m_non_zero_potential != other.m_non_zero_potential) {
            valid = false;
        }
        if (m_adj != other.m_adj) {
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
        return !m_non_zero_potential[u].test(v);
    };

    [[nodiscard]] const auto &unused_pairs(int u) const {
        return m_non_zero_potential[u];
    }

    [[nodiscard]] const auto &adj(int u) const {
        assert(0 <= u && u < n);
        return m_adj[u];
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

            for (size_t j = i + 1; j < star.nodes.size(); ++j) {
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

        bound += ((int) star.nodes.size() - 2) * weight;


        /*for (int i = 1; i < star.nodes.size(); ++i) {
            int u = star.nodes[i];
            stars_for_edge[make_pair(star.center(), u)] = star;
        }*/
        for (auto u : star.nodes) {
            for (auto v : star.nodes) {
                if (u == v)
                    continue;
                if (u == star.center() || v == star.center()) {
                    assert(has_positive_potential(u, v));
                    potential[u][v] -= weight;
                    assert(potential[u][v] >= 0);
                    if (potential[u][v] == 0) {
                        assert(u != v);
                        m_non_zero_potential[u][v] = false;
                    }
                } else {
                    assert(has_negative_potential(u, v));
                    potential[u][v] += weight;
                    assert(potential[u][v] <= 0);
                    if (potential[u][v] == 0) {
                        assert(u != v);
                        m_non_zero_potential[u][v] = false;
                    }
                }
            }
        }

        return ((int) star.nodes.size() - 2) * weight;
    }

    int unapply_star_to_potential(Star star, int weight) {
        assert(0 < weight && weight < INF);
        bound -= ((int) star.nodes.size() - 2) * weight;
        for (auto u : star.nodes) {
            for (auto v : star.nodes) {
                if (u == v)
                    continue;
                if (u == star.center() || v == star.center()) {
                    if (potential[u][v] == 0) {
                        assert(u != v);
                        m_non_zero_potential[u][v] = true;
                    }
                    assert(potential[u][v] >= 0);
                    potential[u][v] += weight;
                    assert(potential[u][v] > 0);

                } else {
                    assert(potential[u][v] <= 0);
                    if (potential[u][v] == 0) {
                        assert(u != v);
                        m_non_zero_potential[u][v] = true;
                    }
                    potential[u][v] -= weight;
                    assert(potential[u][v] < 0);
                }
            }
        }
        return -((int) star.nodes.size() - 2) * weight;
    };

    [[nodiscard]] int get_bound() const {
        return bound;
    }

    [[nodiscard]] const Edges &get_potential() const {
        return potential;
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
    std::vector<robin_hood::unordered_set<Star>> stars;
    robin_hood::unordered_map<pair<int, int>, robin_hood::unordered_set<Star>, pair_hash> stars_for_edge;
    robin_hood::unordered_map<Star, int> stars_in_bound;
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

    [[nodiscard]] bool has_star(const Star &star) const {
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

    [[nodiscard]] bool can_add(const std::vector<int> &nodes) const {
        for (auto v1 : nodes)
            for (auto v2 : nodes)
                if (v1 != v2 && potential.pair_used(v1, v2))
                    return false;
        return true;
    };

    [[nodiscard]] bool can_add_candidate(const Candidate &candidate) const {
        assert(std::all_of(candidate.nodes.begin(), candidate.nodes.end(),
                           [&](int u) { return 0 <= u && u < n; }));
        if (candidate.type == P3)
            return can_add(candidate.nodes);

        auto star = Star(vector<int>(candidate.nodes.begin(), candidate.nodes.end() - 1));
        auto ext = candidate.nodes.back();
        if (!has_star(star))
            return false;
        return std::all_of(star.nodes.begin(), star.nodes.end(), [&](int node) { return !potential.pair_used(node, ext); });
    };

    [[nodiscard]] int add_star(Star star, int weight = -1) {
        // `insert` might invalidates iterators to `stars[...]` and `stars_for_edge[...]`
        // `erase` might invalidate iteration over `free_edges_g[...]`
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
        // `insert` might invalidate iteration over `free_edges_g[...]`
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

    [[nodiscard]] std::vector<Candidate> get_candidates(int u, int v) const {
        assert(0 <= u && u < n && 0 <= v && v < n);
        vector<Candidate> star_extensions;

        vector<array<int, 3>> p3s;
        if (potential[u][v] > 0) {

            {
                const auto &v_row_adj = potential.adj(v);
                const auto &v_row_unused = potential.unused_pairs(v);
                const auto &u_row_adj = potential.adj(u);
                const auto &u_row_unused = potential.unused_pairs(u);
                auto get_block_y = [&](auto i) {
                    // assert(potential.free_edges(v) == (v_row_adj.block(i) & v_row_unused.block(i)));
                    return v_row_adj.block(i) & v_row_unused.block(i) & ~u_row_adj.block(i) & u_row_unused.block(i);
                };
                ds::AdjacencyRow::for_each(v_row_adj.num_blocks(), get_block_y, [&](int y) {
                    assert(potential[u][y] < 0);
                    assert(y != u);
                    p3s.push_back({v, u, y});
                });
            }
            {
                const auto &u_row_adj = potential.adj(u);
                const auto &u_row_unused = potential.unused_pairs(u);
                const auto &v_row_adj = potential.adj(v);
                const auto &v_row_unused = potential.unused_pairs(v);
                auto get_block_y = [&](auto i) {
                    // assert(potential.free_edges(u) == (u_row_adj.block(i) & u_row_unused.block(i)));
                    return u_row_adj.block(i) & u_row_unused.block(i) & ~v_row_adj.block(i) & v_row_unused.block(i);
                };
                ds::AdjacencyRow::for_each(u_row_adj.num_blocks(), get_block_y, [&](int y) {
                    assert(potential[v][y] < 0);
                    assert(y != v);
                    p3s.push_back({u, v, y});
                });
            }

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
            {
                const auto &u_row_adj = potential.adj(u);
                const auto &u_row_unused = potential.unused_pairs(u);
                const auto &v_row_adj = potential.adj(v);
                const auto &v_row_unused = potential.unused_pairs(v);
                auto get_block = [&](auto i) {
                    // assert((potential.free_edges(u).block(i) & potential.free_edges(v).block(i)) == (u_row_adj.block(i) & u_row_unused.block(i) & v_row_adj.block(i) & v_row_unused.block(i)));
                    return u_row_adj.block(i) & u_row_unused.block(i) & v_row_adj.block(i) & v_row_unused.block(i);
                };
                ds::AdjacencyRow::for_each(u_row_adj.num_blocks(), get_block, [&](int y) {
                    p3s.push_back({y, v, u});
                });
            }

            for (auto p : {pair{u, v}, {v, u}}) {
                auto ext = p.first;
                auto existing = p.second;
                const auto &ext_row_adj = potential.adj(ext);
                const auto &ext_row_unused = potential.unused_pairs(ext);
                auto get_block_center = [&](auto i) {
                    // assert(potential.free_edges(ext).block(i) == (ext_row_adj.block(i) & ext_row_unused.block(i)));
                    return ext_row_adj.block(i) & ext_row_unused.block(i);
                };
                ds::AdjacencyRow::for_each(ext_row_adj.num_blocks(), get_block_center, [&](auto center) {
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
                });
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
                    auto p = global_star_bound_config.min_degree_candidate_vs_random_candidate_probability;
                    Candidate replacement = (uni_dist(gen) < p)
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

    [[nodiscard]] const auto &get_stars_in_bound() const {
        return stars_in_bound;
    }
};

vector<int> degeneracyOrdering(const robin_hood::unordered_map<int, vector<int>> &g) {
    vector<int> answer;
    robin_hood::unordered_map<int, size_t> deg;
    vector<robin_hood::unordered_set<int>> nodesByDeg;

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

vector<vector<int>> coloring(const robin_hood::unordered_map<int, vector<int>> &g) {
    int maxColor = -1;
    robin_hood::unordered_map<int, int> color;
    vector<vector<int>> ans;

    for (auto u : degeneracyOrdering(g)) {
        set<int> neighborColors;
        if (auto g_u = g.find(u); g_u != g.end()) {
            for (auto v : g_u->second) {
                if (auto color_it = color.find(v); color_it != color.end())
                    neighborColors.insert(color_it->second);
            }
        }
        int c = 0;
        for (int neighborColor : neighborColors) {
            if (c != neighborColor)
                break;
            ++c;
        }
        if (c > maxColor) {
            maxColor = c;
            ans.emplace_back();
        }
        ans[c].push_back(u);
        color[u] = c;
    }
    return ans;
}

auto star_bound_packing(const Instance &inst, int limit) {
    auto bound = StarBound(inst.edges);

    int n = size(inst.edges);

    // Calculate degrees
    auto node_sizes = vector<pair<int, int>>();
    for (int i = 0; i < n; ++i) {
        int num_neighbors = 0;
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            num_neighbors += inst.edges[i][j] > 0;
        }
        node_sizes.emplace_back(num_neighbors, i);
    }

    sort(node_sizes.begin(), node_sizes.end(), greater<>());

    for (auto[_, u] : node_sizes) {
        // Calc neighbor candidates
        const auto &u_row_adj = bound.potential.adj(u);
        const auto &u_row_unused = bound.potential.unused_pairs(u);
        auto get_block_candidate = [&](auto i) {
            // assert(bound.potential.free_edges(u).block(i) == (u_row_adj.block(i) & u_row_unused.block(i)));
            return u_row_adj.block(i) & u_row_unused.block(i);
        };
        // Build neighborhood graph
        robin_hood::unordered_map<int, vector<int>> neighbor_graph;
        ds::AdjacencyRow::for_each(u_row_adj.num_blocks(), get_block_candidate, [&](auto cand1) {
            neighbor_graph[cand1] = vector<int>();
        });
        ds::AdjacencyRow::for_each(u_row_adj.num_blocks(), get_block_candidate, [&](auto cand1) {
            ds::AdjacencyRow::for_each(u_row_adj.num_blocks(), get_block_candidate,[&](auto cand2) {
                if (cand1 != cand2)
                    if (bound.potential.pair_used(cand1, cand2) || bound.potential[cand1][cand2] >= 0) // > or >= ?
                        neighbor_graph[cand1].push_back(cand2);
            });
        });
        // Build stars based on neighborhood coloring
        auto stars = coloring(neighbor_graph);
        for (const auto &star_leaves : stars) {
            if (size(star_leaves) < 2)
                continue;
            vector new_nodes = {u};
            std::copy(star_leaves.begin(), star_leaves.end(), std::back_inserter(new_nodes));
            bound.add_star(Star(new_nodes));
            if (bound.potential.get_bound() > limit) {
                assert(bound.is_consistent());
                return bound;
            }
        }
    }
    int old_bound = 0;
    size_t num_unchanged = 0;
    size_t max_num_unchanged = global_star_bound_config.max_num_unchanged;
    while (old_bound < bound.potential.get_bound() || num_unchanged < max_num_unchanged) {
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
                return bound;
            }
        }

        assert(bound.is_consistent());
    }
    return bound;
}

int star_bound(const Instance &inst, int limit) {
    return star_bound_packing(inst, limit).potential.get_bound();
}

std::optional<Instance> forcedChoicesStarBound(const Instance& inst, int upper_bound, bool verbose) {
    auto bound = star_bound_packing(inst, upper_bound);
    Edges potential = bound.potential.get_potential();

    int n=size(inst.edges);

    vector<tuple<int,int,int>> forced;
    for(int u=0; u<n; ++u) {
        for(int v=u+1; v<n; ++v) {
            if(inst.edges[u][v]==-INF) continue;

            auto old_edge = inst.edges[u][v];
            // we think that not removing the stars overlapping with the modified edge
            // never gives us triples that we are not allowed to take

            // see what happens if i set uv to permanent or forbidden
            for(auto choice : {INF, -INF}) {

                // remove stuff that uses this edge (we omit this and only check weights)
                auto new_bound = bound.potential.get_bound();

                if (!((old_edge > 0 && choice == INF) || (old_edge < 0 && choice == -INF))) {
                    new_bound -= abs(old_edge - potential[u][v]);
                }

                auto uv_cost = potential[u][v];
                if(choice==INF && uv_cost<0) new_bound += uv_cost;
                if(choice==-INF && uv_cost>0) new_bound += uv_cost;
                potential[u][v] = potential[v][u] = choice;
                for(int x=0; x<n; ++x) {
                    auto t = Triple(u, v, x, potential);
                    if (t.valid) new_bound += t.cost;
                }
                potential[u][v] = potential[v][u] = uv_cost;

                // reapply triple structures (we omit this)

                if(new_bound+inst.spendCost > upper_bound)
                    forced.emplace_back(u,v,-choice);
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

std::optional<Instance> forcedChoicesSingleMerge(const Instance& inst, int upper_bound, bool verbose) {

    int n = size(inst.edges);

    // pay stuff is good
    auto similarity = inst.edges;
    for(int u=0; u<n; ++u)
        for(int v=0; v<n; ++v)
            similarity[u][v] = max(0, similarity[u][v]);

    // create new triples is good
    for(int u=0; u<n; ++u)
        for(int v=u+1; v<n; ++v)
            for(int w=0; w<n; ++w)
                if(inst.edges[u][w]>0 && inst.edges[v][w]>0)
                    similarity[u][v] += min(inst.edges[u][w], inst.edges[v][w]);

    // destroy old triples is bad
    for(int u=0; u<n; ++u)
        for(int v=u+1; v<n; ++v)
            for(int w=0; w<n; ++w)
                if(inst.edges[u][w]>0 != inst.edges[v][w]>0) {
                    auto t = Triple(u,v,w,inst.edges);
                    if(t.valid) similarity[u][v] -= t.cost;
                }

    vector<array<int,3>> best5;
    for(int u=0; u<n; ++u)
        for(int v=u+1; v<n; ++v) {
            best5.push_back({similarity[u][v],u,v});
            sort(begin(best5), end(best5), greater<>());
            if(size(best5)>10) best5.pop_back();
        }

    auto previous_lower = star_bound_packing(inst, INF);
    if(verbose) cout << endl << upper_bound -  inst.spendCost - previous_lower.potential.get_bound() << endl;

    auto potential = previous_lower.potential.get_potential();
    for(auto [sim, u,v] : best5) {
        auto subinst = inst;
        subinst.edges[u][v] = subinst.edges[v][u] = -INF;
        subinst.spendCost += max(0,inst.edges[u][v]);
        auto lower = star_bound(subinst, upper_bound);


        auto new_bound = previous_lower.potential.get_bound();
        auto old_edge = inst.edges[u][v];
        if (old_edge >= 0) new_bound -= abs(old_edge - potential[u][v]);
        auto uv_cost = potential[u][v];
        if(uv_cost>0) new_bound += uv_cost;
        potential[u][v] = potential[v][u] = -INF;
        for(int x=0; x<n; ++x) {
            auto t = Triple(u, v, x, potential);
            if (t.valid) new_bound += t.cost;
        }
        potential[u][v] = potential[v][u] = uv_cost;

        if(verbose) cout << "brute-force gap: " << upper_bound - subinst.spendCost - lower;
        if(verbose) cout << "\tgreedy_gap: " << upper_bound - new_bound-inst.spendCost << endl;

        if(subinst.spendCost + lower > upper_bound)
            return merge(inst, u, v);
    }

    return {};
}
