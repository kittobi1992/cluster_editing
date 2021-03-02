
#include <cluster_editing/exact/star_bound.h>

#include <cassert>
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>
#include <map>
#include <set>

using namespace std;

class StarBound {
    struct Star {
        std::vector<int> nodes;

        Star(std::vector<int> nodes)
        : nodes(nodes)
        {
            sort(nodes.begin()+1, nodes.end());
        }

        int center() const {
            return nodes[0];
        }
    };

    enum CandidateType {
        P3,
        STAR_EXTENSION
    };


    struct Candidate {
        CandidateType type;
        std::vector<int> nodes;
        int shared_p3;
    };
    //Candidate = collections.namedtuple('Candidate', ['candidate_type',
    //        'nodes', 'shared_p3'])

    int bound = 0;
    Edges g;
    std::set<std::pair<int, int>> used_pairs;
    std::vector<std::vector<bool>> used_by_bound;
    std::vector<std::set<Star>> stars;
    std::vector<std::set<int>> free_edges_g;
    std::map<pair<int, int>, Star> stars_for_edge;
    std::vector<std::vector<int>> p3_count;

    std::mt19937_64 gen;

    StarBound(Edges g)
    : g(g)
    , stars(g.size())
    , free_edges_g(g.size())
    , p3_count(g.size(), vector<int>(g.size(), 0))
    {
        //TODO Make this weighted somehow?
        for (int u = 0; u < g.size(); ++u) {
            for (int v = 0; v < g.size(); ++v)
                if (v != u && g[u][v] > 0)
                    free_edges_g[u].insert(v);
        }

        //TODO This is cubic and unweighted
        for (int u = 0; u < g.size(); ++u) {
            for (int v = 0; v < g.size(); ++v) {
                if (v==u || g[u][v] < 0)
                    continue;
                for (int w = 0; w < g.size(); ++w) {
                    if (w>=v || w==u || g[u][w] < 0 || g[v][w] > 0)
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

    std::set<int> free_neighbors(int u) {
        return free_edges_g[u];
    };

    bool pair_used(int u, int v) {
        return used_pairs.find({u, v}) != used_pairs.end();
    };

    bool has_star(const Star & star) {
        return stars[star.center()].find(star) != stars[star.center()].end();
    };

    std::vector<Star> stars_in_random_order() {
        std::vector<Star> ans;
        for (auto & u_stars : stars) {
            std::copy(u_stars.begin(), u_stars.end(), std::back_inserter(ans));
        }
        std::shuffle(ans.begin(), ans.end(), gen);
        return ans;

    };

    bool can_add(const Star & star) {
        for (auto v1 : star.nodes)
            for (auto v2 : star.nodes)
                if (used_pairs.find({v1, v2}) != used_pairs.end())
                    return false;
        return true;
    };

    bool can_add_candidate(Candidate candidate) {
        if (candidate.type == P3)
            return can_add(Star(candidate.nodes));

        auto star = Star(vector<int>(candidate.nodes.begin(), candidate.nodes.end()-1));
        auto ext = candidate.nodes.back();
        if (!has_star(star))
            return false;
        for (auto node : star.nodes) {
            if (used_pairs.find({node, ext}) != used_pairs.end())
                return false;
        }
        return true;
    };

    void add_star(Star star) {
        assert(star.nodes.size() > 2);
        for (int i = 1; i < star.nodes.size(); ++i) {
            int u = star.nodes[i];
            assert(g[u][star.center()] > 0);
            for (int j = i+1; j < star.nodes.size(); ++j) {
                int v = star.nodes[j];
                assert(g[u][v] < 0);
            }
        }

        stars[star.center()].insert(star);
        bound += star.nodes.size() - 2;


        for (int i = 1; i < star.nodes.size(); ++i) {
            int u = star.nodes[i];
            stars_for_edge[{star.center(), u}] = star;
        }
        for (auto u : star.nodes) {
            for (auto v : star.nodes) {
                if (u==v)
                    continue;
                assert(used_pairs.find({u, v}) == used_pairs.end());
                used_pairs.insert({u, v});
                if (g[u][v] > 0)
                    free_edges_g[u].erase(v);
            }
        }
    };

    void remove_star(Star star) {
        stars[star.center()].erase(star);
        bound -= star.nodes.size() - 2;
        for (int i = 1; i < star.nodes.size(); ++i) {
            int u = star.nodes[i];
            stars_for_edge.erase({star.center(), u});
        }
        for (auto u : star.nodes) {
            for (auto v : star.nodes) {
                if (u == v)
                    continue;
                used_pairs.erase({u, v});
                if (g[u][v] > 0)
                    free_edges_g[u].insert(v);
            }
        }
    };

    void add_candidate(const Candidate &candidate) {
        if (candidate.type == CandidateType::P3) {
            add_star(Star(candidate.nodes));
        } else {
            auto existing_star = candidate.nodes;
            existing_star.pop_back();
            remove_star(Star(candidate.nodes));
            add_star(Star(existing_star));
        }
    }

    void remove_candidate(const Candidate &candidate) {
        if (candidate.type == CandidateType::P3) {
            remove_star(Star(candidate.nodes));
        } else {
            auto existing_star = candidate.nodes;
            existing_star.pop_back();
            remove_star(Star(candidate.nodes));
            add_star(Star(existing_star));
        }
    }

    std::vector<Candidate> get_candidates(int u, int v) {
        vector<Candidate> star_extensions;

        vector<array<int, 3>> p3s;
        if (g[u][v] > 0) {
            for (int y : free_edges_g[v])
                if (y != u && g[u][y] < 0 && used_pairs.find({u, y}) == used_pairs.end())
                    p3s.push_back({v, u, y});
            for (int y : free_edges_g[u])
                if (y != v && g[v][y] < 0 && used_pairs.find({v, y}) == used_pairs.end())
                    p3s.push_back({u, v, y});

            for (auto [center, ext] : vector<pair<int, int>>({{u, v}, {v, u}})) {
                for (auto star : stars[center]) {
                    bool valid = true;
                    for (int i = 1; i < star.nodes.size(); ++i) {
                        int x = star.nodes[i];
                        if (pair_used(ext, x) || g[ext][x] > 0) {
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
                        star_extensions.push_back(Candidate{STAR_EXTENSION, new_nodes, count});
                    }
                }
            }
        } else if (g[u][v] < 0) { //TODO or <= 0?
            for (auto y : free_edges_g[u])
                if (free_edges_g[v].find(y) != free_edges_g[v].end())
                    p3s.push_back({y, v, u});


            for (auto [ext, existing] : vector<pair<int, int>>({{u, v}, {v, u}})) {
                for (auto center : free_edges_g[ext]) {
                    if (stars_for_edge.find({center, existing}) != stars_for_edge.end()) {
                        auto star = stars_for_edge[{center, existing}];

                        bool valid = true;
                        for (int i = 1; i < star.nodes.size(); ++i) {
                            int x = star.nodes[i];
                            if (pair_used(ext, x) || g[ext][x] > 0) {
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
                            star_extensions.push_back(Candidate{STAR_EXTENSION, new_nodes, count});
                        }
                    }
                }
            }
        }

        for (auto nodes : p3s) {
            int count = 0;
            auto [a, b, c] = nodes;
            for (auto [x, y] : vector<pair<int, int >>{{a, b}, {a, c}, {b, c}})
                count += p3_count[x][y];
            star_extensions.push_back(Candidate{P3, {a,b,c}, count});
        }
        return star_extensions;
    };

    void try_improve(Star star) {
        // TODO: implement star merging

        Star candidate_star = star;
        std::vector<int> nodes(candidate_star.nodes.begin() + 1 , candidate_star.nodes.end());
        for (auto v : nodes) {
            remove_star(candidate_star);

            bool is_p3 = star.nodes.size() == 3;

            if (!is_p3) {
                assert(star.nodes.size() > 3);
                candidate_star.nodes.erase(std::remove(candidate_star.nodes.begin(), candidate_star.nodes.end(), v), candidate_star.nodes.end());
                add_star(candidate_star);
            }

            std::map<std::pair<int, int>, std::vector<Candidate>> candidates_per_pair;
            if (is_p3) {
                int a = candidate_star.nodes[0];
                int b = candidate_star.nodes[1];
                int c = candidate_star.nodes[2];
                candidates_per_pair[{a, b}] = get_candidates(a, b);
                candidates_per_pair[{a, c}] = get_candidates(a, c);
                candidates_per_pair[{b, c}] = get_candidates(b, c);
            } else {
                for (auto x : candidate_star.nodes) {
                    // assert v not in bound.g[x] or x == star[0]
                    candidates_per_pair[{v, x}] = get_candidates(v, x);
                }
            }

            bool two_found = false;
            for (const auto &[pair, candidates] : candidates_per_pair) {
                assert(!two_found);
                for (auto candidate : candidates) {
                    assert(!two_found);
                    add_candidate(candidate);

                    for (const auto &[pair2, candidates2] : candidates_per_pair) {
                        if (pair == pair2)
                            continue;
                        for (auto candidate2 : candidates2) {
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
                for (auto [pair, candidates] : candidates_per_pair)
                    all_candidates.insert(all_candidates.end(), candidates.begin(), candidates.end());
                if (!all_candidates.empty()) {
                    std::uniform_real_distribution<float> dist(0, 1);
                    auto idx_dist = std::uniform_int_distribution<size_t>(0, all_candidates.size() - 1);

                    Candidate replacement = (dist(gen) < 0.8)
                            ? *std::min_element(all_candidates.begin(), all_candidates.end(), [](const Candidate &a, const Candidate &b) {
                            return a.shared_p3 < b.shared_p3;
                                })
                            : all_candidates[idx_dist(gen)];

                    auto r_nodes = replacement.nodes;
                    auto s_nodes = candidate_star.nodes;
                    s_nodes.push_back(v);
                    std::sort(r_nodes.begin(), r_nodes.end());
                    std::sort(s_nodes.begin(), s_nodes.end());
                    if (r_nodes == s_nodes) {
                        candidate_star.nodes.push_back(v);
                    }
                } else {
                    if (!is_p3) {
                        remove_star(candidate_star);
                        candidate_star.nodes.push_back(v);
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

vector<int> degeneracyOrdering(map<int, vector<int>> g) {
    vector<int> answer;
    map<int, int> deg;
    vector<set<int>> nodesByDeg;

    for (auto & [u, neighbors] : g) {
        auto d = deg[u] = size(neighbors);
        while (nodesByDeg.size() <= d)
            nodesByDeg.push_back({});
        nodesByDeg[d].insert(u);
    }
    for (int d = 0; d < nodesByDeg.size(); ++d) {
        auto curNodes = vector<int>(nodesByDeg[d].begin(), nodesByDeg[d].end());
        for (int i = 0; i < curNodes.size(); ++i) {
            auto u = curNodes[i];
            for (auto v : g[u]) {
                if (deg.find(v) != deg.end()) {
                    auto dv = deg[v];
                    if (dv > d) {
                        nodesByDeg[dv].erase(v);
                        if (dv > d+1)
                            nodesByDeg[dv-1].insert(v);
                        else
                            curNodes.push_back(v);
                    }
                    deg[v]--;
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

vector<vector<int>> coloring(map<int, vector<int>> g) {
    int maxColor = -1;
    auto color = map<int, int>();
    auto ans = vector<vector<int>>();

    for (auto u : degeneracyOrdering(g)) {
        set<int> neighborColors;
        for (auto v : g[u]) {
            if (color.find(v) != color.end())
                neighborColors.insert(color[v]);
        }
        int c = 0;
        while (neighborColors.find(c) != neighborColors.end())
            c++;
        if (c > maxColor) {
            maxColor = c;
            ans.push_back({});
        }
        ans[c].push_back(u);
        color[u] = c;
    }
    return ans;
}

int star_bound(const Instance& inst, int limit) {
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

    return bound;
}