
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

    std::set<std::pair<int, int>> used_pairs;
    std::vector<std::vector<bool>> used_by_bound;
    std::vector<std::set<Star>> stars;
    std::vector<std::vector<std::vector<Star>>> stars_for_edge;
    std::vector<std::vector<int>> p3_count;

    StarBound(Edges g)
    : stars(g.size())
    {

    }

    std::vector<int> free_neighbors(int u);

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
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(ans.begin(), ans.end(), g);
        return ans;

    };

    bool can_add(const Star & star) {
        for (auto v1 : star.nodes)
            for (auto v2 : star.nodes)
                if (used_pairs.find({v1, v2}) != used_pairs.end())
                    return false;
        return true;
    };

    bool can_add_candidate(Candidate candidate);

    void add_star(Star);
    void remove_star(Star);
    void add_candidate(Candidate);

    void remove_candidate(Candidate);

    std::vector<Candidate> get_candidates(int u, int v);

    void try_improve(Star star) {
        // TODO: implement star merging

        Star candidate_star = star;
        std::vector<int> nodes(candidate_star.begin() + 1 , candidate_star.end());
        for (auto v : nodes) {
            remove_star(candidate);

            bool is_p3 = star.nodes.size() == 3;

            if (!is_p3) {
                assert(star.nodes.size() > 3);
                candidate.nodes.erase(std::remove(candidate.nodes.begin(), candidate.nodes.end(), v), candidate.leaves.end());
                add_star(candidate);
            }

            std::map<std::pair<int, int>, std::vector<Candidate>> candidates_per_pair;
            if (is_p3) {
                int a = candidate.nodes[0];
                int b = candidate.nodes[1];
                int c = candidate.nodes[2];
                candidates_per_pair[{a, b}] = get_candidates(a, b);
                candidates_per_pair[{a, c}] = get_candidates(a, c);
                candidates_per_pair[{b, c}] = get_candidates(b, c);
            } else {
                for (auto x : candidate.nodes) {
                    // assert v not in bound.g[x] or x == star[0]
                    candidates_per_pair[{v, x}] = get_candidates(v, x);
                }
            }

            bool two_found = false;
            for (auto [pair, candidates] : candidates_per_pair) {
                assert(!two_found);

                for (auto candidate : candidates) {

                }
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