#include "gmock/gmock.h"

#include <vector>
#include <fstream>

#include <cluster_editing/exact/instance.h>
#include <cluster_editing/data_path.h>
#include <cluster_editing/exact/solver.h>
#include <cluster_editing/exact/reductions.h>

using namespace std;

vector<vector<int>> read_graph(int num) {
    auto suf = to_string(num);
    while(size(suf)<3) suf = '0'+suf;
    auto file_name = EXACT_DATA_DIR + ("exact" + suf + ".gr");
    ifstream in(file_name);

    std::istringstream sstream;
    auto getline = [&]() {
        std::string line;
        do {
            std::getline(in, line);
        } while (line[0] == 'c');
        sstream = std::istringstream(line);
    };

    getline();
    std::string skip;
    int n, m;
    sstream >> skip >> skip >> n >> m;

    vector<vector<int>> res(n, vector<int>(n, -1));

    for (int i = 0; i < m; ++i) {
        getline();
        int u, v;
        sstream >> u >> v;
        --u;
        --v;
        res[u][v] = 1;
        res[v][u] = 1;
    }

    return res;
}

void verifySolution(const Edges& edges, const Solution& sol) {
    if(!sol.worked) return;

    int n = size(edges);
    vector<int> cluster_id(n,-1);
    for(int i=0; i<size(sol.cliques); ++i) {
        for(auto v : sol.cliques[i]) {
            ASSERT_EQ(cluster_id[v],-1);
            cluster_id[v] = i;
        }
    }
    for(auto id : cluster_id)
        ASSERT_NE(id, -1);
    long long cost = 0;
    for(auto u=0; u<n; ++u)
        for(int v=u+1; v<n; ++v)
            if((cluster_id[u]==cluster_id[v]) == (edges[u][v]<0))
                cost += abs(edges[u][v]);
    ASSERT_EQ(cost, sol.cost);
}

class SolverTest : public testing::TestWithParam<int> {
public:
    const map<int,int> optimalSolutions {
            { 1, 3},
            { 65, 128},
            { 77, 78},
            { 97, 95},
            { 137, 16},
            { 153, 6},
            { 155, 63},
            { 173, 100}
    };
};

TEST(ExactTest, canLoadGraphs) {
    for(int i=1; i<200; i+=2) {
        auto edges = read_graph(i);
        EXPECT_GT(size(edges), 0);
    }
}

// test is parametrized by instance number (see instantiation below)
TEST_P(SolverTest, canSolve) {
    auto num = GetParam();
    auto edges = read_graph(num);
    auto inst = Instance(size(edges));
    inst.edges = edges;

    auto sol = solve_exact(inst);
    // check if solution has same value as claimed and is valid partitioning
    verifySolution(edges, sol);
    // check if solution is optimal
    auto it = optimalSolutions.find(num);
    ASSERT_NE(it, end(optimalSolutions));
    EXPECT_EQ(sol.cost, it->second);
}

auto testNameFunctor = [](const auto& info) { return "exact" + to_string(info.param); };
INSTANTIATE_TEST_SUITE_P(EasyInstances, SolverTest, testing::Values(1, 137, 155, 173), testNameFunctor);
#ifdef NDEBUG
INSTANTIATE_TEST_SUITE_P(MediumInstances, SolverTest, testing::Values(65, 77, 97, 153), testNameFunctor);
#endif

