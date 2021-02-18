#include "gmock/gmock.h"

#include <vector>
#include <fstream>

#include <cluster_editing/exact/instance.h>
#include <cluster_editing/data_path.h>
#include <cluster_editing/exact/solver.h>
#include <cluster_editing/exact/reductions.h>

using namespace std;

const map<int,int> optimalSolutions {
        {1, 3},
        {3, 42},
        {5, 46},
        {7, 86},
        {9, 90},
        {11, 81},
        {13, 181},
        {15, 164},
        {17, 236},
        {19, 298},
        {21, 322},
        {23, 281},
        {25, 439},
        {27, 432},
        {29, 509},
        {31, 285},
        {33, 672},
        {35, 385},
        {37, 703},
        {39, 665},
        {41, 184},
        //{43, 899},
        //{45, 1066},
        {47, 749},
        {49, 854},
        //{51, 1219},
        //{53, 1389},
        {55, 1410},
        {57, 122},
        {59, 976},
        {61, 116},
        {63, 665},
        {65, 128},
        {67, 792},
        //{69, 1562},
        {71, 2131},
        {73, 1612},
        {75, 132},
        {77, 78},
        {79, 48},
        //{81, 1882},
        //{83, 2243},
        {85, 1047},
        {87, 3243},
        {89, 534},
        //{91, 2419},
        //{93, 3373},
        {95, 358},
        {97, 95},
        //{99, 644},
        //{101, 3490},
        //{103, 3658},
        {105, 5320},
        {107, 274},
        {109, 2436},
        {111, 2412},
        {113, 999},
        {115, 306},
        {117, 2777},
        {119, 184},
        {121, 742},
        {123, 416},
        {125, 1413},
        {127, 130},
        {129, 2211},
        {131, 386},
        {133, 4836},
        {135, 772},
        {137, 16},
        {139, 5514},
        //{141, 682},
        {143, 1475},
        {145, 1761},
        {147, 3975},
        {149, 5807},
        {151, 1258},
        {153, 6},
        {155, 63},
        {157, 1117},
        {159, 1483},
        {161, 333},
        {163, 1538},
        {165, 6094},
        //{167, 7235},
        //{169, 7864},
        {171, 2901},
        {173, 100},
        {175, 557},
        {177, 4089},
        //{179, 620},
        //{181, 1445},
        //{183, 2609},
        {185, 203},
        //{187, 3654},
        {189, 8563},
        //{191, 15757},
        //{193, 20403},
        //{195, 2068},
        //{197, 0},
        //{199, 0},
};

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

TEST(ExactTest, canLoadGraphs) {
    for(int i=1; i<200; i+=2) {
        auto edges = read_graph(i);
        EXPECT_GT(size(edges), 0);
    }
}

// test is parametrized by instance number (see instantiation below)
struct SolverTest : public testing::TestWithParam<int> { };
TEST_P(SolverTest, canSolve) {
    auto num = GetParam();
    auto edges = read_graph(num);
    auto inst = Instance(size(edges));
    inst.edges = edges;

    auto sol = solve_exact(inst, INF, 500);
    if(!sol.worked) return; // quit if not solve in timelimit
    // check if solution has same value as claimed and is valid partitioning
    verifySolution(edges, sol);
    // check if solution is optimal
    auto it = optimalSolutions.find(num);
    ASSERT_NE(it, end(optimalSolutions));
    ASSERT_EQ(sol.cost, it->second);
}

auto testNameFunctor = [](const auto& info) { return "exact" + to_string(info.param); };
INSTANTIATE_TEST_SUITE_P(SolveIn1s, SolverTest, testing::Range(1,200,2), testNameFunctor);
