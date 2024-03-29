/*******************************************************************************
 * This file is part of KaPoCE.
 *
 * KaPoCE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaPoCE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaPoCE.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "gmock/gmock.h"

#include <vector>
#include <fstream>

#include <cluster_editing/exact/instance.h>
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
        auto inst = load_exact_instance(i);
        EXPECT_GT(size(inst.edges), 0);
    }
}

// test is parametrized by instance number (see instantiation below)
struct SolverTest : public testing::TestWithParam<int> { };
TEST_P(SolverTest, canSolve) {
    auto num = GetParam();
    auto inst = load_exact_instance(num);

    auto sol = solve_exact(inst, INF, 500);
    if(!sol.worked) return; // quit if not solve in timelimit
    // check if solution has same value as claimed and is valid partitioning
    verifySolution(inst.edges, sol);
    // check if solution is optimal
    auto it = optimalSolutions.find(num);
    ASSERT_NE(it, end(optimalSolutions));
    ASSERT_EQ(sol.cost, it->second);
}

auto testNameFunctor = [](const auto& info) { return "exact" + to_string(info.param); };
INSTANTIATE_TEST_SUITE_P(SolveIn1s, SolverTest, testing::Range(1,200,2), testNameFunctor);


TEST(CompleTwin, edgeCase) {
    Instance inst(3);
    inst.edges = {
            {-1,1,-INF},
            {1,-1,1},
            {-INF,1,-1},
    };
    inst.orig = {
            {-1,1,-1},
            {1,-1,1},
            {-1,1,-1},
    };
    if(auto opt=complexTwin(inst,true); opt)
        ASSERT_LT(opt->spendCost, INF);
}