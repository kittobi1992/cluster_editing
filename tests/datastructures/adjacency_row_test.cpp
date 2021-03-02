#include "gmock/gmock.h"

#include <chrono>
#include <random>

#include "cluster_editing/datastructures/adjacency_row.h"

using ::testing::Test;

namespace cluster_editing::ds {
    class AdjacencyRowTest : public Test {
    public:
        AdjacencyRowTest() = default;
    };

    TEST_F(AdjacencyRowTest, CheckStats) {
        AdjacencyRow row(193);

        ASSERT_EQ(193, row.num_nodes());
        ASSERT_EQ(4, row.num_blocks());
    }

    TEST_F(AdjacencyRowTest, BitModification) {
        AdjacencyRow row(193);

        row.set();

        ASSERT_EQ(true, row[0]);
        ASSERT_EQ(true, row[45]);
        ASSERT_EQ(true, row[46]);
        ASSERT_EQ(true, row[63]);
        ASSERT_EQ(true, row[64]);
        ASSERT_EQ(true, row[192]);

        row.reset(46);
        for (NodeID u = 53; u < 73; u++) {
            row.reset(u);
        }

        ASSERT_EQ(true, row[0]);
        ASSERT_EQ(true, row[45]);
        ASSERT_EQ(false, row[46]);
        ASSERT_EQ(false, row[63]);
        ASSERT_EQ(false, row[64]);
        ASSERT_EQ(true, row[192]);

        row.flip();
        row.flip(63);

        ASSERT_EQ(false, row[0]);
        ASSERT_EQ(false, row[45]);
        ASSERT_EQ(true, row[46]);
        ASSERT_EQ(false, row[63]);
        ASSERT_EQ(true, row[64]);
        ASSERT_EQ(false, row[192]);

        row.reset();
        row.set(45);

        ASSERT_EQ(false, row[0]);
        ASSERT_EQ(true, row[45]);
        ASSERT_EQ(false, row[46]);
        ASSERT_EQ(false, row[63]);
        ASSERT_EQ(false, row[64]);
        ASSERT_EQ(false, row[192]);
    }

    TEST_F(AdjacencyRowTest, CheckNodeIterator) {
        AdjacencyRow row(193);

        row.set(0);
        for (NodeID u = 10; u < 14; ++u) {
            row.set(u);
        }
        row.set(63);
        row.set(64);
        row.set(127);
        row.set(128);
        row.set(192);
        std::vector<NodeID> nodes{row.begin(), row.end()};

        ASSERT_EQ(nodes.size(), 10);
        ASSERT_EQ(nodes[0], 0);
        ASSERT_EQ(nodes[1], 10);
        ASSERT_EQ(nodes[2], 11);
        ASSERT_EQ(nodes[3], 12);
        ASSERT_EQ(nodes[4], 13);
        ASSERT_EQ(nodes[5], 63);
        ASSERT_EQ(nodes[6], 64);
        ASSERT_EQ(nodes[7], 127);
        ASSERT_EQ(nodes[8], 128);
        ASSERT_EQ(nodes[9], 192);
    }

    TEST_F(AdjacencyRowTest, Move) {
        ASSERT_EQ(std::is_nothrow_move_assignable_v<AdjacencyRow>, true);
        ASSERT_EQ(std::is_nothrow_move_constructible_v<AdjacencyRow>, true);

        AdjacencyRow a(100), b(100);

        a = std::move(b);

        ASSERT_EQ(a.num_nodes(), 100);
        ASSERT_EQ(b.num_nodes(), 0);
    }

    TEST_F(AdjacencyRowTest, Swap) {
        AdjacencyRow a(100), b(100);
        a[10] = true;
        b[20] = true;

        a.swap(b);

        ASSERT_EQ(a[10], false);
        ASSERT_EQ(a[20], true);
        ASSERT_EQ(b[10], true);
        ASSERT_EQ(b[20], false);
    }

    TEST_F(AdjacencyRowTest, ViolateInvariants) {
        auto f = []() {
            AdjacencyRow row(100);
            row[60] = true;
            row.block(1) |= row.block(0);
        };
        ASSERT_DEBUG_DEATH(f(), "");
    }

    TEST_F(AdjacencyRowTest, IndexOutOfBounds) {
        AdjacencyRow row(100);

        ASSERT_EQ(row[0], false);
        ASSERT_EQ(row[50], false);
        ASSERT_DEBUG_DEATH(row[100] = false, "");
        ASSERT_DEBUG_DEATH(row[150] = false, "");
        ASSERT_DEBUG_DEATH(row.set(150), "");
        ASSERT_DEBUG_DEATH(row.reset(150), "");
        ASSERT_DEBUG_DEATH(row.flip(150), "");
        ASSERT_DEBUG_DEATH([[maybe_unused]] auto _ = row.test(150), "");

        ASSERT_EQ(row.block(0), 0);
        ASSERT_EQ(row.block(1), 0);
        ASSERT_DEBUG_DEATH([[maybe_unused]] auto _ = row.block(2), "");

        ASSERT_DEBUG_DEATH([[maybe_unused]] auto _ = row.node_range(50, 10), "");
        ASSERT_DEBUG_DEATH([[maybe_unused]] auto _ = row.node_range(50, 110), "");
    }

    TEST_F(AdjacencyRowTest, DifferentSizes) {
        AdjacencyRow a(50), b(100);

        ASSERT_DEBUG_DEATH(a &= b, "");
        ASSERT_DEBUG_DEATH(a |= b, "");
        ASSERT_DEBUG_DEATH(a ^= b, "");
        ASSERT_DEBUG_DEATH(a -= b, "");
    }

    TEST_F(AdjacencyRowTest, BitRangeOperations) {
        AdjacencyRow x(256), y(256), result(256);
        NodeID a = 0, b = 64, c = 128, d = 192;

        x[a] = true;
        x[b] = true;
        y[a] = true;
        y[c] = true;

        result = x;
        result &= y;
        ASSERT_EQ(result[a], true);
        ASSERT_EQ(result[b], false);
        ASSERT_EQ(result[c], false);
        ASSERT_EQ(result[d], false);

        result = x;
        result |= y;
        ASSERT_EQ(result[a], true);
        ASSERT_EQ(result[b], true);
        ASSERT_EQ(result[c], true);
        ASSERT_EQ(result[d], false);

        result = x;
        result ^= y;
        ASSERT_EQ(result[a], false);
        ASSERT_EQ(result[b], true);
        ASSERT_EQ(result[c], true);
        ASSERT_EQ(result[d], false);

        result = x;
        result -= y;
        ASSERT_EQ(result[a], false);
        ASSERT_EQ(result[b], true);
        ASSERT_EQ(result[c], false);
        ASSERT_EQ(result[d], false);
    }

    TEST_F(AdjacencyRowTest, BitReference) {
        AdjacencyRow row(100);

        row[20] = true;
        row[21] = false;
        ASSERT_EQ(row[20], true);
        ASSERT_EQ(row[21], false);

        row[20] = true;
        row[21] = false;
        row[22] = row[20];
        row[23] = row[21];
        ASSERT_EQ(row[20], row[22]);
        ASSERT_EQ(row[21], row[23]);

        row[20] = true;
        row[21] = false;
        ASSERT_EQ(~row[20], false);
        ASSERT_EQ(~row[21], true);

        row[20] = true;
        row[21] = false;
        row[20].flip();
        row[21].flip();
        ASSERT_EQ(row[20], false);
        ASSERT_EQ(row[21], true);

        row[20] = true;
        row[21] = true;
        row[22] = false;
        row[23] = false;
        row[20] &= true;
        row[21] &= false;
        row[22] &= true;
        row[23] &= false;
        ASSERT_EQ(row[20], true);
        ASSERT_EQ(row[21], false);
        ASSERT_EQ(row[22], false);
        ASSERT_EQ(row[23], false);

        row[20] = true;
        row[21] = true;
        row[22] = false;
        row[23] = false;
        row[20] |= true;
        row[21] |= false;
        row[22] |= true;
        row[23] |= false;
        ASSERT_EQ(row[20], true);
        ASSERT_EQ(row[21], true);
        ASSERT_EQ(row[22], true);
        ASSERT_EQ(row[23], false);

        row[20] = true;
        row[21] = true;
        row[22] = false;
        row[23] = false;
        row[20] ^= true;
        row[21] ^= false;
        row[22] ^= true;
        row[23] ^= false;
        ASSERT_EQ(row[20], false);
        ASSERT_EQ(row[21], true);
        ASSERT_EQ(row[22], true);
        ASSERT_EQ(row[23], false);

        row[20] = true;
        row[21] = true;
        row[22] = false;
        row[23] = false;
        row[20] -= true;
        row[21] -= false;
        row[22] -= true;
        row[23] -= false;
        ASSERT_EQ(row[20], false);
        ASSERT_EQ(row[21], true);
        ASSERT_EQ(row[22], false);
        ASSERT_EQ(row[23], false);
    }

    TEST_F(AdjacencyRowTest, EmptyNodeIterator) {
        AdjacencyRow row(70);

        ASSERT_EQ(row.begin(), row.end());
    }

    TEST_F(AdjacencyRowTest, NodeIterator) {
        AdjacencyRow row(70);
        row[20] = true;
        row[26] = true;

        auto it = row.begin();
        auto it2 = it++;
        ++it;
        ASSERT_EQ(it, row.end());
        ASSERT_DEBUG_DEATH(++it, "");
        ASSERT_EQ(it2, row.begin());
    }


    TEST_F(AdjacencyRowTest, NodeRange) {
        AdjacencyRow row(70);

        std::vector<NodeID> nodes{12, 14, 20, 63, 64, 69};

        for (NodeID u : nodes) {
            row[u] = true;
        }

        auto range1 = row.node_range(0, row.num_nodes());
        ASSERT_EQ(row.begin(), range1.begin());
        ASSERT_EQ(row.end(), range1.end());
        ASSERT_EQ(std::vector<NodeID>(range1.begin(), range1.end()), nodes);

        auto range2 = row.node_range(40, row.num_nodes());
        ASSERT_EQ(std::vector<NodeID>(range2.begin(), range2.end()), std::vector<NodeID>({63, 64, 69}));

        auto range3 = row.node_range(30, 50);
        ASSERT_EQ(range3.begin(), range3.end());

        auto range4 = row.node_range(30, 63);
        ASSERT_EQ(range4.begin(), range4.end());
    }

    TEST_F(AdjacencyRowTest, ForEach) {
        AdjacencyRow row(70);

        std::vector<NodeID> nodes{12, 14, 20, 63, 64, 69};

        for (NodeID u : nodes) {
            row[u] = true;
        }

        auto it = row.begin();
        row.for_each([&](NodeID u) {
            ASSERT_EQ(u, *it++);
        });

        std::vector<NodeID> nodes2;
        row.for_each([&](NodeID u) {
            nodes2.push_back(u);
        });

        ASSERT_EQ(nodes, nodes2);
    }

    TEST_F(AdjacencyRowTest, TriangleIterationWithForbiddenNodePairs) {
        NodeID n = 16;
        // neighbors u: 0..7
        // neighbors v: 0..3, 8..11
        AdjacencyRow u_row(n), v_row(n);

        for (NodeID w = 0; w < n / 2; ++w) {
            u_row[w] = true;
        }

        for (NodeID w = 0; w < n / 4; ++w) {
            v_row[w] = true;
        }

        for (NodeID w = n / 2; w < 3 * n / 4; ++w) {
            v_row[w] = true;
        }

        // forbidden node pairs u: 0, 4, 8, 12
        // forbidden node pairs v: 2, 6, 10, 14
        AdjacencyRow u_forbidden_row(n), v_forbidden_row(n);
        for (NodeID w = 0; w < n; ++w) {
            if (w % 4 == 0)
                u_forbidden_row[w] = true;
        }
        for (NodeID w = 0; w < n; ++w) {
            if (w % 4 == 2)
                v_forbidden_row[w] = true;
        }


        AdjacencyRow w_row = u_row;
        for (unsigned int i = 0; i < w_row.num_blocks(); ++i) {
            w_row.block(i) &= v_row.block(i) & ~u_forbidden_row.block(i) & ~v_forbidden_row.block(i);
        }

        std::vector<NodeID> nodes;
        w_row.for_each([&](NodeID u) {
            nodes.push_back(u);
        });

        ASSERT_EQ(nodes.size(), 2);
        ASSERT_EQ(nodes[0], 1);
        ASSERT_EQ(nodes[1], 3);
    }

    TEST_F(AdjacencyRowTest, Benchmark) {
        NodeID n = 300;
        std::uniform_int_distribution<NodeID> dist(0, n - 1);
        std::mt19937_64 gen(0);

        int num_iter = 10000;

        std::cout << "n\tdensity\tnum_iter\tNodeIterator [s]\tfor_each [s]\tindex operator [s]\n";

        for (auto alpha : {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {
            AdjacencyRow row(n);

            for (NodeID i = 0; i < n * alpha; ++i) {
                row[dist(gen)] = true;
            }

            std::cout << n << "\t" << alpha << "\t" << num_iter << "\t";

            {
                NodeID sum = 0;
                auto start = std::chrono::steady_clock::now();
                for (int i = 0; i < num_iter; ++i) {
                    for (auto u : row) {
                        sum += u;
                    }
                }
                std::cout << std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::steady_clock::now() - start).count() * 1e-6 << "\t";
            }

            {
                NodeID sum = 0;
                auto start = std::chrono::steady_clock::now();
                for (int i = 0; i < num_iter; ++i) {
                    row.for_each([&](auto u) {
                        sum += u;
                    });
                }
                std::cout << std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::steady_clock::now() - start).count() * 1e-6 << "\t";
            }

            {
                NodeID sum = 0;
                auto start = std::chrono::steady_clock::now();
                for (int i = 0; i < num_iter; ++i) {
                    for (NodeID u = 0; u < n; ++u) {
                        if (row[u]) {
                            sum += u;
                        }
                    }
                }
                std::cout << std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::steady_clock::now() - start).count() * 1e-6 << "\n";
            }
        }
    }
} // namespace cluster_editing::ds
