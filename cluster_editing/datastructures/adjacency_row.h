
#ifndef CLUSTEREDITING_ADJACENCY_ROW_H
#define CLUSTEREDITING_ADJACENCY_ROW_H

#include <vector>
#include <cassert>

#include "gmock/gmock.h"

#include "graph_common.h"

namespace cluster_editing {
namespace ds {
    namespace detail {
        [[nodiscard]] constexpr int ctz(unsigned long long n) noexcept {
            assert(n != 0);
            return __builtin_ctzll(n);
        }
    }

    class AdjacencyRow {
    public:
        using block_type = unsigned long long;
        using block_width_type = unsigned int;
        static constexpr block_width_type bits_per_block = std::numeric_limits<block_type>::digits;
        using block_index_type = unsigned int;  // typename std::vector<block_type>::size_type;

    private:
        NodeID m_num_nodes{};
        std::vector<block_type> m_blocks;

        class reference {
        private:
            block_type &m_block;
            const block_type m_mask;

        public:
            void operator&() = delete; //NOLINT

            reference(block_type &b, block_width_type bit_pos) noexcept:
                    m_block(b), m_mask(block_type(1) << bit_pos) {
                assert(bit_pos < bits_per_block);
            }

            [[nodiscard]] operator bool() const noexcept { // NOLINT(google-explicit-constructor)
                return (m_block & m_mask) != 0;
            }

            [[nodiscard]] bool operator~() const noexcept {
                return (m_block & m_mask) == 0;
            }

            reference &flip() noexcept {
                m_block ^= m_mask;
                return *this;
            }

            reference &operator=(bool x) noexcept {
                if (x) {
                    m_block |= m_mask;
                } else {
                    m_block &= ~m_mask;
                }
                return *this;
            }

            reference &operator=(const reference &rhs) noexcept {
                *this = bool(rhs);
                return *this;
            }

            reference &operator|=(bool x) noexcept {
                if (x) {
                    m_block |= m_mask;
                }
                return *this;
            }

            reference &operator&=(bool x) noexcept {
                if (!x) {
                    m_block &= ~m_mask;
                }
                return *this;
            }

            reference &operator^=(bool x) noexcept {
                if (x) {
                    m_block ^= m_mask;
                }
                return *this;
            }

            reference &operator-=(bool x) noexcept {
                if (x) {
                    m_block &= ~m_mask;
                }
                return *this;
            }
        };

        using const_reference = bool;

        class NodeIterator {
        private:
            const block_type *m_blocks{nullptr};
            NodeID m_node{};
            NodeID m_end_node{};

        public:
            using value_type = NodeID;
            using difference_type = std::ptrdiff_t;
            using pointer = void;
            using reference = value_type;
            using iterator_category = std::forward_iterator_tag;

            constexpr NodeIterator() noexcept = default;

            NodeIterator(const block_type *blocks, NodeID start_node, NodeID end_node) :
                    m_blocks(blocks), m_node(start_node), m_end_node(end_node) {
                assert(m_node <= m_end_node);
                if (m_node < m_end_node && (m_blocks[block_index(m_node)] & bit_mask(m_node)) == block_type(0)) {
                    ++*this;
                }
            }

            constexpr NodeIterator(const NodeIterator &other) noexcept = default;

            constexpr NodeIterator &operator=(const NodeIterator &other) noexcept = default;

            [[nodiscard]] constexpr value_type operator*() const {
                return m_node;
            }

            NodeIterator &operator++() noexcept;

            [[nodiscard]] NodeIterator operator++(int) noexcept { // NOLINT(cert-dcl21-cpp)
                NodeIterator copy(*this);
                ++*this;
                return copy;
            }

            [[nodiscard]] constexpr bool operator==(const NodeIterator &other) const noexcept {
                assert(m_blocks == other.m_blocks);
                assert(m_end_node == other.m_end_node);
                return m_node == other.m_node;
            }

            [[nodiscard]] constexpr bool operator!=(const NodeIterator &other) const noexcept {
                assert(m_blocks == other.m_blocks);
                assert(m_end_node == other.m_end_node);
                return m_node != other.m_node;
            }
        };

        class NodeRange {
            NodeIterator m_begin;
            NodeIterator m_end;
        public:
            NodeRange(const AdjacencyRow &row, NodeID start_node, NodeID end_node) :
                    m_begin(row.m_blocks.data(), start_node, end_node),
                    m_end(row.m_blocks.data(), end_node, end_node) {}

            [[nodiscard]] constexpr NodeIterator begin() const noexcept { return m_begin; }

            [[nodiscard]] constexpr NodeIterator end() const noexcept { return m_end; }
        };

    public:

        AdjacencyRow() noexcept = default;

        explicit AdjacencyRow(NodeID num_nodes) : m_num_nodes(num_nodes), m_blocks(calc_num_blocks(m_num_nodes)) {}

        AdjacencyRow(const AdjacencyRow &b) = default;

        ~AdjacencyRow() noexcept {
            assert(check_invariants());
        }

        void swap(AdjacencyRow &other) noexcept {
            using std::swap;
            swap(m_blocks, other.m_blocks);
            swap(m_num_nodes, other.m_num_nodes);
        }

        AdjacencyRow &operator=(const AdjacencyRow &b) = default;

        AdjacencyRow(AdjacencyRow &&src) noexcept: m_num_nodes(src.m_num_nodes), m_blocks(std::move(src.m_blocks)) {
            src.m_num_nodes = 0;
        }

        AdjacencyRow &operator=(AdjacencyRow &&src) noexcept {
            m_num_nodes = src.m_num_nodes;
            m_blocks = std::move(src.m_blocks);
            src.m_num_nodes = 0;
            return *this;
        }

        bool operator==(const AdjacencyRow &rhs) const {
            return m_num_nodes == rhs.m_num_nodes && m_blocks == rhs.m_blocks;
        }

        bool operator!=(const AdjacencyRow &rhs) const {
            return !(*this == rhs);
        }

        AdjacencyRow &operator&=(const AdjacencyRow &rhs) {
            assert(num_nodes() == rhs.num_nodes());
            for (block_type i = 0; i < num_blocks(); ++i)
                m_blocks[i] &= rhs.m_blocks[i];
            return *this;
        }

        AdjacencyRow &operator|=(const AdjacencyRow &rhs) {
            assert(num_nodes() == rhs.num_nodes());
            for (block_index_type i = 0; i < num_blocks(); ++i)
                m_blocks[i] |= rhs.m_blocks[i];
            return *this;
        }

        AdjacencyRow &operator^=(const AdjacencyRow &rhs) {
            assert(num_nodes() == rhs.num_nodes());
            for (block_index_type i = 0; i < num_blocks(); ++i)
                m_blocks[i] ^= rhs.m_blocks[i];
            return *this;

        }

        AdjacencyRow &operator-=(const AdjacencyRow &rhs) {
            assert(num_nodes() == rhs.num_nodes());
            for (block_index_type i = 0; i < num_blocks(); ++i)
                m_blocks[i] &= ~rhs.m_blocks[i];
            return *this;
        }

        AdjacencyRow &set() noexcept {
            std::fill(m_blocks.begin(), m_blocks.end(), ~static_cast<block_type>(0));
            zero_unused_bits();
            return *this;
        }

        AdjacencyRow &set(NodeID node) {
            assert(node < m_num_nodes);
            m_blocks[block_index(node)] |= bit_mask(node);
            return *this;
        }

        AdjacencyRow &reset() noexcept {
            std::fill(m_blocks.begin(), m_blocks.end(), block_type(0));
            return *this;
        }

        AdjacencyRow &reset(NodeID node) {
            assert(node < m_num_nodes);
            m_blocks[block_index(node)] &= ~bit_mask(node);
            return *this;
        }

        AdjacencyRow &flip() noexcept {
            for (block_index_type i = 0; i < num_blocks(); ++i)
                m_blocks[i] = ~m_blocks[i];
            zero_unused_bits();
            return *this;
        }

        AdjacencyRow &flip(NodeID node) {
            assert(node < m_num_nodes);
            m_blocks[block_index(node)] ^= bit_mask(node);
            return *this;
        }

        [[nodiscard]] bool test(NodeID node) const {
            assert(node < m_num_nodes);
            return (m_blocks[block_index(node)] & bit_mask(node)) != 0;
        }

        [[nodiscard]] reference operator[](NodeID node) {
            assert(node < m_num_nodes);
            return reference(m_blocks[block_index(node)], bit_index(node));
        }

        [[nodiscard]] const_reference operator[](NodeID node) const { return test(node); }

    private:
        // Disabled NodeIterator for public use. See AdjacencyRowTest.Benchmark.

        FRIEND_TEST(AdjacencyRowTest, CheckNodeIterator);

        FRIEND_TEST(AdjacencyRowTest, IndexOutOfBounds);

        FRIEND_TEST(AdjacencyRowTest, EmptyNodeIterator);

        FRIEND_TEST(AdjacencyRowTest, NodeIterator);

        FRIEND_TEST(AdjacencyRowTest, NodeRange);

        FRIEND_TEST(AdjacencyRowTest, ForEach);

        FRIEND_TEST(AdjacencyRowTest, Benchmark);

        [[nodiscard]] NodeIterator begin() const noexcept { return NodeIterator(m_blocks.data(), 0, m_num_nodes); }

        [[nodiscard]] NodeIterator end() const noexcept {
            return NodeIterator(m_blocks.data(), m_num_nodes, m_num_nodes);
        }

        [[nodiscard]] NodeRange node_range(NodeID start_node, NodeID end_node) const {
            assert(start_node <= end_node);
            assert(end_node <= m_num_nodes);
            return NodeRange(*this, start_node, end_node);
        }

    public:

        [[nodiscard]] constexpr NodeID num_nodes() const noexcept { return m_num_nodes; }

        [[nodiscard]] block_index_type num_blocks() const noexcept { return m_blocks.size(); }

        [[nodiscard]] block_type &block(block_index_type index) {
            assert(index < num_blocks());
            return m_blocks[index];
        }

        [[nodiscard]] const block_type &block(block_index_type index) const {
            assert(index < num_blocks());
            return m_blocks[index];
        }

        void zero_unused_bits() noexcept;

        /**
         * Make sure that \p f does mutate the AdjacencyRows where \a block(i) is called on.
         *
         * Additionally \p get_block must be a subset of at least one \a block(i) from an AdjacencyRow with the
         * desired size. AdjacencyRow holds up the invariant that the bits in the last block beyond the last valid bit
         * are zeroed out. Disregarding this might lead \p f being called with nodes beyond the last valid bit.
         *
         * @param num_blocks
         * @param get_block
         * @param f
         * @return
         */
        template<class B, class F>
        constexpr static void for_each(block_index_type num_blocks, B get_block, F f) {
            static_assert(std::is_invocable_r_v<block_type, B, block_index_type>);
            static_assert(std::is_invocable_r_v<void, F, NodeID>);

            for (block_index_type i = 0; i < num_blocks; i++) {
                auto block = get_block(i);
                while (block) {
                    const block_width_type bit_pos = detail::ctz(block);
                    block ^= (block_type(1) << bit_pos);
                    const NodeID v = i * bits_per_block + bit_pos;
                    f(v);
                }
            }
        }

        template<class F>
        constexpr void for_each(F f) const {
            return for_each(num_blocks(), [&](auto i) { return m_blocks[i]; }, f);
        }

    private:

        [[nodiscard]] bool check_invariants() const noexcept;

        [[nodiscard]] constexpr block_width_type count_extra_bits() const noexcept { return bit_index(num_nodes()); }

        [[nodiscard]] constexpr static block_type block_index(NodeID node) noexcept { return node / bits_per_block; }

        [[nodiscard]] constexpr static block_width_type bit_index(NodeID node) noexcept {
            return static_cast<block_width_type>(node % bits_per_block);
        }

        [[nodiscard]] constexpr static block_type bit_mask(NodeID node) noexcept {
            return block_type(1) << bit_index(node);
        }

        [[nodiscard]] constexpr static block_index_type calc_num_blocks(NodeID num_nodes) noexcept {
            return num_nodes / bits_per_block
                   + static_cast<block_index_type>(num_nodes % bits_per_block != 0);
        }

        [[nodiscard]] block_type &highest_block() {
            assert(num_nodes() > 0 && num_blocks() > 0);
            return m_blocks.back();
        }

        [[nodiscard]] const block_type &highest_block() const {
            assert(num_nodes() > 0 && num_blocks() > 0);
            return m_blocks.back();
        }
    };
}
}


#endif //CLUSTEREDITING_ADJACENCY_ROW_H
