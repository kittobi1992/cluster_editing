
#include "adjacency_row.h"

namespace cluster_editing::ds {

    namespace detail {
        [[nodiscard]] constexpr int ctz(unsigned long long n) noexcept {
            assert(n != 0);
            return __builtin_ctzll(n);
        }
    }

    AdjacencyRow::NodeIterator &AdjacencyRow::NodeIterator::operator++() noexcept {
        assert(m_node < m_end_node);

        ++m_node;

        if (m_node >= m_end_node) {
            m_node = m_end_node;
            return *this;
        }

        auto block_idx = block_index(m_node);
        const auto idx = bit_index(m_node);

        assert(block_idx < calc_num_blocks(m_end_node));
        const auto remaining_block = m_blocks[block_idx] >> idx;
        if (remaining_block != block_type(0)) {
            m_node += static_cast<NodeID>(detail::ctz(remaining_block));
            if (m_node > m_end_node)
                m_node = m_end_node;
            return *this;
        }

        const auto num_blocks = calc_num_blocks(m_end_node);

        ++block_idx;
        while (block_idx < num_blocks && (m_blocks[block_idx] == block_type(0)))
            ++block_idx;

        m_node = block_idx < num_blocks
                 ? block_idx * bits_per_block + static_cast<NodeID>(detail::ctz(m_blocks[block_idx]))
                 : m_end_node;

        assert(m_node <= m_end_node);
        return *this;
    }

    void AdjacencyRow::zero_unused_bits() noexcept {
        assert(num_blocks() == calc_num_blocks(m_num_nodes));

        // if != 0 this is the number of bits used in the last block
        const block_width_type extra_bits = count_extra_bits();

        if (extra_bits != 0)
            highest_block() &= (block_type(1) << extra_bits) - 1;
    }

    bool AdjacencyRow::check_invariants() const noexcept {
        const block_width_type extra_bits = count_extra_bits();
        if (extra_bits > 0) {
            const block_type mask = (~block_type(0)) << extra_bits;
            if ((highest_block() & mask) != 0)
                return false;
        }
        if (m_blocks.size() > m_blocks.capacity() || num_blocks() != calc_num_blocks(num_nodes()))
            return false;

        return true;
    }
}