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

#pragma once

#include <vector>
#include <algorithm>

#include "cluster_editing/macros.h"
#include "cluster_editing/datastructures/graph_common.h"
#include "cluster_editing/datastructures/spin_lock.h"
#include "cluster_editing/utils/range.h"
#include "cluster_editing/utils/randomize.h"

namespace cluster_editing {
namespace ds {

// Forward
class GraphFactory;

class Graph {

 private:
  class Node {

   public:
    Node() :
      _begin(0),
      _clique(INVALID_CLIQUE) { }

    size_t firstEntry() const {
      return _begin;
    }

    void setFirstEntry(const size_t begin) {
      _begin = begin;
    }

    CliqueID clique() const {
      return _clique;
    }

    void setClique(const CliqueID clique) {
      _clique = clique;
    }

   private:
    size_t _begin;
    CliqueID _clique;
  };

  /*!
   * Iterator for GraphElements (Nodes/Edges)
   */
  template <typename IDType>
  class GraphElementIterator :
    public std::iterator<std::forward_iterator_tag, // iterator_category
                         IDType,                    // value_type
                         std::ptrdiff_t,            // difference_type
                         const IDType*,             // pointer
                         IDType> {                  // reference
   public:
    /*!
     * Construct a GraphElementIterator
     * See Graph::nodes() or Graph::edges() for usage.
     *
     * \param id The index of the element the pointer points to
     * \param max_id The maximum index allowed
     */
    GraphElementIterator(IDType id, IDType max_id) :
      _id(id),
      _max_id(max_id) { }

    // ! Returns the id of the element the iterator currently points to.
    IDType operator* () const {
      return _id;
    }

    // ! Prefix increment. The iterator advances to the next valid element.
    GraphElementIterator & operator++ () {
      if ( _id < _max_id ) {
        ++_id;
      }
      return *this;
    }

    // ! Postfix increment. The iterator advances to the next valid element.
    GraphElementIterator operator++ (int) {
      GraphElementIterator copy = *this;
      operator++ ();
      return copy;
    }

    bool operator!= (const GraphElementIterator& rhs) {
      return _id != rhs._id;
    }

    bool operator== (const GraphElementIterator& rhs) {
      return _id == rhs._id;
    }

   private:
    // Handle to the HypergraphElement the iterator currently points to
    IDType _id = 0;
    // Maximum allowed index
    IDType _max_id = 0;
  };

  // ! Iterator to iterate over the nodes
  using NodeIterator = GraphElementIterator<NodeID>;
  // ! Iterator to iterate over the edges
  using EdgeIterator = GraphElementIterator<EdgeID>;
  // ! Iterator to iterate over incident edges
  using NeighborIterator = typename std::vector<NodeID>::const_iterator;

 public:
  Graph() :
    _num_nodes(0),
    _num_edges(0),
    _total_weight(0),
    _max_degree(0),
    _nodes(),
    _edges(),
    _best_cliques(),
    _best_quality(std::numeric_limits<EdgeWeight>::max()) { }

  // ####################### General Graph Stats #######################

  NodeID numNodes() const {
    return _num_nodes;
  }

  EdgeID numEdges() const {
    return _num_edges;
  }

  NodeWeight totalWeight() const {
    return _total_weight;
  }

  NodeID maxDegree() const {
    return _max_degree;
  }

  // ########################## Iterators ##########################

  // ! Returns a range of all nodes of the graph
  IteratorRange<NodeIterator> nodes() const {
    return IteratorRange<NodeIterator>(
      NodeIterator(0, _num_nodes), NodeIterator(_num_nodes, _num_nodes));
  }

  // ! Returns a range of all edges of the graph
  IteratorRange<EdgeIterator> edges() const {
    return IteratorRange<EdgeIterator>(
      EdgeIterator(0, _num_edges), EdgeIterator(_num_edges, _num_edges));
  }

  // ! Returns a range to iterate over all neighbors of a vertex
  IteratorRange<NeighborIterator> neighbors(const NodeID u) const {
    return IteratorRange<NeighborIterator>(
      _edges.cbegin() + node(u).firstEntry(),
      _edges.cbegin() + node(u + 1).firstEntry());
  }

  void sortNeighborsByCliqueID(const NodeID u) {
    ASSERT(u <= _num_nodes, "Node" << u << "does not exist");
    std::sort(_edges.begin() + node(u).firstEntry(),
              _edges.begin() + node(u + 1).firstEntry(),
              [&](const NodeID& lhs, const NodeID& rhs) {
                return clique(lhs) < clique(rhs);
              });
  }

  NodeID randomNeighbor(const NodeID u) const {
    if ( degree(u) > 0 ) {
      const int start_idx = node(u).firstEntry();
      const int end_idx = node(u + 1).firstEntry();
      const int random_tie_breaking =
        utils::Randomize::instance().getRandomInt(start_idx, end_idx - 1);
      return _edges[random_tie_breaking];
    } else {
      return INVALID_NODE;
    }
  }

  // ####################### Node Information #######################

  size_t degree(const NodeID u) const {
    return node(u + 1).firstEntry() - node(u).firstEntry();
  }

  CliqueID clique(const NodeID u) const {
    return node(u).clique();
  }

  CliqueID bestClique(const NodeID u) const {
    return _best_cliques[u];
  }

  void setClique(const NodeID u, const CliqueID id) {
    node(u).setClique(id);
  }

  // ####################### Checkpointing #######################

  EdgeWeight bestEdits() const {
    return _best_quality;
  }

  void applyBestCliques() {
    copy_lock.lock();
    for ( const NodeID& u : nodes() ) {
      setClique(u, _best_cliques[u]);
    }
    copy_lock.unlock();
  }

  void checkpoint(const EdgeWeight quality) {
    copy_lock.lock();
    if ( quality < _best_quality ) {
      for ( const NodeID& u : nodes() ) {
        _best_cliques[u] = clique(u);
      }
      _best_quality = quality;
    }
    copy_lock.unlock();
  }

  Graph copyBestSolution() {
    copy_lock.lock();
    Graph cpy(*this);
    cpy.applyBestCliques();
    copy_lock.unlock();
    return cpy;
  }

  // ####################### Reset #######################
  void reset();

 private:
  friend class GraphFactory;

  // ####################### Node Information #######################

  // ! Accessor for node-related information
  ATTRIBUTE_ALWAYS_INLINE const Node& node(const NodeID u) const {
    ASSERT(u <= _num_nodes, "Node" << u << "does not exist");
    return _nodes[u];
  }

  // ! To avoid code duplication we implement non-const version in terms of const version
  ATTRIBUTE_ALWAYS_INLINE Node& node(const NodeID u) {
    return const_cast<Node&>(static_cast<const Graph&>(*this).node(u));
  }

  NodeID _num_nodes;
  EdgeID _num_edges;
  NodeWeight _total_weight;
  NodeID _max_degree;

  // ! Index Array
  std::vector<Node> _nodes;
  // ! Indicence Array
  std::vector<NodeID> _edges;

  // Checkpointing
  std::vector<CliqueID> _best_cliques;
  EdgeWeight _best_quality;
  SpinLock copy_lock;
};

} // namespace ds
} // namespace cluster_editing