#pragma once

#include <vector>

#include "cluster_editing/macros.h"
#include "cluster_editing/datastructures/graph_common.h"
#include "cluster_editing/datastructures/spin_lock.h"
#include "cluster_editing/utils/range.h"

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

  class Edge {
   public:
    Edge() :
      _source(INVALID_NODE),
      _target(INVALID_NODE) { }

    NodeID source() const {
      return _source;
    }

    void setSource(const NodeID source) {
      _source = source;
    }

    NodeID target() const {
      return _target;
    }

    void setTarget(const NodeID target) {
      _target = target;
    }

   private:
    NodeID _source;
    NodeID _target;
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

  class NeighborIterator :
    public std::iterator<std::forward_iterator_tag, // iterator_category
                         Neighbor,                  // value_type
                         std::ptrdiff_t,            // difference_type
                         const Neighbor*,           // pointer
                         Neighbor> {                // reference
   public:
    NeighborIterator(const Edge* start_edge, EdgeID start_id) :
      _current_edge(start_edge),
      _current_neighbor(Neighbor { start_id, start_edge->target() }) { }

    // ! Returns the id of the element the iterator currently points to.
    Neighbor operator* () const {
      return _current_neighbor;
    }

    // ! Prefix increment.
    NeighborIterator & operator++ () {
      ++_current_neighbor.id;
      ++_current_edge;
      _current_neighbor.target = _current_edge->target();
      return *this;
    }

    // ! Postfix increment.
    NeighborIterator operator++ (int) {
      NeighborIterator copy = *this;
      operator++ ();
      return copy;
    }

    bool operator!= (const NeighborIterator& rhs) {
      return _current_neighbor.id != rhs._current_neighbor.id;
    }

    bool operator== (const NeighborIterator& rhs) {
      return _current_neighbor.id == rhs._current_neighbor.id;
    }

   private:
    const Edge* _current_edge;
    Neighbor _current_neighbor;
  };

  // ! Iterator to iterate over the nodes
  using NodeIterator = GraphElementIterator<NodeID>;
  // ! Iterator to iterate over the edges
  using EdgeIterator = GraphElementIterator<EdgeID>;

 public:
  Graph() :
    _num_nodes(0),
    _num_edges(0),
    _total_weight(0),
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
      NeighborIterator(_edges.data() + node(u).firstEntry(), node(u).firstEntry()),
      NeighborIterator(_edges.data() + node(u + 1).firstEntry(), node(u + 1).firstEntry()));
  }

  // ####################### Node Information #######################

  NodeWeight nodeWeight(const NodeID) const {
    return 1;
  }

  NodeWeight selfloopWeight(const NodeID) const {
    return 0;
  }

  size_t degree(const NodeID u) const {
    return node(u + 1).firstEntry() - node(u).firstEntry();
  }

  EdgeWeight weightedDegree(const NodeID u) const {
    return degree(u);
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

  // ####################### Edge Information #######################

  NodeID source(const EdgeID e) const {
    return edge(e).source();
  }

  NodeID target(const EdgeID e) const {
    return edge(e).target();
  }

  EdgeWeight edgeWeight(const EdgeID) const {
    return 1;
  }

  // ####################### Checkpointing #######################

  void applyBestCliques() {
    copy_lock.lock();
    for ( const NodeID& u : nodes() ) {
      setClique(u, _best_cliques[u]);
    }
    copy_lock.unlock();
  }

  void checkpoint(const EdgeWeight quality) {
    copy_lock.lock();
    if ( quality <= _best_quality ) {
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

  // ####################### Edge Information #######################

  // ! Accessor for edge-related information
  ATTRIBUTE_ALWAYS_INLINE const Edge& edge(const EdgeID e) const {
    ASSERT(e <= _num_edges, "Edge" << e << "does not exist");
    return _edges[e];
  }

  // ! To avoid code duplication we implement non-const version in terms of const version
  ATTRIBUTE_ALWAYS_INLINE Edge& edge(const EdgeID e) {
    return const_cast<Edge&>(static_cast<const Graph&>(*this).edge(e));
  }

  NodeID _num_nodes;
  EdgeID _num_edges;
  NodeWeight _total_weight;

  // ! Index Array
  std::vector<Node> _nodes;
  // ! Indicence Array
  std::vector<Edge> _edges;

  // Checkpointing
  std::vector<CliqueID> _best_cliques;
  EdgeWeight _best_quality;
  SpinLock copy_lock;
};

} // namespace ds
} // namespace cluster_editing