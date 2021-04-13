#include "graph.h"

#include <vector>

#include "cluster_editing/macros.h"
#include "cluster_editing/datastructures/sparse_map.h"
#include "cluster_editing/utils/timer.h"

namespace cluster_editing {
namespace ds {

  void Graph::reset() {
    for ( const NodeID& u : nodes() ) {
      _nodes[u].setClique(u);
    }
  }
}
} // namespace cluster_editing::ds