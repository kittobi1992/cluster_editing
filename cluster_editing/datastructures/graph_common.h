#pragma once

#include <cstdint>
#include <limits>

namespace cluster_editing {

  using NodeID = uint32_t;
  using CliqueID = uint32_t;
  using EdgeID = uint32_t;
  using NodeWeight = int32_t;
  using EdgeWeight = int32_t;

  // Constant Declarations
  static constexpr NodeID INVALID_NODE = std::numeric_limits<NodeID>::max();
  static constexpr CliqueID INVALID_CLIQUE = std::numeric_limits<CliqueID>::max();
  static constexpr EdgeID INVALID_EDGE = std::numeric_limits<EdgeID>::max();

  struct Neighbor {
    EdgeID id;
    NodeID target;
  };

} // namespace cluster_editing