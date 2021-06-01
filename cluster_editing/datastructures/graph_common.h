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
} // namespace cluster_editing