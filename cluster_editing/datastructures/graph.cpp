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

#include "graph.h"

#include <vector>

#include "cluster_editing/macros.h"
#include "cluster_editing/datastructures/sparse_map.h"
#include "cluster_editing/utils/timer.h"

namespace cluster_editing::ds {

  void Graph::reset() {
    for ( const NodeID& u : nodes() ) {
      _nodes[u].setClique(u);
    }
  }

} // namespace cluster_editing::ds