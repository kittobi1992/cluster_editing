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

#include "cluster_editing/definitions.h"

namespace cluster_editing {

class IRefiner {

 public:
  IRefiner(const IRefiner&) = delete;
  IRefiner(IRefiner&&) = delete;
  IRefiner & operator= (const IRefiner &) = delete;
  IRefiner & operator= (IRefiner &&) = delete;

  virtual ~IRefiner() = default;

  void initialize(Graph& graph) {
    initializeImpl(graph);
  }

  EdgeWeight refine(Graph& graph,
                    const EdgeWeight current_edits,
                    const EdgeWeight target_edits = 0) {
    return refineImpl(graph, current_edits, target_edits);
  }

 protected:
  IRefiner() = default;

 private:
  virtual void initializeImpl(Graph& graph) = 0;
  virtual EdgeWeight refineImpl(Graph& graph,
                                const EdgeWeight current_edits,
                                const EdgeWeight target_edits) = 0;
};

} // namespace cluster_editing