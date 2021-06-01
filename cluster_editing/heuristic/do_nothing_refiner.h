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

#include "cluster_editing/heuristic/i_refiner.h"

namespace cluster_editing {
class DoNothingRefiner final : public IRefiner {
 public:
  template <typename ... Args>
  explicit DoNothingRefiner(Args&& ...) noexcept { }
  DoNothingRefiner(const DoNothingRefiner&) = delete;
  DoNothingRefiner(DoNothingRefiner&&) = delete;
  DoNothingRefiner & operator= (const DoNothingRefiner &) = delete;
  DoNothingRefiner & operator= (DoNothingRefiner &&) = delete;

 private:
  void initializeImpl(Graph&) override final { }

  EdgeWeight refineImpl(Graph&,
                        const EdgeWeight current_edits
                        const EdgeWeight) override final { return current_edits; }
};
}  // namespace cluster_editing
