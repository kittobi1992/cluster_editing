#pragma once

#include "cluster_editing/refinement/i_refiner.h"

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
