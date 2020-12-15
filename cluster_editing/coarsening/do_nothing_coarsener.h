#pragma once

#include "cluster_editing/coarsening/i_coarsener.h"

namespace cluster_editing {
class DoNothingCoarsener final : public ICoarsener {
 public:
  template <typename ... Args>
  explicit DoNothingCoarsener(Args&& ...) noexcept { }
  DoNothingCoarsener(const DoNothingCoarsener&) = delete;
  DoNothingCoarsener(DoNothingCoarsener&&) = delete;
  DoNothingCoarsener & operator= (const DoNothingCoarsener &) = delete;
  DoNothingCoarsener & operator= (DoNothingCoarsener &&) = delete;

 private:
  void coarsenImpl() override final {
    // do nothing
  }

  void uncoarsenImpl(std::unique_ptr<IRefiner>&) override final {
    // Call refiner:
    // refiner->initialize();
    // refiner->refine(graph);
  }
};
}  // namespace cluster_editing
