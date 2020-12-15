#pragma once

#include <memory>

#include "cluster_editing/definitions.h"
#include "cluster_editing/refinement/i_refiner.h"

namespace cluster_editing {

class ICoarsener {

 public:
  ICoarsener(const ICoarsener&) = delete;
  ICoarsener(ICoarsener&&) = delete;
  ICoarsener & operator= (const ICoarsener &) = delete;
  ICoarsener & operator= (ICoarsener &&) = delete;

  void coarsen() {
    coarsenImpl();
  }

  void uncoarsen(std::unique_ptr<IRefiner>& refiner) {
    return uncoarsenImpl(refiner);
  }

  virtual ~ICoarsener() = default;

 protected:
  ICoarsener() = default;

 private:
  virtual void coarsenImpl() = 0;
  virtual void uncoarsenImpl(std::unique_ptr<IRefiner>& refiner) = 0;
};

} // namespace cluster_editing