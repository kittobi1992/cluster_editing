#include "multilevel.h"

#include "cluster_editing/macros.h"
#include "cluster_editing/coarsening/do_nothing_coarsener.h"
#include "cluster_editing/refinement/do_nothing_refiner.h"

namespace cluster_editing {
namespace multilevel {

std::unique_ptr<ICoarsener> instantiateCoarsener(Graph& graph, const Context& context) {
  switch(context.coarsening.algorithm) {
    case CoarseningAlgorithm::do_nothing: return std::make_unique<DoNothingCoarsener>(graph, context);
    case CoarseningAlgorithm::UNDEFINED:
      ERROR("No valid coarsening algorithm");
      return nullptr;
  }
  return nullptr;
}

std::unique_ptr<IRefiner> instantiateRefiner(const Context& context) {
  switch(context.refinement.algorithm) {
    case RefinementAlgorithm::do_nothing: return std::make_unique<DoNothingRefiner>();
    case RefinementAlgorithm::UNDEFINED:
      ERROR("No valid refinement algorithm");
      return nullptr;
  }
  return nullptr;
}

void solve(Graph& graph, const Context& context) {
  std::unique_ptr<ICoarsener> coarsener = instantiateCoarsener(graph, context);
  std::unique_ptr<IRefiner> refiner = instantiateRefiner(context);

  coarsener->coarsen();

  // Do some other stuff before we start uncoarsening

  coarsener->uncoarsen(refiner);
}

} // namespace multilevel
} // namespace cluster_editing