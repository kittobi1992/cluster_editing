#include "multilevel.h"

#include "cluster_editing/macros.h"
#include "cluster_editing/coarsening/do_nothing_coarsener.h"
#include "cluster_editing/coarsening/lp_coarsener.h"
#include "cluster_editing/refinement/do_nothing_refiner.h"
#include "cluster_editing/refinement/lp_refiner.h"
#include "cluster_editing/utils/timer.h"
#include "cluster_editing/io/output.h"

namespace cluster_editing {
namespace multilevel {

std::unique_ptr<ICoarsener> instantiateCoarsener(Graph& graph, const Context& context) {
  switch(context.coarsening.algorithm) {
    case CoarseningAlgorithm::do_nothing: return std::make_unique<DoNothingCoarsener>(graph, context);
    case CoarseningAlgorithm::lp_coarsener: return std::make_unique<LabelPropagationCoarsener>(graph, context);
    case CoarseningAlgorithm::UNDEFINED:
      ERROR("No valid coarsening algorithm");
      return nullptr;
  }
  return nullptr;
}

void solve(Graph& graph, const Context& context) {
  std::unique_ptr<ICoarsener> coarsener = instantiateCoarsener(graph, context);

  io::printCoarseningBanner(context);
  utils::Timer::instance().start_timer("coarsening", "Coarsening");
  coarsener->coarsen();
  utils::Timer::instance().stop_timer("coarsening");

  // Do some other stuff before we start uncoarsening

  std::unique_ptr<IRefiner> lp_refiner;
  if ( context.refinement.use_lp_refiner ) {
    lp_refiner = std::make_unique<LabelPropagationRefiner>(graph, context);
  } else {
    lp_refiner = std::make_unique<DoNothingRefiner>(graph, context);
  }

  io::printUncoarseningBanner(context);
  utils::Timer::instance().start_timer("uncoarsening", "Unoarsening");
  coarsener->uncoarsen(lp_refiner);
  utils::Timer::instance().stop_timer("uncoarsening");
}

} // namespace multilevel
} // namespace cluster_editing