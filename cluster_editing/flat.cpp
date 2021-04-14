#include "flat.h"

#include "cluster_editing/macros.h"
#include "cluster_editing/refinement/lp_refiner.h"
#include "cluster_editing/refinement/fm_refiner.h"
#include "cluster_editing/refinement/stopping_rule.h"
#include "cluster_editing/utils/timer.h"
#include "cluster_editing/io/output.h"

namespace cluster_editing {
namespace flat {

void solve(Graph& graph, const Context& context) {

  io::printFlatBanner(context);

  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();

  // Label Propagation
  if ( context.refinement.use_lp_refiner ) {
    LabelPropagationRefiner lp_refiner(graph, context);
    lp_refiner.initialize(graph);
    lp_refiner.refine(graph);
  }

  // Boundary FM
  if ( context.refinement.use_boundary_fm_refiner ) {
    FMRefiner<FruitlessMovesStoppingRule> fm_refiner(graph, context, FMType::boundary);
    fm_refiner.initialize(graph);
    fm_refiner.refine(graph);
  }

  // Localized FM
  if ( context.refinement.use_localized_fm_refiner ) {
    FMRefiner<AdaptiveStoppingRule> fm_refiner(graph, context, FMType::localized);
    fm_refiner.initialize(graph);
    fm_refiner.refine(graph);
  }

  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_seconds(end - start);
  if ( context.general.verbose_output ) {
    io::printObjectives(graph, elapsed_seconds);
  }
}

} // namespace flat
} // namespace cluster_editing