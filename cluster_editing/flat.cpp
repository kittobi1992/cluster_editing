#include "flat.h"

#include "cluster_editing/macros.h"
#include "cluster_editing/refinement/evolutionary.h"
#include "cluster_editing/refinement/lp_refiner.h"
#include "cluster_editing/refinement/fm_refiner.h"
#include "cluster_editing/refinement/exact_refiner.h"
#include "cluster_editing/refinement/stopping_rule.h"
#include "cluster_editing/utils/timer.h"
#include "cluster_editing/utils/randomize.h"
#include "cluster_editing/metrics.h"
#include "cluster_editing/io/output.h"

namespace cluster_editing {
namespace flat {

void solve(Graph& graph, const Context& context) {

  io::printFlatBanner(context);

  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();

  EdgeWeight current_edits = static_cast<EdgeWeight>(graph.numEdges()) / 2;

  // Label Propagation
  if ( context.refinement.use_lp_refiner ) {
    LabelPropagationRefiner lp_refiner(graph, context);
    lp_refiner.initialize(graph);
    current_edits = lp_refiner.refine(graph, current_edits);
  }

  // Boundary FM
  if ( context.refinement.use_boundary_fm_refiner ) {
    FMRefiner<FruitlessMovesStoppingRule> fm_refiner(graph, context, FMType::boundary);
    fm_refiner.initialize(graph);
    current_edits = fm_refiner.refine(graph, current_edits);
  }

  // Localized FM
  if ( context.refinement.use_localized_fm_refiner ) {
    FMRefiner<AdaptiveStoppingRule> fm_refiner(graph, context, FMType::localized);
    fm_refiner.initialize(graph);
    current_edits = fm_refiner.refine(graph, current_edits);
  }

  // Evolutionary
  if ( context.refinement.use_evo ) {
    Evolutionary evo(graph, context);
    evo.initialize(graph);
    current_edits = evo.refine(graph, current_edits);
  }

  // Exact Refiner
  if ( context.refinement.use_exact_refiner ) {
    ExactRefiner exact(graph, context);
    exact.initialize(graph);
    current_edits = exact.refine(graph, current_edits);
  }

  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_seconds(end - start);
  if ( context.general.verbose_output ) {
    io::printObjectives(graph, elapsed_seconds);
    io::printCliqueInfo(graph, context);
  }
}

} // namespace flat
} // namespace cluster_editing