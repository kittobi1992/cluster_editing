#include "flat.h"

#include "cluster_editing/macros.h"
#include "cluster_editing/refinement/lp_refiner.h"
#include "cluster_editing/utils/timer.h"
#include "cluster_editing/io/output.h"

namespace cluster_editing {
namespace flat {

void solve(Graph& graph, const Context& context) {

  io::printFlatBanner(context);

  if ( context.refinement.use_lp_refiner ) {
    LabelPropagationRefiner lp_refiner(graph, context);
    lp_refiner.initialize(graph);
    lp_refiner.refine(graph);
  }

}

} // namespace flat
} // namespace cluster_editing