#include "flat.h"

#include "cluster_editing/macros.h"
#include "cluster_editing/refinement/lp_refiner.h"
#include "cluster_editing/refinement/fm_refiner.h"
#include "cluster_editing/refinement/stopping_rule.h"
#include "cluster_editing/utils/timer.h"
#include "cluster_editing/utils/randomize.h"
#include "cluster_editing/metrics.h"
#include "cluster_editing/io/output.h"

namespace cluster_editing {
namespace flat {

void solve(Graph& graph, const Context& context) {

  io::printFlatBanner(context);

  EdgeWeight best_metric = std::numeric_limits<EdgeWeight>::max();
  std::vector<CliqueID> best_cliques(graph.numNodes(), INVALID_CLIQUE);

  auto is_best_clique = [&](const EdgeWeight current_metric) {
    if ( current_metric < best_metric ) {
      best_metric = current_metric;
      for ( const NodeID& u : graph.nodes() ) {
        best_cliques[u] = graph.clique(u);
      }
      ASSERT(current_metric == metrics::edits(graph));
    }
  };

  // Label Propagation
  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  if ( context.refinement.use_lp_refiner ) {
    LabelPropagationRefiner lp_refiner(graph, context);
    lp_refiner.initialize(graph);
    const EdgeWeight current_metric = lp_refiner.refine(graph);
    is_best_clique(current_metric);
  }
  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> lp_time(end - start);

  // Spend more repititions for LP refiner, if it is fast (< 10 seconds)
  if ( context.refinement.use_lp_refiner && lp_time.count() < 10 ) {
    Context lp_context(context);
    for ( int i = 0; i < context.refinement.lp.maximum_lp_repititions - 1; ++i ) {
      // Choose some random node ordering for LP refiner
      lp_context.refinement.lp.node_order =
        static_cast<NodeOrdering>(utils::Randomize::instance().getRandomInt(0, 3));
      if ( lp_context.refinement.lp.node_order == NodeOrdering::random_shuffle ) {
        lp_context.refinement.lp.random_shuffle_each_round =
          utils::Randomize::instance().flipCoin();
      }

      // Perform LP refinement
      graph.reset();
      LabelPropagationRefiner lp_refiner(graph, lp_context);
      lp_refiner.initialize(graph);
      const EdgeWeight current_metric = lp_refiner.refine(graph);
      is_best_clique(current_metric);
    }
  }

  if ( best_metric != std::numeric_limits<EdgeWeight>::max() ) {
    for ( const NodeID& u : graph.nodes() ) {
      graph.setClique(u, best_cliques[u]);
    }
    if ( context.general.verbose_output ) {
      LOG << "Best LP Refiner Solution is" << GREEN << best_metric << END;
    }
    ASSERT(best_metric == metrics::edits(graph));
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

  end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_seconds(end - start);
  if ( context.general.verbose_output ) {
    io::printObjectives(graph, elapsed_seconds);
  }
}

} // namespace flat
} // namespace cluster_editing