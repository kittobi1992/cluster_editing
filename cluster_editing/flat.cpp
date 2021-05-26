#include "flat.h"

#include "cluster_editing/macros.h"
#include "cluster_editing/refinement/evolutionary.h"
#include "cluster_editing/refinement/localized_evo.h"
#include "cluster_editing/refinement/lp_refiner.h"
#include "cluster_editing/refinement/fm_refiner.h"
#include "cluster_editing/refinement/stopping_rule.h"
#include "cluster_editing/datastructures/fast_reset_flag_array.h"
#include "cluster_editing/utils/timer.h"
#include "cluster_editing/utils/randomize.h"
#include "cluster_editing/metrics.h"
#include "cluster_editing/io/output.h"

namespace cluster_editing {
namespace flat {


bool isSpecialInstance(const Graph& graph) {
  if ( graph.maxDegree() < 20 && graph.numEdges() / 2 > 1000000 ) {
    utils::CommonOperations::instance(graph).computeClusterSizes(graph);
    std::vector<NodeID> cluster_sizes =
            utils::CommonOperations::instance(graph)._cluster_sizes;
    std::sort(cluster_sizes.begin(), cluster_sizes.end(), std::greater<NodeID>());
    while ( cluster_sizes.back() == 0 ) {
      cluster_sizes.pop_back();
    }
    size_t percentile = 0.99 * cluster_sizes.size();
    const NodeID p_99 = cluster_sizes[cluster_sizes.size() - percentile];
    const NodeID max_cluster_size = cluster_sizes[0];
    return max_cluster_size <= 5 && p_99 <= 3;
  }
  return false;
}

void solve(Graph& graph, const Context& context) {

  io::printFlatBanner(context);

  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();

  EdgeWeight current_edits = static_cast<EdgeWeight>(graph.numEdges()) / 2;

  if ( context.general.read_from_file ) {
    io::readSolutionFile(graph, context.general.output_file);
    current_edits = metrics::edits(graph);
  }

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
  Evolutionary evo(graph, context);
  if ( context.refinement.use_evo ) {
    evo.initialize(graph);
    current_edits = evo.createInitialPopulation(graph, current_edits);
  } else {
    evo.setDone();
  }

  if ( current_edits != graph.numEdges()/2 ) {
    // have computed an initial solution --> do special instance check
    utils::CommonOperations::instance(graph)._is_special_instance = isSpecialInstance(graph);
  }

  // Localized Evolutionary
  LocalizedEvolutionary localized_evo(graph, context);
  if (!context.refinement.use_localized_evo) {
    localized_evo.setDone();
  }

  double burst_time_limit = 20;   // 20 seconds
  EdgeWeight evo_edits = -1, localized_evo_edits = -1;
  double evo_improvement_per_time = 0.0, localized_evo_improvement_per_time = 0.0;

  // TODO reset measurements at some point (in case either one makes 0 improvement)

  while (!context.isTimeLimitReached() && (!evo.done() || !localized_evo.done())) {
    if ( context.refinement.use_evo
        && (evo_edits == -1 || evo_improvement_per_time > localized_evo_improvement_per_time
                            || !context.refinement.use_localized_evo)) {
      auto start_time = std::chrono::high_resolution_clock::now();
      evo_edits = evo.performTimeLimitedEvoSteps(graph, burst_time_limit, current_edits);
      double elapsed = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();
      evo_improvement_per_time = double(current_edits - evo_edits) / elapsed;
      current_edits = evo_edits;
      if ( context.general.verbose_output)
        LOG << V(elapsed) << V(evo_improvement_per_time) << V(current_edits);
    }

    if ( context.refinement.use_localized_evo
         && (localized_evo_edits == -1 || evo_improvement_per_time <= localized_evo_improvement_per_time
                                        || !context.refinement.use_evo)) {
      auto start_time = std::chrono::high_resolution_clock::now();
      localized_evo.initialize(graph);
      localized_evo_edits = localized_evo.performTimeLimitedEvoSteps(graph, burst_time_limit, current_edits);
      double elapsed = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();
      localized_evo_improvement_per_time = double(current_edits - localized_evo_edits) / elapsed;
      current_edits = localized_evo_edits;
      if ( context.general.verbose_output)
        LOG << V(elapsed) << V(localized_evo_improvement_per_time) << V(current_edits);
    }
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