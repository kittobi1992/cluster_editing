/*******************************************************************************
 * This file is part of KaPoCE.
 *
 * KaPoCE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaPoCE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaPoCE.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "clustering.h"

#include "cluster_editing/macros.h"
#include "cluster_editing/heuristic/evolutionary.h"
#include "cluster_editing/heuristic/localized_evo.h"
#include "cluster_editing/heuristic/lp_refiner.h"
#include "cluster_editing/heuristic/fm_refiner.h"
#include "cluster_editing/heuristic/stopping_rule.h"
#include "cluster_editing/datastructures/fast_reset_flag_array.h"
#include "cluster_editing/utils/timer.h"
#include "cluster_editing/utils/randomize.h"
#include "cluster_editing/metrics.h"
#include "cluster_editing/io/output.h"

namespace cluster_editing {

namespace {

bool isSpecialInstance(const Graph& graph) {
  if ( graph.maxDegree() < 20 && graph.numEdges() / 2 > 1000000 ) {
    utils::CommonOperations::instance(graph).computeClusterSizes(graph);
    std::vector<NodeID> cluster_sizes =
            utils::CommonOperations::instance(graph)._cluster_sizes;
    std::sort(cluster_sizes.begin(), cluster_sizes.end(), std::greater<NodeID>());
    while ( cluster_sizes.back() == 0 ) {
      cluster_sizes.pop_back();
    }
    size_t percentile_75 = 0.75 * cluster_sizes.size();
    size_t percentile_99 = 0.99 * cluster_sizes.size();
    const NodeID p_75 = cluster_sizes[cluster_sizes.size() - percentile_75];
    const NodeID p_99 = cluster_sizes[cluster_sizes.size() - percentile_99];
    const NodeID max_cluster_size = cluster_sizes[0];
    return p_75 <= 2 && p_99 <= 4 && max_cluster_size <= 6;
  }
  return false;
}

void printEvoProgress(const Context& context,
                      const double last_burst_time,
                      const double evo_improvement_per_seconds,
                      const double localized_evo_improvement_per_seconds,
                      const EdgeWeight last_evo_delta,
                      const EdgeWeight last_localized_evo_delta) {
  if ( context.general.verbose_output ) {
    LOG << "Last Burst Elapsed Time (seconds)                =" << last_burst_time;
    LOG << "Evo. Improvement Per Sec. (last burst)           ="
        << evo_improvement_per_seconds << "( Last Delta =" << last_evo_delta << ")";
    LOG << "Localized Evo. Improvement Per Sec. (last burst) ="
        << localized_evo_improvement_per_seconds<< "( Last Delta =" << last_localized_evo_delta << ")";
  }
}

} // namespace

void solve(Graph& graph, const Context& context) {

  io::printInitialSolutionBanner(context);

  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();

  EdgeWeight current_edits = static_cast<EdgeWeight>(graph.numEdges()) / 2;

  if ( context.general.read_from_file && context.general.graph_filename != "" ) {
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
    HighResClockTimepoint start_init = std::chrono::high_resolution_clock::now();
    HighResClockTimepoint end_init = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds(end_init - start_init);
    for ( size_t i = 0; i < 100 && elapsed_seconds.count() <= 10 && !context.isTimeLimitReached(); ++i ) {
      graph.reset();
      current_edits = static_cast<EdgeWeight>(graph.numEdges()) / 2;
      evo.initialize(graph);
      current_edits = evo.createInitialPopulation(graph, current_edits);
      graph.checkpoint(current_edits);
      end_init = std::chrono::high_resolution_clock::now();
      elapsed_seconds = end_init - start_init;
      if ( context.general.verbose_output ) io::printStripe();
    }
    graph.applyBestCliques();
    current_edits = graph.bestEdits();
  } else {
    evo.setDone();
  }

  if ( current_edits != EdgeWeight(graph.numEdges()/2) ) {
    // have computed an initial solution --> do special instance check
    utils::CommonOperations::instance(graph)._is_special_instance = isSpecialInstance(graph);
  }

  // Localized Evolutionary
  LocalizedEvolutionary localized_evo(graph, context);
  if (!context.refinement.use_localized_evo) {
    localized_evo.setDone();
  }

  double burst_time_limit = 20;   // 20 seconds
  EdgeWeight evo_edits = -1;
  EdgeWeight evo_last_delta = 0;
  EdgeWeight localized_evo_edits = -1;
  EdgeWeight localized_evo_last_delta = 0;
  double evo_improvement_per_time = 0.0;
  double localized_evo_improvement_per_time = 0.0;
  size_t evo_not_run = 0;
  size_t localized_evo_not_run = 0;
  while (!context.isTimeLimitReached() && (!evo.done() || !localized_evo.done())) {
    if ( context.refinement.use_evo &&
         ( evo_edits == -1 ||
           evo_improvement_per_time > localized_evo_improvement_per_time ||
           !context.refinement.use_localized_evo)) {
      io::printGlobalEvoBanner(context);
      auto start_time = std::chrono::high_resolution_clock::now();
      evo_edits = evo.performTimeLimitedEvoSteps(graph, burst_time_limit, current_edits);
      double elapsed = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();
      evo_last_delta = current_edits - evo_edits;
      evo_improvement_per_time = static_cast<double>(evo_last_delta) / elapsed;
      current_edits = evo_edits;
      printEvoProgress(context, elapsed,
        evo_improvement_per_time, localized_evo_improvement_per_time,
        evo_last_delta, localized_evo_last_delta);
      evo_not_run = 0;
    } else {
      evo_not_run++;
    }

    if ( context.refinement.use_localized_evo &&
        ( localized_evo_edits == -1 ||
          evo_improvement_per_time < localized_evo_improvement_per_time ||
          !context.refinement.use_evo)) {
      io::printLocalizedEvoBanner(context);
      auto start_time = std::chrono::high_resolution_clock::now();
      localized_evo.initialize(graph);
      localized_evo_edits = localized_evo.performTimeLimitedEvoSteps(graph, burst_time_limit, current_edits);
      double elapsed = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count();
      localized_evo_last_delta = current_edits - localized_evo_edits;
      localized_evo_improvement_per_time = static_cast<double>(localized_evo_last_delta) / elapsed;
      current_edits = localized_evo_edits;
      printEvoProgress(context, elapsed,
        evo_improvement_per_time, localized_evo_improvement_per_time,
        evo_last_delta, localized_evo_last_delta);
      localized_evo_not_run = 0;
    } else {
      localized_evo_not_run++;
    }

    if (localized_evo_not_run == 7 && localized_evo_improvement_per_time < 0.05) {
      localized_evo_edits = -1;   // run again
    }
    if (evo_not_run == 7 && evo_improvement_per_time < 0.05) {
      evo_edits = -1;             // run again
    }

    if (evo_improvement_per_time == localized_evo_improvement_per_time) {
      // in case both had the same improvement (most likely zero -.-) randomly decide which one to run again. do global evo more often.
      if (utils::Randomize::instance().getRandomFloat(0.0, 1.0) < 0.75) {
        evo_edits = -1;
      } else {
        localized_evo_edits = -1;
      }
    }
  }

  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_seconds(end - start);
  if ( context.general.verbose_output ) {
    LOG << "";
    io::printStripe();
    LOG << "";
    io::printObjectives(graph, elapsed_seconds);
    io::printCliqueInfo(graph, context);
  }
}

} // namespace cluster_editing