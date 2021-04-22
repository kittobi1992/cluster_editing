#include "flat.h"

#include "cluster_editing/macros.h"
#include "cluster_editing/refinement/evolutionary.h"
#include "cluster_editing/refinement/lp_refiner.h"
#include "cluster_editing/refinement/fm_refiner.h"
#include "cluster_editing/refinement/stopping_rule.h"
#include "cluster_editing/utils/timer.h"
#include "cluster_editing/utils/randomize.h"
#include "cluster_editing/metrics.h"
#include "cluster_editing/io/output.h"

namespace cluster_editing {
namespace flat {

namespace {

  struct InitialSolution {
    EdgeWeight edits;
    std::vector<CliqueID>* solution;
  };

  void initial_partition(Graph& graph, Context context) {
    std::vector<InitialSolution> initial_solutions(
      context.refinement.lp.maximum_lp_iterations, InitialSolution { static_cast<EdgeWeight>(graph.numEdges()), nullptr });
    std::vector<std::vector<CliqueID>> initial_cliques(context.refinement.lp.maximum_lp_iterations);

    auto store_cliques = [&](std::vector<CliqueID>& clique) {
      clique.resize(graph.numNodes());
      for ( const NodeID& u : graph.nodes() ) {
        clique[u] = graph.clique(u);
      }
    };

    auto partition = [&](const int i) {
      graph.reset();
      if ( initial_solutions[i].solution != nullptr ) {
        const std::vector<CliqueID>& clique = *initial_solutions[i].solution;
        for ( const NodeID& u : graph.nodes() ) {
          graph.setClique(u, clique[u]);
        }
      } else {
        initial_solutions[i].solution = &initial_cliques[i];
      }

      // Choose some random node ordering for LP refiner
      context.refinement.lp.node_order =
        static_cast<NodeOrdering>(utils::Randomize::instance().getRandomInt(0, 3));

      // Perform LP Refinement
      LabelPropagationRefiner lp_refiner(graph, context);
      lp_refiner.initialize(graph);
      initial_solutions[i].edits = lp_refiner.refine(graph);
      store_cliques(*initial_solutions[i].solution);
    };

    // Create Initial Solution Pool
    int solution_pool_size = context.refinement.evo.solution_pool_size;
    int current_lp_iterations = context.refinement.evo.lp_iterations;
    while ( solution_pool_size > 0 ) {
      context.refinement.lp.maximum_lp_iterations = current_lp_iterations;
      for ( int i = 0; i < solution_pool_size; ++i ) {
        partition(i);
      }
      std::sort(initial_solutions.begin(), initial_solutions.end(),
        [&](const InitialSolution& sol_1, const InitialSolution& sol_2) {
          return sol_1.edits < sol_2.edits;
        });
      solution_pool_size = ( solution_pool_size == 1 ? 0 : ( solution_pool_size / 2 ) + ( solution_pool_size % 2 != 0) );
    }

    // Apply best solution
    if ( context.general.verbose_output ) {
      LOG << GREEN << "Best initial partition has" << initial_solutions[0].edits << "edits" << END;
    }
    const std::vector<CliqueID>& best_clique = *initial_solutions[0].solution;
    for ( const NodeID& u : graph.nodes() ) {
      graph.setClique(u, best_clique[u]);
    }
  }
}

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

  // Evolutionary
  if ( context.refinement.use_evo ) {
    Evolutionary evo(graph, context);
    evo.initialize(graph);
    evo.refine(graph);
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