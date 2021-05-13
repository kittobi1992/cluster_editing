#include "evolutionary.h"

#include "cluster_editing/io/output.h"
#include "cluster_editing/utils/timer.h"
#include "cluster_editing/utils/progress_bar.h"
#include "cluster_editing/metrics.h"

namespace cluster_editing {

void Evolutionary::initializeImpl(Graph& graph) {
  for ( int i = 0; i < _context.refinement.evo.solution_pool_size; ++i ) {
    _solutions[i].assign(graph.numNodes(), INVALID_CLIQUE);
    _population[i].edits = graph.numEdges();
    _population[i].solution = &_solutions[i];
  }
}

EdgeWeight Evolutionary::refineImpl(Graph& graph,
                                    const EdgeWeight current_edits,
                                    const EdgeWeight) {
  utils::Timer::instance().start_timer("evolutionary", "Evolutionary");
  if ( _context.isTimeLimitReached() ) {
    return 0;
  }

  // Create Initial Population
  utils::Timer::instance().start_timer("initial_population", "Create Initial Population");
  createInitialPopulation(graph, current_edits);
  utils::Timer::instance().stop_timer("initial_population");

  sortSolutions();
  utils::ProgressBar evo_progress(
    _context.refinement.evo.evolutionary_steps, _population[0].edits,
    _context.general.verbose_output && !debug && !_context.refinement.evo.enable_detailed_output);
  utils::Timer::instance().start_timer("evo_steps", "Evolutionary Steps");
  for ( int i = 0; i < _context.refinement.evo.evolutionary_steps; ++i ) {
    evolutionaryStep(graph);
    if ( evo_progress.isEnabled() ) {
      sortSolutions();
      evo_progress.setObjective(_population[0].edits);
      evo_progress += 1;
    }

    if ( i == _context.refinement.evo.enable_all_mutations_after_steps ) {
      const std::string enable_all_mutations(
        static_cast<int>(Mutation::NUM_MUTATIONS), '1');
      _mutator.activateMutations(enable_all_mutations);
    }

    if ( _context.isTimeLimitReached() ) {
      break;
    }
  }
  utils::Timer::instance().stop_timer("evo_steps");
  evo_progress += (_context.refinement.evo.evolutionary_steps - evo_progress.count());

  // Apply best solution to graph
  sortSolutions();
  applySolution(graph, *_population[0].solution);
  ASSERT(metrics::edits(graph) == _population[0].edits);

  utils::Timer::instance().stop_timer("evolutionary");
  return _population[0].edits;
}

void Evolutionary::createInitialPopulation(Graph& graph, const EdgeWeight current_edits) {
  for ( int i = 0; i < _context.refinement.evo.solution_pool_size; ++i ) {
    _population[i].edits = refineSolution(graph, current_edits,
      _context.refinement.evo.initial_lp_iterations,
      _context.refinement.evo.use_random_node_ordering,
      _context.general.verbose_output);
    storeSolution(graph, *_population[i].solution);
  }
  sortSolutions();

  if ( _show_detailed_output ) {
    std::cout << "Create Initial Solution Pool:";
    for ( size_t i = 0; i < _population.size(); ++i ) {
      std::cout << " " << _population[i].edits;
    }
    std::cout << std::endl;
  }
}

void Evolutionary::evolutionaryStep(Graph& graph) {
  // Pick random solution
  int i = utils::Randomize::instance().getRandomInt(0, _population.size() - 1);
  applySolution(graph, *_population[i].solution);
  const EvoAction action = _evo_action_selector.chooseAction(_show_detailed_output);

  EdgeWeight delta = 0;
  if ( action == EvoAction::INTESIVATE ) {
    if ( _show_detailed_output ) LOG << "Evolutionary Action: INTENSIVATE";
    delta = intensivate(graph, _population[i]);
  } else if ( action == EvoAction::MUTATE ) {
    if ( _show_detailed_output ) LOG << "Evolutionary Action: MUTATE";
    delta = mutate(graph, _population[i]);
  }
  _evo_action_selector.notifyImprovement(action, delta);
}

EdgeWeight Evolutionary::intensivate(Graph& graph, SolutionStats& stats) {
  // Refine
  const EdgeWeight initial_edits = stats.edits;
  stats.edits = refineSolution(graph, stats.edits,
    _context.refinement.evo.intensivate_lp_iterations,
    _context.refinement.evo.use_random_node_ordering,
    _show_detailed_output);
  storeSolution(graph, *stats.solution);
  if ( stats.edits < initial_edits ) {
    if ( _show_detailed_output ) {
      LOG << GREEN << "INTENSIVATE: Improve solution from"
          << initial_edits << "to" << stats.edits
          << "( Delta: " << (stats.edits - initial_edits)
          << ")" << END;
    }
  }
  return std::max(0, initial_edits - stats.edits);
}

EdgeWeight Evolutionary::mutate(Graph& graph, SolutionStats& stats) {
  const EdgeWeight initial_edits = stats.edits;

  // Mutate
  const Mutation mutation = _mutator.mutate(graph);
  if ( _show_detailed_output ) {
    const EdgeWeight mutate_edits = metrics::edits(graph);
    LOG << "Mutation changed solution from" << initial_edits
        << "to" << mutate_edits << "edits";
  }

  // Refine
  _context.refinement.lp.maximum_lp_iterations =
    _context.refinement.evo.lp_iterations_after_mutate;
  const EdgeWeight edits = refineSolution(graph,
    metrics::edits(graph),
    _context.refinement.evo.lp_iterations_after_mutate,
    _context.refinement.evo.use_random_node_ordering,
    _show_detailed_output,
    initial_edits);
  _mutator.updateProbs(mutation, initial_edits, edits);
  if ( edits < initial_edits ) {
    if ( _show_detailed_output ) {
      LOG << GREEN << "MUTATE: Improved solution quality from"
          << initial_edits << "to" << edits << "( Delta:"
          << (edits - initial_edits) << ")" << END;
    }
    evictSolution(graph, edits);
  } else if ( _show_detailed_output && edits > initial_edits ) {
    LOG << RED << "MUTATE: Worsen solution quality (Initial:"
        << initial_edits << ", After: " << edits << ", Delta:"
        << (edits - initial_edits) << ")" << END;
  }
  return std::max(0, initial_edits - edits);
}

EdgeWeight Evolutionary::refineSolution(Graph& graph,
                                        const EdgeWeight current_edits,
                                        const int lp_iterations,
                                        const bool use_random_node_order,
                                        const bool show_detailed_output,
                                        const EdgeWeight target_edits) {
  _lp_refiner.initialize(graph);
  // Choose some random node ordering for LP refiner
  _context.refinement.lp.maximum_lp_iterations = lp_iterations;
  if ( use_random_node_order ) {
    _context.refinement.lp.node_order =
      static_cast<NodeOrdering>(utils::Randomize::instance().getRandomInt(0, 3));
    } else {
    _context.refinement.lp.node_order = _original_context.refinement.lp.node_order;
  }
  _context.general.verbose_output = show_detailed_output;
  const EdgeWeight edits = _lp_refiner.refine(graph, current_edits, target_edits);
  _context.general.verbose_output = _original_context.general.verbose_output;
  return edits;
}

} // namespace cluster_editing