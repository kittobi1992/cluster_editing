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
  _step = 0;
}

EdgeWeight Evolutionary::refineImpl(Graph& graph,
                                    const EdgeWeight current_edits,
                                    const EdgeWeight) {
  utils::Timer::instance().start_timer("evolutionary", "Evolutionary");
  if ( _context.isTimeLimitReached() ) {
    return current_edits;
  }
  _start = std::chrono::high_resolution_clock::now();

  // Create Initial Population
  utils::Timer::instance().start_timer("initial_population", "Create Initial Population");
  createInitialPopulation(graph, current_edits);
  utils::Timer::instance().stop_timer("initial_population");

  sortSolutions();
  if ( _context.isTimeLimitReached() ) {
    return _population[0].edits;
  }

  if ( _context.general.verbose_output ) LOG << "Evolutionary Algorithm:";
  performTimeLimitedEvoSteps(graph, _context.refinement.evo.time_limit, current_edits);
  utils::Timer::instance().stop_timer("evolutionary");
  return _population[0].edits;
}

EdgeWeight Evolutionary::performTimeLimitedEvoSteps(Graph& graph, double time_limit, EdgeWeight current_edits) {
  utils::Timer::instance().start_timer("evo_steps", "Evolutionary Steps");
  if (current_edits < _population[0].edits) {
    storeSolution(graph, *_population[0].solution);
    _population[0].edits = current_edits;
  }
  _is_special_instance = utils::CommonOperations::instance(graph)._is_special_instance;
  utils::ProgressBar evo_progress(
          _max_steps, _population[0].edits,
          _context.general.verbose_output && !debug && !_context.refinement.evo.enable_detailed_output);
  if (evo_progress.isEnabled()) {
    evo_progress += _step;
  }
  auto start_time = std::chrono::high_resolution_clock::now();
  for ( ; _step < _max_steps; ++_step ) {
    evolutionaryStep(graph);
    if ( evo_progress.isEnabled() ) {
      sortSolutions();
      evo_progress.setObjective(_population[0].edits);
      evo_progress += 1;
    }

    if (_step == static_cast<size_t>(_context.refinement.evo.enable_all_mutations_after_steps) ) {
      const std::string enable_all_mutations(
              static_cast<int>(Mutation::NUM_MUTATIONS), '1');
      _mutator.activateMutations(enable_all_mutations);
    }

    std::chrono::duration<double> elapsed_seconds(std::chrono::high_resolution_clock::now() - start_time);
    if ( _context.isTimeLimitReached() || elapsed_seconds.count() > time_limit ) {
      break;
    }
  }

  evo_progress += (_max_steps - evo_progress.count());

  // Apply best solution to graph
  sortSolutions();
  applySolution(graph, *_population[0].solution);
  ASSERT(metrics::edits(graph) == _population[0].edits);
  utils::Timer::instance().stop_timer("evo_steps");
  return _population[0].edits;
}

EdgeWeight Evolutionary::createInitialPopulation(Graph& graph, const EdgeWeight current_edits) {
  for ( int i = 0; i < _context.refinement.evo.solution_pool_size; ++i ) {
    _population[i].edits = refineSolution(graph, current_edits,
      _context.refinement.evo.initial_lp_iterations,
      _context.refinement.evo.initial_node_swapper_iterations,
      _context.refinement.evo.use_random_node_ordering,
      _context.general.verbose_output,
      0, true);
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
  return _population[0].edits;
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
  graph.checkpoint(_population[i].edits);
}

EdgeWeight Evolutionary::intensivate(Graph& graph, SolutionStats& stats) {
  // Refine
  const EdgeWeight initial_edits = stats.edits;
  stats.edits = refineSolution(graph, stats.edits,
    _context.refinement.evo.intensivate_lp_iterations,
    _context.refinement.evo.node_swapper_iterations,
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
    _context.refinement.evo.node_swapper_iterations,
    _context.refinement.evo.use_random_node_ordering,
    _show_detailed_output,
    initial_edits);
  _mutator.notifyResult(mutation, initial_edits, edits);
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
                                        const int node_swapper_iterations,
                                        const bool use_random_node_order,
                                        const bool show_detailed_output,
                                        const EdgeWeight target_edits,
                                        const bool is_initial_partition) {
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
  if ( show_detailed_output ) LOG << "Label Propagation Refiner:";
  EdgeWeight edits = _lp_refiner.refine(graph, current_edits, target_edits);
  const bool was_aborted = utils::CommonOperations::instance(graph)._lp_aborted_flag;

  if ( !was_aborted ) {
    if ( _is_special_instance || is_initial_partition ) {
      // Node Swapper Refiner
      _context.refinement.evo.node_swapper_iterations = node_swapper_iterations;
      _node_swapper.initialize(graph);
      if ( show_detailed_output ) LOG << "Node Swapper Refiner:";
      edits = _node_swapper.refine(graph, edits, target_edits);
      _context.refinement.evo.node_swapper_iterations =
        _original_context.refinement.evo.node_swapper_iterations;
    }

    if ( !_is_special_instance && ( is_initial_partition || utils::Randomize::instance().flipCoin() ) ) {
      // Clique Remover Refiner
      _clique_remover.initialize(graph);
      if ( show_detailed_output ) LOG << "Clique Remover Refiner:";
      edits = _clique_remover.refine(graph, edits, target_edits);

      // Clique Splitter Refiner
      _clique_splitter.initialize(graph);
      if ( show_detailed_output ) LOG << "Clique Splitter Refiner:";
      edits = _clique_splitter.refine(graph, edits, target_edits);
    }
  }

  _context.general.verbose_output = _original_context.general.verbose_output;
  return edits;
}

} // namespace cluster_editing