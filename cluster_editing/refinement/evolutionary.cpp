#include "evolutionary.h"

#include "cluster_editing/io/output.h"
#include "cluster_editing/utils/timer.h"
#include "cluster_editing/utils/progress_bar.h"
#include "cluster_editing/metrics.h"
#include "cluster_editing/refinement/mutation.h"

namespace cluster_editing {

void Evolutionary::initializeImpl(Graph& graph) {
  for ( int i = 0; i < _context.refinement.evo.solution_pool_size; ++i ) {
    _solutions[i].assign(graph.numNodes(), INVALID_CLIQUE);
    _population[i].edits = graph.numEdges();
    _population[i].solution = &_solutions[i];
  }
}

EdgeWeight Evolutionary::refineImpl(Graph& graph) {
  utils::Timer::instance().start_timer("evolutionary", "Evolutionary");

  // Create Initial Population
  utils::Timer::instance().start_timer("initial_population", "Create Initial Population");
  createInitialPopulation(graph);
  utils::Timer::instance().stop_timer("initial_population");

  utils::Timer::instance().start_timer("evo_steps", "Evolutionary Steps");
  for ( int i = 0; i < _context.refinement.evo.evolutionary_steps; ++i ) {
    evolutionaryStep(graph);
  }
  utils::Timer::instance().stop_timer("evo_steps");

  // Apply best solution to graph
  sortSolutions();
  applySolution(graph, *_population[0].solution);
  ASSERT(metrics::edits(graph) == _population[0].edits);

  utils::Timer::instance().stop_timer("evolutionary");
  return 0;
}

void Evolutionary::createInitialPopulation(Graph& graph) {
  for ( int i = 0; i < _context.refinement.evo.solution_pool_size; ++i ) {
    graph.reset();
    _population[i].edits = refineSolution(graph);
    storeSolution(graph, *_population[i].solution);
  }
  sortSolutions();

  if ( _context.general.verbose_output ) {
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
  const Action action = selectAction(_population[i].no_intensivate);

  if ( action == Action::INTESIVATE ) {
    if ( _context.general.verbose_output ) LOG << "Evolutionary Action: INTENSIVATE";
    intensivate(graph, _population[i]);
  } else if ( action == Action::MUTATE ) {
    if ( _context.general.verbose_output ) LOG << "Evolutionary Action: MUTATE";
    mutate(graph, _population[i]);
  } else if ( action == Action::COMBINE ) {
    if ( _context.general.verbose_output ) LOG << "Evolutionary Action: COMBINE";
  }
}

void Evolutionary::intensivate(Graph& graph, SolutionStats& stats) {
  // Refine
  const EdgeWeight initial_edits = stats.edits;
  _context.refinement.lp.maximum_lp_iterations =
    _context.refinement.evo.lp_iterations;
  stats.edits = refineSolution(graph);
  storeSolution(graph, *stats.solution);
  if ( stats.edits < initial_edits ) {
    stats.reset();
    if ( _context.general.verbose_output ) {
      LOG << GREEN << "INTENSIVATE: Improve solution from"
          << initial_edits << "to" << stats.edits
          << "( Delta: " << (stats.edits - initial_edits)
          << ")" << END;
    }
  } else {
    stats.no_intensivate = true;
  }
}

namespace {

enum class Mutation : uint8_t {
  LARGE_CLIQUE_ISOLATOR = 0,
  LARGE_CLIQUE_WITH_NEIGHBOR_ISOLATOR = 1,
  RANDOM_NODE_ISOLATOR = 2,
  RANDOM_NODE_MOVER = 3,
  NUM_MUTATIONS = 4
};

} // namespace

void Evolutionary::mutate(Graph& graph, SolutionStats& stats) {
  // Mutate
  const EdgeWeight initial_edits = stats.edits;
  const Mutation mutation = static_cast<Mutation>(
    utils::Randomize::instance().getRandomInt(0,
      static_cast<int>(Mutation::NUM_MUTATIONS) - 1));
  if ( mutation == Mutation::LARGE_CLIQUE_ISOLATOR ) {
    if ( _context.general.verbose_output )
      LOG << "Mutation Action: LARGE_CLIQUE_ISOLATOR";
    LargeCliqueIsolator::mutate(graph, _context);
  } else if ( mutation == Mutation::LARGE_CLIQUE_WITH_NEIGHBOR_ISOLATOR ) {
    if ( _context.general.verbose_output )
      LOG << "Mutation Action: LARGE_CLIQUE_WITH_NEIGHBOR_ISOLATOR";
    LargeCliqueWithNeighborIsolator::mutate(graph, _context);
  } else if ( mutation == Mutation::RANDOM_NODE_ISOLATOR ) {
    if ( _context.general.verbose_output )
      LOG << "Mutation Action: RANDOM_NODE_ISOLATOR";
    RandomNodeIsolator::mutate(graph, _context);
  } else if ( mutation == Mutation::RANDOM_NODE_MOVER ) {
    if ( _context.general.verbose_output )
      LOG << "Mutation Action: RANDOM_NODE_MOVER";
    RandomNodeMover::mutate(graph, _context);
  }

  if ( _context.general.verbose_output ) {
    const EdgeWeight mutate_edits = metrics::edits(graph);
    LOG << "Mutation changed solution from" << initial_edits
        << "to" << mutate_edits << "edits";
  }

  // Refine
  _context.refinement.lp.maximum_lp_iterations =
    _context.refinement.evo.lp_iterations_after_mutate;
  const EdgeWeight edits = refineSolution(graph);
  if ( edits < initial_edits ) {
    if ( _context.general.verbose_output ) {
      LOG << GREEN << "MUTATE: Improved solution quality from"
          << initial_edits << "to" << edits << "( Delta:"
          << (edits - initial_edits) << ")" << END;
    }
    evictSolution(graph, edits);
  } else if ( _context.general.verbose_output && edits > initial_edits ) {
    LOG << RED << "MUTATE: Worsen solution quality (Initial:"
        << initial_edits << ", After: " << edits << ", Delta:"
        << (edits - initial_edits) << ")" << END;
  }
}

EdgeWeight Evolutionary::refineSolution(Graph& graph) {
  _lp_refiner.initialize(graph);
  // Choose some random node ordering for LP refiner
  _context.refinement.lp.node_order =
    static_cast<NodeOrdering>(utils::Randomize::instance().getRandomInt(0, 3));
  return _lp_refiner.refine(graph);
}

} // namespace cluster_editing