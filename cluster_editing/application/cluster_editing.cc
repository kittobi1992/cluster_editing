#include <iostream>
#include <set>
#include <signal.h>

#include "cluster_editing/definitions.h"
#include "cluster_editing/preprocessing.h"
#include "cluster_editing/flat.h"
#include "cluster_editing/multilevel.h"
#include "cluster_editing/datastructures/spin_lock.h"
#include "cluster_editing/utils/timer.h"
#include "cluster_editing/utils/randomize.h"
#include "cluster_editing/metrics.h"
#include "cluster_editing/io/output.h"
#include "cluster_editing/io/graph_io.h"

using namespace cluster_editing;

SpinLock terminate_lock;
Context context;
Graph graph;
HighResClockTimepoint start, end;

void printResult(Graph& best) {
  // Print RESULT line
  if ( context.general.print_result_line ) {
    // Deletions
    std::vector<std::vector<NodeID>> cliques(best.numNodes());
    for ( const NodeID u : best.nodes() ) {
      cliques[best.clique(u)].push_back(u);
      for ( const auto& v : best.neighbors(u) ) {
        if ( u < v.target && best.clique(u) != best.clique(v.target) ) {
          std::cout << (u + 1) << " " << (v.target + 1) << std::endl;
        }
      }
    }

    // Insertions
    for ( CliqueID c = 0; c < best.numNodes(); ++c ) {
      if ( cliques[c].size() > 0 ) {
        std::set<NodeID> clique(cliques[c].begin(), cliques[c].end());

        for ( const NodeID u : cliques[c] ) {
          for ( const auto& v : best.neighbors(u) ) {
            if ( best.clique(u) == best.clique(v.target) ) {
              clique.erase(v.target);
            }
          }

          for ( const NodeID v : clique ) {
            if ( u < v ) {
              std::cout << (u + 1) << " " << (v + 1) << std::endl;
            }
          }

          for ( const auto& v : best.neighbors(u) ) {
            if ( best.clique(u) == best.clique(v.target) ) {
              clique.insert(v.target);
            }
          }
        }
      }
    }
  }

  // Print Stats
  std::chrono::duration<double> elapsed_seconds(end - start);
  io::printClusterEditingResults(best, context, elapsed_seconds);

  // Print RESULT line
  if ( context.general.print_result_line ) {
    io::printResultLine(best, context, elapsed_seconds);
  }

  if ( context.general.print_csv ) {
    std::cout << context.general.graph_filename.substr(context.general.graph_filename.find_last_of('/') + 1)
              << "," << metrics::edits(best) << "," << elapsed_seconds.count() << std::endl;
  }
}

void terminate(int) {
  if ( terminate_lock.tryLock() ) {
    end = std::chrono::high_resolution_clock::now();
    // Copy graph since refinement algorithm can checks
    // cliques while we try to output solution
    Graph cpy_graph = graph.copyBestSolution();
    printResult(cpy_graph);
    std::exit(0);
  }
}

int main() {
  // Register signal handler
  struct sigaction action;
  memset(&action, 0, sizeof(struct sigaction));
  action.sa_handler = terminate;
  sigaction(SIGTERM, &action, NULL);

  context.general.verbose_output = true;
  context.general.print_result_line = false;
  context.general.seed = 0;
  context.general.use_multilevel = false;
  context.general.num_repititions = 1;
  context.general.num_fruitless_repititions = 10;
  context.refinement.use_lp_refiner = true;
  context.refinement.lp.maximum_lp_iterations = 1000;
  context.refinement.lp.activate_all_cliques_after_rounds = 10;
  context.refinement.lp.random_shuffle_each_round = false;
  context.refinement.lp.node_order = NodeOrdering::none;
  context.refinement.lp.min_improvement = 1;
  context.refinement.use_fm_refiner = true;
  context.refinement.fm.maximum_fm_iterations = 25;
  context.refinement.fm.fraction_of_fruitless_moves = 0.05;
  utils::Randomize::instance().setSeed(context.general.seed);

  utils::Timer::instance().start_timer("import_graph", "Import Graph");
  graph = io::readGraphFile();
  utils::Timer::instance().stop_timer("import_graph");
  context.computeParameters(graph.numNodes());

  if ( context.general.verbose_output ) {
    io::printBanner();
    // Print context description
    LOG << context;
  }
  io::printInputInfo(graph, context);

  if ( graph.numNodes() <= 10000 ) {
    context.refinement.lp.node_order = NodeOrdering::random_shuffle;
    context.refinement.lp.random_shuffle_each_round = true;
  } else if ( graph.numNodes() <= 100000 ) {
    context.refinement.lp.node_order = NodeOrdering::random_shuffle;
    context.refinement.lp.random_shuffle_each_round = true;
  } else {
    context.refinement.lp.node_order = NodeOrdering::none;
    context.refinement.lp.random_shuffle_each_round = false;
  }

  // Preprocessing
  io::printPreprocessingBanner(context);
  start = std::chrono::high_resolution_clock::now();
  utils::Timer::instance().start_timer("preprocessing", "Preprocessing");
  Preprocessor preprocessor(graph, context);
  preprocessor.preprocess();
  utils::Timer::instance().stop_timer("preprocessing");

  // Multilevel Solver
  utils::Timer::instance().start_timer("solver", "Solver");
  int fruitless_repititions = 0;
  std::vector<CliqueID> best_cliques(graph.numNodes(), INVALID_CLIQUE);
  size_t best_objective = std::numeric_limits<size_t>::max();
  auto check_for_new_best_solution = [&] {
    const size_t current_objective = metrics::edits(graph);
    if ( current_objective < best_objective ) {
      if ( context.general.verbose_output ) {
        LOG << GREEN << "Improved best solution from"
            << best_objective << "to" << current_objective << END;
      }
      graph.checkpoint(current_objective);
      best_objective = current_objective;
      fruitless_repititions = 0;
    } else {
      ++fruitless_repititions;
    }
  };
  auto solve = [&] {
    graph.reset();
    if ( context.general.use_multilevel ) {
      multilevel::solve(graph, context);
    } else {
      flat::solve(graph, context);
    }

    // Check if solution is better than best solution found so far
    check_for_new_best_solution();
  };

  // First Iteration
  HighResClockTimepoint s = std::chrono::high_resolution_clock::now();
  solve();
  HighResClockTimepoint e = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> solve_time(e - s);
  double time = solve_time.count();
  if ( time <= 1.0 ) {
    context.general.num_repititions = 100;
  } else if ( time <= 2.0 ) {
    context.general.num_repititions = 25;
  } else if ( time <= 5.0 ) {
    context.general.num_repititions = 10;
  } else if ( time <= 10 ) {
    context.general.num_repititions = 5;
  } else {
    context.general.num_repititions = 1;
  }

  if ( context.general.num_repititions == 1 ) {
    context.refinement.lp.node_order = NodeOrdering::none;
    context.refinement.lp.random_shuffle_each_round = false;
  } else {
    context.refinement.lp.random_shuffle_each_round = true;
  }

  for ( int i = 0;
        i < context.general.num_repititions &&
        fruitless_repititions < context.general.num_fruitless_repititions ; ++i ) {
    context.refinement.lp.node_order =
      static_cast<NodeOrdering>(utils::Randomize::instance().getRandomInt(0,3));
    solve();
  }

  utils::Timer::instance().stop_timer("solver");

  // Undo Preprocessing
  io::printUndoPreprocessingBanner(context);
  utils::Timer::instance().start_timer("undo_preprocessing", "Undo Preprocessing");
  preprocessor.undoPreprocessing();
  utils::Timer::instance().stop_timer("undo_preprocessing");
  end = std::chrono::high_resolution_clock::now();

  if ( terminate_lock.tryLock() ) {
    graph.applyBestCliques();
    printResult(graph);
  }

  return 0;
}
