#include <iostream>
#include <set>
#include <signal.h>

#include "cluster_editing/definitions.h"
#include "cluster_editing/flat.h"
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
  } else {
    // Print Stats
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds(end - start);
    io::printClusterEditingResults(best, context, elapsed_seconds);

    // Print RESULT line
    io::printResultLine(best, context, elapsed_seconds);

    if ( context.general.print_csv ) {
      std::cout << context.general.graph_filename.substr(context.general.graph_filename.find_last_of('/') + 1)
                << "," << metrics::edits(best) << "," << elapsed_seconds.count() << std::endl;
    }
  }
}

void terminate(int) {
  if ( terminate_lock.tryLock() ) {
    end = std::chrono::high_resolution_clock::now();
    // Copy graph since refinement algorithm can change
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

  // General Options
  context.general.verbose_output = false;
  context.general.print_result_line = true;
  context.general.seed = 2;
  context.general.num_repititions = 1;
  context.general.num_fruitless_repititions = 10;
  utils::Randomize::instance().setSeed(context.general.seed);

  // Initial Partitioning Options
  context.refinement.use_ip = true;
  context.refinement.ip.initial_solution_pool_size = 5;
  context.refinement.ip.initial_lp_iterations = 25;
  context.refinement.ip.scale_lp_iteration_factor = 2.0;

  // LP Refiner Options
  context.refinement.use_lp_refiner = true;
  context.refinement.lp.maximum_lp_iterations = 1000;
  context.refinement.lp.random_shuffle_each_round = false;
  context.refinement.lp.node_order = NodeOrdering::none;
  context.refinement.lp.min_improvement = 5;
  context.refinement.lp.early_exit_window = 100;

  // FM Refiner Options
  context.refinement.use_boundary_fm_refiner = true;
  context.refinement.use_localized_fm_refiner = false;
  context.refinement.fm.maximum_fm_iterations = 100;
  context.refinement.fm.fraction_of_fruitless_moves = 0.05;
  context.refinement.fm.num_seed_nodes = 100;
  context.refinement.fm.min_improvement = 2;
  context.refinement.fm.early_exit_window = 10;

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

  // Multilevel Solver
  start = std::chrono::high_resolution_clock::now();
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
    flat::solve(graph, context);
    // Check if solution is better than best solution found so far
    check_for_new_best_solution();
  };

  // First Iteration
  HighResClockTimepoint s = std::chrono::high_resolution_clock::now();
  solve();
  HighResClockTimepoint e = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> solve_time(e - s);
  double time = solve_time.count();
  if ( time < 240 ) {
    context.general.num_repititions = std::min(std::max(static_cast<int>((240.0 - time) / time), 1), 1000);
  } else {
    context.general.num_repititions = 0;
  }

  for ( int i = 0;
        i < context.general.num_repititions &&
        fruitless_repititions < context.general.num_fruitless_repititions ; ++i ) {
    solve();
  }
  utils::Timer::instance().stop_timer("solver");

  if ( terminate_lock.tryLock() ) {
    graph.applyBestCliques();
    printResult(graph);
  }

  return 0;
}
