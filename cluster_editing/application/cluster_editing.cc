#include <iostream>
#include <set>

#include "cluster_editing/definitions.h"
#include "cluster_editing/preprocessing.h"
#include "cluster_editing/flat.h"
#include "cluster_editing/multilevel.h"
#include "cluster_editing/utils/timer.h"
#include "cluster_editing/utils/randomize.h"
#include "cluster_editing/metrics.h"
#include "cluster_editing/io/output.h"
#include "cluster_editing/io/graph_io.h"

using namespace cluster_editing;

int main() {
  Context context;
  context.general.verbose_output = false;
  context.general.print_result_line = true;
  context.general.seed = 0;
  context.general.use_multilevel = false;
  context.general.num_repititions = 1;
  context.refinement.use_lp_refiner = true;
  context.refinement.lp.maximum_lp_iterations = 100;
  context.refinement.lp.activate_all_cliques_after_rounds = 10;
  context.refinement.lp.random_shuffle_each_round = false;
  context.refinement.lp.node_order = NodeOrdering::none;
  utils::Randomize::instance().setSeed(context.general.seed);

  if ( context.general.verbose_output ) {
    io::printBanner();
    // Print context description
    LOG << context;
  }

  utils::Timer::instance().start_timer("import_graph", "Import Graph");
  Graph graph = io::readGraphFile();
  utils::Timer::instance().stop_timer("import_graph");
  io::printInputInfo(graph, context);

  // Preprocessing
  io::printPreprocessingBanner(context);
  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  utils::Timer::instance().start_timer("preprocessing", "Preprocessing");
  Preprocessor preprocessor(graph, context);
  preprocessor.preprocess();
  utils::Timer::instance().stop_timer("preprocessing");

  // Multilevel Solver
  utils::Timer::instance().start_timer("solver", "Solver");
  std::vector<CliqueID> best_cliques(graph.numNodes(), INVALID_CLIQUE);
  size_t best_objective = std::numeric_limits<size_t>::max();
  for ( int i = 0; i < context.general.num_repititions; ++i ) {
    graph.reset();
    if ( context.general.use_multilevel ) {
      multilevel::solve(graph, context);
    } else {
      flat::solve(graph, context);
    }

    // Check if solution is better than best solution found so far
    const size_t current_objective = metrics::edits(graph);
    if ( current_objective < best_objective ) {
      if ( context.general.verbose_output ) {
        LOG << GREEN << "Improved best solution from"
            << best_objective << "to" << current_objective << END;
      }
      for ( const NodeID& u : graph.nodes() ) {
        best_cliques[u] = graph.clique(u);
      }
      best_objective = current_objective;
    }
  }

  // Apply best found solution
  for ( const NodeID& u : graph.nodes() ) {
    graph.setClique(u, best_cliques[u]);
  }

  utils::Timer::instance().stop_timer("solver");

  // Undo Preprocessing
  io::printUndoPreprocessingBanner(context);
  utils::Timer::instance().start_timer("undo_preprocessing", "Undo Preprocessing");
  preprocessor.undoPreprocessing();
  utils::Timer::instance().stop_timer("undo_preprocessing");
  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();

  // Print RESULT line
  if ( context.general.print_result_line ) {
    // Deletions
    std::vector<std::vector<NodeID>> cliques(graph.numNodes());
    for ( const NodeID u : graph.nodes() ) {
      cliques[graph.clique(u)].push_back(u);
      for ( const auto& v : graph.neighbors(u) ) {
        if ( u < v.target && graph.clique(u) != graph.clique(v.target) ) {
          std::cout << (u + 1) << " " << (v.target + 1) << std::endl;
        }
      }
    }

    // Insertions
    for ( CliqueID c = 0; c < graph.numNodes(); ++c ) {
      if ( cliques[c].size() > 0 ) {
        std::set<NodeID> clique(cliques[c].begin(), cliques[c].end());

        for ( const NodeID u : cliques[c] ) {
          for ( const auto& v : graph.neighbors(u) ) {
            if ( graph.clique(u) == graph.clique(v.target) ) {
              clique.erase(v.target);
            }
          }

          for ( const NodeID v : clique ) {
            if ( u < v ) {
              std::cout << (u + 1) << " " << (v + 1) << std::endl;
            }
          }

          for ( const auto& v : graph.neighbors(u) ) {
            if ( graph.clique(u) == graph.clique(v.target) ) {
              clique.insert(v.target);
            }
          }
        }
      }
    }
  }

  // Print Stats
  std::chrono::duration<double> elapsed_seconds(end - start);
  io::printClusterEditingResults(graph, context, elapsed_seconds);

  if ( context.general.print_csv ) {
    std::cout << context.general.graph_filename.substr(context.general.graph_filename.find_last_of('/') + 1)
              << "," << metrics::edits(graph) << "," << elapsed_seconds.count() << std::endl;
  }

  return 0;
}
