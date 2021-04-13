#include "output.h"

#include <algorithm>
#include <vector>
#include <numeric>

#include "cluster_editing/utils/timer.h"
#include "cluster_editing/utils/math.h"
#include "cluster_editing/metrics.h"

namespace cluster_editing {
namespace io {

namespace internal {
  struct Statistic {
    uint64_t min = 0;
    uint64_t q1 = 0;
    uint64_t med = 0;
    uint64_t q3 = 0;
    uint64_t max = 0;
    double avg = 0.0;
    double sd = 0.0;
  };

  template <typename T>
  Statistic createStats(const std::vector<T>& vec, const double avg, const double stdev) {
    internal::Statistic stats;
    if (!vec.empty()) {
      const auto quartiles = math::firstAndThirdQuartile(vec);
      stats.min = vec.front();
      stats.q1 = quartiles.first;
      stats.med = math::median(vec);
      stats.q3 = quartiles.second;
      stats.max = vec.back();
      stats.avg = avg;
      stats.sd = stdev;
    }
    return stats;
  }


  void printStats(const Statistic& node_deg_stats) {
    LOG << "Node Degrees"
        << "| min =" << node_deg_stats.min
        << "| Q1 =" << node_deg_stats.q1
        << "| med =" << node_deg_stats.med
        << "| Q3 =" << node_deg_stats.q3
        << "| max =" << node_deg_stats.max
        << "| avg =" << node_deg_stats.avg
        << "| sd =" << node_deg_stats.sd;
  }
}  // namespace internal

  void printContext(const Context& context) {
    if (context.general.verbose_output) {
      LOG << context;
    }
  }

  void printGraphInfo(const Graph& graph, const Context& context, const std::string& name) {
    unused(graph);
    if (context.general.verbose_output) {
      LOG << "Graph:" << context.general.graph_filename;
      LOG << "Name:" << name << "-"
          << "#Nodes:" << graph.numNodes() << "-"
          << "#Edges:" << graph.numEdges() / 2 << "-"
          << "Total Weight:" << graph.totalWeight();

      std::vector<NodeID> node_degrees(graph.numNodes());
      double avg_node_degree = 0.0;
      double stdev_node_degree = 0.0;
      for ( const NodeID& u : graph.nodes() ) {
        node_degrees[u] = graph.degree(u);
        avg_node_degree += node_degrees[u];
      }
      avg_node_degree /= static_cast<double>(graph.numNodes());

      for ( const NodeID& u : graph.nodes() ) {
        stdev_node_degree += (graph.degree(u) - avg_node_degree) * (graph.degree(u) - avg_node_degree);
      }
      stdev_node_degree = std::sqrt(stdev_node_degree / (graph.numNodes() - 1));

      std::sort(node_degrees.begin(), node_degrees.end());

      internal::printStats(internal::createStats(node_degrees, avg_node_degree, stdev_node_degree));
    }
  }

  void printInputInfo(const Graph& graph, const Context& context) {
    unused(graph);
    if (context.general.verbose_output) {
      LOG << "\n********************************************************************************";
      LOG << "*                                    Input                                     *";
      LOG << "********************************************************************************";
      printGraphInfo(graph, context, "Input Graph");
      printStripe();
    }
  }

  void printPreprocessingBanner(const Context& context) {
    if (context.general.verbose_output) {
      LOG << "\n********************************************************************************";
      LOG << "*                              Preprocessing...                                *";
      LOG << "********************************************************************************";
    }
  }

  void printUndoPreprocessingBanner(const Context& context) {
    if (context.general.verbose_output) {
      LOG << "\n********************************************************************************";
      LOG << "*                            Undo Preprocessing...                             *";
      LOG << "********************************************************************************";
    }
  }

  void printCoarseningBanner(const Context& context) {
    if (context.general.verbose_output) {
      LOG << "\n********************************************************************************";
      LOG << "*                                Coarsening...                                 *";
      LOG << "********************************************************************************";
    }
  }

  void printFlatBanner(const Context& context) {
    if (context.general.verbose_output) {
      LOG << "\n********************************************************************************";
      LOG << "*                             Flat Clustering...                               *";
      LOG << "********************************************************************************";
    }
  }

  void printUncoarseningBanner(const Context& context) {
    if (context.general.verbose_output) {
      LOG << "\n********************************************************************************";
      LOG << "*                              Uncoarsening...                                 *";
      LOG << "********************************************************************************";
    }
  }

  void printObjectives(const Graph& graph,
                       const std::chrono::duration<double>& elapsed_seconds) {
    std::vector<bool> _visited_cliques(graph.numNodes(), false);
    size_t num_cliques = 0;
    for ( const NodeID& u : graph.nodes() ) {
      const CliqueID from = graph.clique(u);
      if ( !_visited_cliques[from] ) {
        ++num_cliques;
        _visited_cliques[from] = true;
      }
    }

    const size_t edge_insertions = metrics::edge_insertions(graph);
    const size_t edge_deletions = metrics::edge_deletions(graph);
    LOG << "Objectives:";
    LOG << " Number of Cliques              =" << num_cliques;
    LOG << " Edge Insertions     (minimize) =" << edge_insertions;
    LOG << " Edge Deletions      (minimize) =" << edge_deletions;
    LOG << " Total Modifications (minimize) =" << (edge_insertions + edge_deletions);
    LOG << " Cluster Editing Time           =" << elapsed_seconds.count() << "s";
  }

  void printClusterEditingResults(const Graph& graph,
                                  const Context& context,
                                  const std::chrono::duration<double>& elapsed_seconds) {
    if (context.general.verbose_output) {
      LOG << "\n********************************************************************************";
      LOG << "*                          Cluster Editing Results                             *";
      LOG << "********************************************************************************";
      printObjectives(graph, elapsed_seconds);

      LOG << "\nTimings:";
      LOG << utils::Timer::instance(true);
    }
  }

  void printResultLine(const Graph& graph,
                       const Context& context,
                       const std::chrono::duration<double>& elapsed_seconds) {
    std::cout << "RESULT"
              << " graph=" << context.general.graph_filename.substr(
                  context.general.graph_filename.find_last_of('/') + 1)
              << " numNodes=" << graph.numNodes()
              << " numEdges=" << graph.numEdges()
              << " seed=" << context.general.seed
              << " totalTime=" << elapsed_seconds.count();

    // Serialize Timings
    utils::Timer::instance().serialize(std::cout);

    // Metrics
    const size_t edge_insertions = metrics::edge_insertions(graph);
    const size_t edge_deletions = metrics::edge_deletions(graph);
    std::cout << " insertions=" << edge_insertions
              << " deletions=" << edge_deletions
              << " modifications=" << (edge_insertions + edge_deletions) << std::endl;
  }

  void printStripe() {
    LOG << "--------------------------------------------------------------------------------";
  }

  void printBanner() {
    LOG << R"(+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++)";
    LOG << R"(+          ___ _           _                __   _ _ _   _                    +)";
    LOG << R"(+         / __\ |_   _ ___| |_ ___ _ __    /__\_| (_) |_(_)_ __   __ _        +)";
    LOG << R"(+        / /  | | | | / __| __/ _ \ '__|  /_\/ _` | | __| | '_ \ / _` |       +)";
    LOG << R"(+       / /___| | |_| \__ \ ||  __/ |    //_| (_| | | |_| | | | | (_| |       +)";
    LOG << R"(+       \____/|_|\__,_|___/\__\___|_|    \__/\__,_|_|\__|_|_| |_|\__, |       +)";
    LOG << R"(+                                                                |___/        +)";
    LOG << R"(+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++)";
  }
}
} // namespace cluster_editing::io