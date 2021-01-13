#include "output.h"

#include <algorithm>
#include <vector>
#include <numeric>

#include "cluster_editing/utils/timer.h"
#include "cluster_editing/utils/math.h"
#include "cluster_editing/metrics.h"

namespace cluster_editing::io {

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


  void printStats(const Statistic& node_deg_stats,
                  const Statistic& node_weight_stats,
                  const Statistic& selfloop_weight_stats,
                  const Statistic& edge_weight_stats) {
    // default double precision is 7
    const uint8_t double_width = 7;
    const uint8_t node_deg_width = std::max(math::digits(node_deg_stats.max), double_width) + 4;
    const uint8_t node_weight_width = std::max(math::digits(node_weight_stats.max), double_width) + 4;
    const uint8_t selfloop_weight_width = std::max(math::digits(selfloop_weight_stats.max), double_width) + 4;
    const uint8_t edge_weight_width = std::max(math::digits(edge_weight_stats.max), double_width) + 4;

    LOG << "Node Degree" << std::right << std::setw(node_deg_width + 8)
        << "Node Weight" << std::right << std::setw(node_weight_width + 12)
        << "Selfloop Weight" << std::right << std::setw(selfloop_weight_width + 4)
        << "Edge Weight" << std::right << std::setw(edge_weight_width + 8);
    LOG << "| min=" << std::left << std::setw(node_deg_width) << node_deg_stats.min
        << " | min=" << std::left << std::setw(node_weight_width) << node_weight_stats.min
        << " | min=" << std::left << std::setw(node_weight_width) << selfloop_weight_stats.min
        << " | min=" << std::left << std::setw(edge_weight_width) << edge_weight_stats.min;
    LOG << "| Q1 =" << std::left << std::setw(node_deg_width) << node_deg_stats.q1
        << " | Q1 =" << std::left << std::setw(node_weight_width) << node_weight_stats.q1
        << " | Q1 =" << std::left << std::setw(node_weight_width) << selfloop_weight_stats.q1
        << " | Q1 =" << std::left << std::setw(edge_weight_width) << edge_weight_stats.q1;
    LOG << "| med=" << std::left << std::setw(node_deg_width) << node_deg_stats.med
        << " | med=" << std::left << std::setw(node_weight_width) << node_weight_stats.med
        << " | med=" << std::left << std::setw(node_weight_width) << selfloop_weight_stats.med
        << " | med=" << std::left << std::setw(edge_weight_width) << edge_weight_stats.med;
    LOG << "| Q3 =" << std::left << std::setw(node_deg_width) << node_deg_stats.q3
        << " | Q3 =" << std::left << std::setw(node_weight_width) << node_weight_stats.q3
        << " | Q3 =" << std::left << std::setw(node_weight_width) << selfloop_weight_stats.q3
        << " | Q3 =" << std::left << std::setw(edge_weight_width) << edge_weight_stats.q3;
    LOG << "| max=" << std::left << std::setw(node_deg_width) << node_deg_stats.max
        << " | max=" << std::left << std::setw(node_weight_width) << node_weight_stats.max
        << " | max=" << std::left << std::setw(node_weight_width) << selfloop_weight_stats.max
        << " | max=" << std::left << std::setw(edge_weight_width) << edge_weight_stats.max;
    LOG << "| avg=" << std::left << std::setw(node_deg_width) << node_deg_stats.avg
        << " | avg=" << std::left << std::setw(node_weight_width) << node_weight_stats.avg
        << " | avg=" << std::left << std::setw(node_weight_width) << selfloop_weight_stats.avg
        << " | avg=" << std::left << std::setw(edge_weight_width) << edge_weight_stats.avg;
    LOG << "| sd =" << std::left << std::setw(node_deg_width) << node_deg_stats.sd
        << " | sd =" << std::left << std::setw(node_weight_width) << node_weight_stats.sd
        << " | sd =" << std::left << std::setw(node_weight_width) << selfloop_weight_stats.sd
        << " | sd =" << std::left << std::setw(edge_weight_width) << edge_weight_stats.sd;
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
      std::vector<NodeWeight> node_weights(graph.numNodes());
      std::vector<NodeWeight> selfloop_weights(graph.numNodes());
      std::vector<EdgeWeight> edge_weights(graph.numEdges());

      double avg_node_degree = 0.0;
      double stdev_node_degree = 0.0;
      double avg_node_weight = 0.0;
      double stdev_node_weight = 0.0;
      double avg_selfloop_weight = 0.0;
      double stdev_selfloop_weight = 0.0;
      for ( const NodeID& u : graph.nodes() ) {
        node_degrees[u] = graph.degree(u);
        avg_node_degree += node_degrees[u];
        node_weights[u] = graph.nodeWeight(u);
        avg_node_weight += node_weights[u];
        selfloop_weights[u] = graph.selfloopWeight(u);
        avg_selfloop_weight += selfloop_weights[u];
      }
      avg_node_degree /= static_cast<double>(graph.numNodes());
      avg_node_weight /= static_cast<double>(graph.numNodes());
      avg_selfloop_weight /= static_cast<double>(graph.numNodes());

      for ( const NodeID& u : graph.nodes() ) {
        stdev_node_degree += (graph.degree(u) - avg_node_degree) * (graph.degree(u) - avg_node_degree);
        stdev_node_weight += (graph.nodeWeight(u) - avg_node_weight) * (graph.nodeWeight(u) - avg_node_weight);
        stdev_selfloop_weight += (graph.selfloopWeight(u) - avg_selfloop_weight) * (graph.selfloopWeight(u) - avg_selfloop_weight);
      }
      stdev_node_degree = std::sqrt(stdev_node_degree / (graph.numNodes() - 1));
      stdev_node_weight = std::sqrt(stdev_node_weight / (graph.numNodes() - 1));
      stdev_selfloop_weight = std::sqrt(stdev_selfloop_weight / (graph.numNodes() - 1));

      double avg_edge_weight = 0.0;
      double stdev_edge_weight = 0.0;
      for ( const EdgeID& e : graph.edges() ) {
        edge_weights[e] = graph.edgeWeight(e);
        avg_edge_weight += edge_weights[e];
      }
      avg_edge_weight /= static_cast<double>(graph.numEdges());

      for ( const EdgeID& e : graph.edges() ) {
        stdev_edge_weight += (graph.edgeWeight(e) - avg_edge_weight) * (graph.edgeWeight(e) - avg_edge_weight);
      }
      stdev_edge_weight = std::sqrt(stdev_edge_weight / (graph.numEdges() - 1));

      std::sort(node_degrees.begin(), node_degrees.end());
      std::sort(node_weights.begin(), node_weights.end());
      std::sort(selfloop_weights.begin(), selfloop_weights.end());
      std::sort(edge_weights.begin(), edge_weights.end());

      internal::printStats(internal::createStats(node_degrees, avg_node_degree, stdev_node_degree),
                          internal::createStats(node_weights, avg_node_weight, stdev_node_weight),
                          internal::createStats(selfloop_weights, avg_selfloop_weight, stdev_selfloop_weight),
                          internal::createStats(edge_weights, avg_edge_weight, stdev_edge_weight));
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
      LOG << "********************************************************************************";
      LOG << "*                                Coarsening...                                 *";
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
    unused(graph);
    const size_t edge_insertions = metrics::edge_insertions(graph);
    const size_t edge_deletions = metrics::edge_deletions(graph);
    LOG << "Objectives:";
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
} // namespace cluster_editing::io