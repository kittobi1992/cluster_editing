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

#include "output.h"

#include <algorithm>
#include <vector>
#include <numeric>
#include <fstream>
#include <cstring>

#include "cluster_editing/datastructures/fast_reset_flag_array.h"
#include "cluster_editing/utils/timer.h"
#include "cluster_editing/utils/math.h"
#include "cluster_editing/utils/common_operations.h"
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

  void printCliqueStats(const Statistic& clique_sizes,
                        const Statistic& intra_edges,
                        const Statistic& inter_edges) {
    LOG << "Clique Sizes"
        << "| min =" << clique_sizes.min
        << "| Q1 =" << clique_sizes.q1
        << "| med =" << clique_sizes.med
        << "| Q3 =" << clique_sizes.q3
        << "| max =" << clique_sizes.max
        << "| avg =" << clique_sizes.avg
        << "| sd =" << clique_sizes.sd;
    LOG << "Intra Clique Edges"
        << "| min =" << intra_edges.min
        << "| Q1 =" << intra_edges.q1
        << "| med =" << intra_edges.med
        << "| Q3 =" << intra_edges.q3
        << "| max =" << intra_edges.max
        << "| avg =" << intra_edges.avg
        << "| sd =" << intra_edges.sd;
    LOG << "Inter Clique Edges"
        << "| min =" << inter_edges.min
        << "| Q1 =" << inter_edges.q1
        << "| med =" << inter_edges.med
        << "| Q3 =" << inter_edges.q3
        << "| max =" << inter_edges.max
        << "| avg =" << inter_edges.avg
        << "| sd =" << inter_edges.sd;
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

  void printCliqueInfo(const Graph& graph, const Context& context) {
    if (context.general.verbose_output) {

      std::vector<NodeID> clique_sizes(graph.numNodes(), 0);
      std::vector<NodeID> intra_edges(graph.numNodes(), 0);
      std::vector<NodeID> inter_edges(graph.numNodes(), 0);
      for ( const NodeID& u : graph.nodes() ) {
        const CliqueID c = graph.clique(u);
        ++clique_sizes[c];
        for ( const NodeID& v : graph.neighbors(u) ) {
          if ( graph.clique(v) == c && u < v ) {
            ++intra_edges[c];
          } else if ( graph.clique(v) != c ) {
            ++inter_edges[c];
          }
        }
      }

      double avg_clique_size = 0.0;
      double avg_intra_edges = 0.0;
      double avg_inter_edges = 0.0;
      double stdev_clique_size = 0.0;
      double stdev_intra_edges = 0.0;
      double stdev_inter_edges = 0.0;
      size_t num_cliques = clique_sizes.size();
      for ( CliqueID c = 0; c < num_cliques; ++c ) {
        if ( clique_sizes[c] == 0 ) {
          std::swap(intra_edges[c], intra_edges[num_cliques - 1]);
          std::swap(inter_edges[c], inter_edges[num_cliques - 1]);
          std::swap(clique_sizes[c--], clique_sizes[--num_cliques]);
          intra_edges.pop_back();
          inter_edges.pop_back();
          clique_sizes.pop_back();
        } else {
          avg_clique_size += clique_sizes[c];
          avg_intra_edges += intra_edges[c];
          avg_inter_edges += intra_edges[c];
        }
      }
      avg_clique_size /= static_cast<double>(num_cliques);
      avg_intra_edges /= static_cast<double>(num_cliques);
      avg_inter_edges /= static_cast<double>(num_cliques);

      for ( CliqueID c = 0; c < clique_sizes.size(); ++c ) {
        stdev_clique_size += (clique_sizes[c] - avg_clique_size) * (clique_sizes[c] - avg_clique_size);
        stdev_intra_edges += (intra_edges[c] - avg_intra_edges) * (intra_edges[c] - avg_intra_edges);
        stdev_inter_edges += (inter_edges[c] - avg_inter_edges) * (inter_edges[c] - avg_inter_edges);
      }
      stdev_clique_size = std::sqrt(stdev_clique_size / (num_cliques - 1));
      stdev_intra_edges = std::sqrt(stdev_intra_edges / (num_cliques - 1));
      stdev_inter_edges = std::sqrt(stdev_inter_edges / (num_cliques - 1));

      std::sort(clique_sizes.begin(), clique_sizes.end());
      std::sort(intra_edges.begin(), intra_edges.end());
      std::sort(inter_edges.begin(), inter_edges.end());

      internal::printCliqueStats(
        internal::createStats(clique_sizes, avg_clique_size, stdev_clique_size),
        internal::createStats(intra_edges, avg_intra_edges, stdev_intra_edges),
        internal::createStats(inter_edges, avg_inter_edges, stdev_inter_edges));
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

  void printClusteringBanner(const Context& context) {
    if (context.general.verbose_output) {
      LOG << "\n********************************************************************************";
      LOG << "*                               Clustering...                                  *";
      LOG << "********************************************************************************";
    }
  }

  void printInitialSolutionBanner(const Context& context) {
    if (context.general.verbose_output) {
      LOG << "\n********************************************************************************";
      LOG << "*                        Create Initial Solution...                            *";
      LOG << "********************************************************************************";
    }
  }

  void printGlobalEvoBanner(const Context& context) {
    if (context.general.verbose_output) {
      LOG << "\n********************************************************************************";
      LOG << "*                       Global Evolutionary Burst...                           *";
      LOG << "********************************************************************************";
    }
  }

  void printLocalizedEvoBanner(const Context& context) {
    if (context.general.verbose_output) {
      LOG << "\n********************************************************************************";
      LOG << "*                      Localized Evolutionary Burst...                         *";
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

  void printEdits(const Graph& graph) {
    utils::CommonOperations::instance(graph).computeNodesOfCliqueWithEmptyCliques(graph);
    const std::vector<std::vector<NodeID>>& cliques =
      utils::CommonOperations::instance(graph)._cliques;
    ds::FastResetFlagArray<> neighbors(graph.numNodes());
    for ( CliqueID c = 0; c < graph.numNodes(); ++c ) {
      if ( cliques[c].size() > 0 ) {
        for ( const NodeID u : cliques[c] ) {
          neighbors.reset();
          neighbors.set(u, true);
          for ( const NodeID& v : graph.neighbors(u) ) {
            if ( graph.clique(u) == graph.clique(v) ) {
              neighbors.set(v, true);
            } else if ( u < v ) {
              // Deletions
              std::cout << (u + 1) << " " << (v + 1) << std::endl;
            }
          }

          // Insertions
          for ( const NodeID v : cliques[c] ) {
            if ( u < v && !neighbors[v] ) {
              std::cout << (u + 1) << " " << (v + 1) << std::endl;
            }
          }
        }
      }
    }
  }

  void printStripe() {
    LOG << "--------------------------------------------------------------------------------";
  }

  void printBanner() {
    LOG << R"(+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++)";
    LOG << R"(+             ██ ▄█▀▄▄▄       ██▓███   ▒█████   ▄████▄  ▓█████                +)";
    LOG << R"(+             ██▄█▒▒████▄    ▓██░  ██▒▒██▒  ██▒▒██▀ ▀█  ▓█   ▀                +)";
    LOG << R"(+            ▓███▄░▒██  ▀█▄  ▓██░ ██▓▒▒██░  ██▒▒▓█    ▄ ▒███                  +)";
    LOG << R"(+            ▓██ █▄░██▄▄▄▄██ ▒██▄█▓▒ ▒▒██   ██░▒▓▓▄ ▄██▒▒▓█  ▄                +)";
    LOG << R"(+            ▒██▒ █▄▓█   ▓██▒▒██▒ ░  ░░ ████▓▒░▒ ▓███▀ ░░▒████▒               +)";
    LOG << R"(+            ▒ ▒▒ ▓▒▒▒   ▓▒█░▒▓▒░ ░  ░░ ▒░▒░▒░ ░ ░▒ ▒  ░░░ ▒░ ░               +)";
    LOG << R"(+            ░ ░▒ ▒░ ▒   ▒▒ ░░▒ ░       ░ ▒ ▒░   ░  ▒    ░ ░  ░               +)";
    LOG << R"(+            ░ ░░ ░  ░   ▒   ░░       ░ ░ ░ ▒  ░           ░                  +)";
    LOG << R"(+            ░  ░        ░  ░             ░ ░  ░ ░         ░  ░               +)";
    LOG << R"(+                                              ░                              +)";
    LOG << R"(+            The Karlsruhe and Potsdam Cluster Editing Framework              +)";
    LOG << R"(+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++)";
  }

  void readSolutionFile(Graph& graph, const std::string& filename) {
  ASSERT(!filename.empty(), "No filename for solution file specified");
  std::ifstream file(filename);
  if (file) {
    CliqueID clique = INVALID_CLIQUE;
    NodeID u = 0;
    while (file >> clique) {
      graph.setClique(u, clique);
      ++u;
    }
    file.close();
  }
}

void writeSolutionFile(const Graph& graph, const std::string& filename) {
  if (filename.empty()) {
    LOG << "No filename for partition file specified";
  } else {
    std::ofstream out_stream(filename.c_str());
    for ( const NodeID& u : graph.nodes() ) {
      out_stream << graph.clique(u) << std::endl;
    }
    out_stream.close();
  }
}
} // namespace cluster_editing::io