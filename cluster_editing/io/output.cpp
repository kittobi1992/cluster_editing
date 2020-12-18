#include "output.h"

#include <vector>

#include "cluster_editing/utils/timer.h"



namespace cluster_editing::io {
  void printContext(const Context& context) {
    if (context.general.verbose_output) {
      LOG << context;
    }
  }

  void printGraphInfo(const Graph& graph, const Context& context, const std::string& name) {
    unused(graph);
    if (context.general.verbose_output) {
      LOG << "Graph:" << context.general.graph_filename << "(" << name << ")";
      INFO("Display Stats of Graph");
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
    LOG << "Objectives:";
    LOG << " Edge Insertions     (minimize) = #total_insertions";
    LOG << " Edge Deletions      (minimize) = #total_deletions";
    LOG << " Total Modifications (minimize) = #total_insertions> + <#total_deletions";
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