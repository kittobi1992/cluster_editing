#include "lp_coarsener.h"

#include "cluster_editing/io/output.h"
#include "cluster_editing/utils/randomize.h"
#include "cluster_editing/utils/timer.h"
#include "cluster_editing/metrics.h"

namespace cluster_editing {

void LabelPropagationCoarsener::coarsenImpl() {

  while ( true ) {
    Graph& current_graph = currentGraph();
    _lp_refiner.initialize(current_graph);

    utils::Timer::instance().start_timer("clustering", "Clustering");
    _lp_refiner.refine(current_graph);
    utils::Timer::instance().stop_timer("clustering");

    if ( _lp_refiner.movedVertices() > 0 ) {
      utils::Timer::instance().start_timer("contraction", "Contraction");
      _hierarchies.emplace_back(current_graph.contract());
      utils::Timer::instance().stop_timer("contraction");

      if ( debug ) {
        io::printGraphInfo(_hierarchies.back().first, _context, "Contracted Graph");
      }

    } else {
      break;
    }
  }

  if ( _context.general.verbose_output ) {
    const EdgeWeight edge_insertions = metrics::edge_insertions(currentGraph());
    const EdgeWeight edge_deletions = metrics::edge_deletions(currentGraph());
    LOG << "Initial Solution:";
    LOG << "Number of Hierarchies:" << _hierarchies.size();
    LOG << "Edge Insertions:" << edge_insertions << "-"
        << "Edge Deletions:" << edge_deletions << "-"
        << "Edge Modifications:" << (edge_deletions + edge_insertions);
    io::printGraphInfo(_hierarchies.back().first, _context, "Contracted Graph");
  }
}

namespace {

  void refine(Graph& graph, std::unique_ptr<IRefiner>& refiner) {
    utils::Timer::instance().start_timer("initialize_refiner", "Initialize Refiner");
    refiner->initialize(graph);
    utils::Timer::instance().stop_timer("initialize_refiner");

    utils::Timer::instance().start_timer("refine", "Refine");
    refiner->refine(graph);
    utils::Timer::instance().stop_timer("refine");
  }

  void project(Graph& original_graph,
               const Graph& coarse_graph,
               const std::vector<NodeID>& clique_to_coarse) {
    utils::Timer::instance().start_timer("project_clique", "Projecting Cliques");
    for ( const NodeID& u : original_graph.nodes() ) {
      const CliqueID& u_c = original_graph.clique(u);
      original_graph.setClique(u, coarse_graph.clique(clique_to_coarse[u_c]));
    }
    utils::Timer::instance().stop_timer("project_clique");
  }

} // namespace

void LabelPropagationCoarsener::uncoarsenImpl(std::unique_ptr<IRefiner>& refiner) {
  // Refine coarsest graph
  refine(currentGraph(), refiner);

  while ( true ) {
    const Graph& coarse_graph = currentGraph();
    if ( &coarse_graph != &_graph ) {
      // Project partition to next level finer graph
      const size_t h_size = _hierarchies.size();
      if ( h_size > 1 ) {
        project(_hierarchies[h_size - 2].first, coarse_graph, _hierarchies.back().second);
      } else {
        project(_graph, coarse_graph, _hierarchies.back().second);
      }
      _hierarchies.pop_back();

      // Refine current graph
      refine(currentGraph(), refiner);
    } else {
      break;
    }
  }
}

} // namespace cluster_editing