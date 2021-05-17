#include "exact_refiner.h"

#include <queue>

#include "cluster_editing/exact/instance.h"
#include "cluster_editing/exact/solver.h"
#include "cluster_editing/utils/common_operations.h"
#include "cluster_editing/utils/timer.h"

namespace cluster_editing {

void ExactRefiner::initializeImpl(Graph&) {

}

EdgeWeight ExactRefiner::refineImpl(Graph& graph,
                                    const EdgeWeight current_edits,
                                    const EdgeWeight) {
  EdgeWeight current_metric = current_edits;

  utils::Timer::instance().start_timer("exact", "ExactRefiner");
  _visited.reset();
  std::random_shuffle(_nodes.begin(), _nodes.end());
  utils::CommonOperations::instance(graph).computeNodesOfCliqueWithEmptyCliques(graph);
  std::vector<CliqueID> empty_cliques = utils::CommonOperations::instance(graph)._empty_cliques;
  std::vector<std::vector<NodeID>> cliques = utils::CommonOperations::instance(graph)._cliques;

  // Flag large cliques as already visited
  for ( const CliqueID c : graph.nodes() ) {
    if ( cliques[c].size() >= _context.refinement.exact.max_subgraph_size ) {
      for ( const NodeID u : cliques[c] ) {
        _visited.set(u, true);
      }
    }
  }

  for ( const NodeID& u : _nodes ) {
    if ( !_visited[u] ) {
      std::vector<NodeID> subgraph = growSubgraph(graph, u, cliques);
      if ( subgraph.size() > _context.refinement.exact.min_subgraph_size ) {
        _mapping.clear();
        for ( NodeID i = 0; i < subgraph.size(); ++i ) {
          _mapping[subgraph[i]] = i;
        }

        // Create instance for exact solver
        size_t insertions = 0;
        size_t deletions = 0;
        _in_queue.reset();
        for ( const NodeID& u : subgraph ) {
          const CliqueID from = graph.clique(u);
          if ( !_in_queue[from] ) {
            insertions += ((cliques[from].size() - 1) * cliques[from].size()) / 2;
            _in_queue.set(from, true);
          }
        }

        Instance instance(subgraph.size());
        for ( const NodeID& u : subgraph ) {
          const CliqueID u_c = graph.clique(u);
          const NodeID mapped_u = _mapping[u];
          for ( const NodeID& v : graph.neighbors(u) ) {
            if ( _mapping.contains(v) ) {
              const CliqueID v_c = graph.clique(v);
              const NodeID mapped_v = _mapping[v];
              instance.edges[mapped_u][mapped_v] = 1;
              instance.orig[mapped_u][mapped_v] = 1;
              if ( u_c == v_c && u < v ) {
                --insertions;
              } else if ( u_c < v_c ) {
                ++deletions;
              }
            }
          }
        }

        // Solve instance
        const size_t subgraph_edits = insertions + deletions;
        ExactSolver solver;
        solver.verbose = false;
        solver.time_limit = std::chrono::steady_clock::now() +
          std::chrono::milliseconds(_context.refinement.exact.time_limit);
        Solution solution = solver.solve(instance);
        if ( subgraph_edits > solution.cost ) {
          LOG << GREEN << V(subgraph_edits) << V(solution.cost) << END;
        } else {
          // LOG << RED << V(subgraph_edits) << V(solution.cost) << END;
        }
      }
    }
  }

  utils::Timer::instance().stop_timer("exact");

  return current_metric;
}

std::vector<NodeID> ExactRefiner::growSubgraph(const Graph& graph,
                                               const NodeID seed,
                                               const std::vector<std::vector<NodeID>>& cliques) {
  ASSERT(!_visited[seed]);
  _in_queue.reset();
  std::vector<NodeID> subgraph;
  std::queue<NodeID> q;
  q.push(seed);
  _in_queue.set(seed, true);

  while ( !q.empty() && subgraph.size() < _context.refinement.exact.max_subgraph_size ) {
    const NodeID u = q.front();
    q.pop();

    const CliqueID from = graph.clique(u);
    if ( !_visited[u] &&
         cliques[from].size() + subgraph.size() <= _context.refinement.exact.max_subgraph_size ) {
      for ( const NodeID& v : cliques[from] ) {
        _visited.set(v, true);
        _in_queue.set(v, true);
        subgraph.push_back(v);
        for ( const NodeID& w : graph.neighbors(v) ) {
          if ( !_in_queue[w] && !_visited[w] ) {
            q.push(w);
            _in_queue.set(w, true);
          }
        }
      }
    }
  }

  return subgraph;
}


} // namespace cluster_editing