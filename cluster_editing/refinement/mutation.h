/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include "cluster_editing/context/context.h"
#include "cluster_editing/definitions.h"
#include "cluster_editing/utils/randomize.h"

namespace cluster_editing {

struct Clique {
  CliqueID c;
  std::vector<NodeID> nodes;

  size_t size() {
    return nodes.size();
  }
};

struct CliqueStats {
  explicit CliqueStats(const Graph& graph) :
    cliques(graph.numNodes()),
    empty_cliques(),
    num_cliques(0) {
    for ( const NodeID& u : graph.nodes() ) {
      cliques[u].c = u;
      cliques[graph.clique(u)].nodes.push_back(u);
    }
    for ( const NodeID& u : graph.nodes() ) {
      if ( cliques[u].nodes.size() == 0 ) {
        empty_cliques.push_back(u);
      } else {
        ++num_cliques;
      }
    }
  }

  std::vector<Clique> cliques;
  std::vector<CliqueID> empty_cliques;
  size_t num_cliques;
};

class LargeCliqueIsolator {

 public:
  static void mutate(Graph& graph, const Context& context) {
    CliqueStats stats(graph);
    for ( const CliqueID& c : graph.nodes() ) {
      if ( stats.cliques[c].size() >= context.refinement.evo.large_clique_threshold ) {
        const float p = utils::Randomize::instance().getRandomFloat(0.0, 1.0);
        if ( p <= context.refinement.evo.clique_isolate_prob ) {
          stats.empty_cliques.push_back(c);
          for ( const NodeID& u : stats.cliques[c].nodes ) {
            ASSERT(!stats.empty_cliques.empty());
            const CliqueID target = stats.empty_cliques.back();
            stats.empty_cliques.pop_back();
            graph.setClique(u, target);
          }
        }
      }
    }
  }

 private:
  LargeCliqueIsolator() { }
};

class LargeCliqueWithNeighborIsolator {

 public:
  static void mutate(Graph& graph, const Context& context) {
    CliqueStats stats(graph);
    for ( const CliqueID& c : graph.nodes() ) {
      if ( stats.cliques[c].size() >= context.refinement.evo.large_clique_threshold ) {
        const float p = utils::Randomize::instance().getRandomFloat(0.0, 1.0);
        if ( p <= context.refinement.evo.neighbor_clique_isolate_prob ) {
          stats.empty_cliques.push_back(c);
          std::vector<NodeID> nodes_to_isolate;
          for ( const NodeID& u : stats.cliques[c].nodes ) {
            nodes_to_isolate.push_back(u);
            for ( const Neighbor& n : graph.neighbors(u) ) {
              const CliqueID target = graph.clique(n.target);
              if ( target != c && stats.cliques[target].size() > 0 ) {
                stats.empty_cliques.push_back(target);
                for ( const NodeID& v : stats.cliques[target].nodes ) {
                  nodes_to_isolate.push_back(v);
                }
                stats.cliques[target].nodes.clear();
              }
            }
          }
          stats.cliques[c].nodes.clear();

          // Isolate all nodes and clique c and all nodes in neighbor cliques
          for ( const NodeID& u : nodes_to_isolate ) {
            ASSERT(!stats.empty_cliques.empty());
            const CliqueID target = stats.empty_cliques.back();
            stats.empty_cliques.pop_back();
            graph.setClique(u, target);
          }
        }
      }
    }
  }

 private:
  LargeCliqueWithNeighborIsolator() { }
};

}  // namespace cluster_editing
