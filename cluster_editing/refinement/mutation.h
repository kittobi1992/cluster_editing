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
#include "cluster_editing/utils/common_operations.h"
#include "cluster_editing/refinement/action_selector.h"

namespace cluster_editing {

class LargeCliqueIsolator {

 public:
  static void mutate(Graph& graph, const Context& context, const bool prob) {
    utils::CommonOperations::instance(graph).computeNodesOfCliqueWithEmptyCliques(graph);
    std::vector<CliqueID>& empty_cliques =
      utils::CommonOperations::instance(graph)._empty_cliques;
    std::vector<std::vector<NodeID>> cliques =
      utils::CommonOperations::instance(graph)._cliques;
    for ( const CliqueID& c : graph.nodes() ) {
      if ( cliques[c].size() >= context.refinement.evo.large_clique_threshold ) {
        const float p = utils::Randomize::instance().getRandomFloat(0.0, 1.0);
        if ( p <= prob ) {
          empty_cliques.push_back(c);
          for ( const NodeID& u : cliques[c] ) {
            ASSERT(!empty_cliques.empty());
            const CliqueID target = empty_cliques.back();
            empty_cliques.pop_back();
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
  static void mutate(Graph& graph, const Context& context, const float prob) {
    utils::CommonOperations::instance(graph).computeNodesOfCliqueWithEmptyCliques(graph);
    std::vector<CliqueID>& empty_cliques =
      utils::CommonOperations::instance(graph)._empty_cliques;
    std::vector<std::vector<NodeID>> cliques =
      utils::CommonOperations::instance(graph)._cliques;
    for ( const CliqueID& c : graph.nodes() ) {
      if ( cliques[c].size() >= context.refinement.evo.large_clique_threshold ) {
        const float p = utils::Randomize::instance().getRandomFloat(0.0, 1.0);
        if ( p <= prob ) {
          empty_cliques.push_back(c);
          std::vector<NodeID> nodes_to_isolate;
          for ( const NodeID& u : cliques[c] ) {
            nodes_to_isolate.push_back(u);
            for ( const NodeID& v : graph.neighbors(u) ) {
              const CliqueID target = graph.clique(v);
              if ( target != c && cliques[target].size() > 0 ) {
                empty_cliques.push_back(target);
                for ( const NodeID& v : cliques[target] ) {
                  nodes_to_isolate.push_back(v);
                }
                cliques[target].clear();
              }
            }
          }
          cliques[c].clear();

          // Isolate all nodes and clique c and all nodes in neighbor cliques
          for ( const NodeID& u : nodes_to_isolate ) {
            ASSERT(!empty_cliques.empty());
            const CliqueID target = empty_cliques.back();
            empty_cliques.pop_back();
            graph.setClique(u, target);
          }
        }
      }
    }
  }

 private:
  LargeCliqueWithNeighborIsolator() { }
};

class RandomNodeIsolator {

 public:
  static void mutate(Graph& graph, const float prob) {
    utils::CommonOperations::instance(graph).computeNodesOfCliqueWithEmptyCliques(graph);
    std::vector<CliqueID>& empty_cliques =
      utils::CommonOperations::instance(graph)._empty_cliques;
    std::vector<std::vector<NodeID>> cliques =
      utils::CommonOperations::instance(graph)._cliques;
    for ( const CliqueID& c : graph.nodes() ) {
      if ( cliques[c].size() > 1 ) {
        size_t current_size = cliques[c].size();
        for ( const NodeID& u : cliques[c] ) {
          const float p = utils::Randomize::instance().getRandomFloat(0.0, 1.0);
          if ( p <= prob ) {
            const CliqueID target = empty_cliques.back();
            empty_cliques.pop_back();
            graph.setClique(u, target);
            --current_size;
          }
          if ( current_size == 1 ) {
            break;
          }
        }
      }
    }
  }

 private:
  RandomNodeIsolator() { }
};

class RandomNodeMover {

 public:
  static void mutate(Graph& graph, const float prob) {
    std::vector<CliqueID> target_cliques;
    for ( const NodeID& u : graph.nodes() ) {
      const float p = utils::Randomize::instance().getRandomFloat(0.0, 1.0);
      if ( p <= prob ) {
        const CliqueID from = graph.clique(u);
        for ( const NodeID& v : graph.neighbors(u) ) {
          const CliqueID to = graph.clique(v);
          if ( from != to ) {
            target_cliques.push_back(to);
          }
        }

        if ( !target_cliques.empty() ) {
          std::sort(target_cliques.begin(), target_cliques.end());
          target_cliques.erase( std::unique(
            target_cliques.begin(), target_cliques.end() ), target_cliques.end() );
          const CliqueID target = target_cliques[
            utils::Randomize::instance().getRandomInt(0, target_cliques.size() - 1)];
          graph.setClique(u, target);
          target_cliques.clear();
        }
      }
    }
  }

 private:
  RandomNodeMover() { }
};

class CliqueSplit {

 public:
  static void mutate(Graph& graph, const float prob) {
    utils::CommonOperations::instance(graph).computeNodesOfCliqueWithEmptyCliques(graph);
    utils::CommonOperations::instance(graph).computeClusterSizes(graph);
    std::vector<CliqueID>& empty_cliques =
      utils::CommonOperations::instance(graph)._empty_cliques;
    std::vector<std::vector<NodeID>>& cliques =
      utils::CommonOperations::instance(graph)._cliques;
    std::vector<NodeID> cluster_sizes =
      utils::CommonOperations::instance(graph)._cluster_sizes;
    for ( const CliqueID& c : graph.nodes() ) {
      if ( cliques[c].size() > 2 && utils::Randomize::instance().getRandomFloat(0.0, 1.0) <= prob ) {
        // Choose vertex from clique that we isolate
        const int seed_idx = utils::Randomize::instance().getRandomInt(0, cliques[c].size() - 1);
        const NodeID seed = cliques[c][seed_idx];
        std::swap(cliques[c][seed_idx], cliques[c].back());
        cliques[c].pop_back();

        // Isolate selected vertex
        const CliqueID to = empty_cliques.back();
        empty_cliques.pop_back();
        graph.setClique(seed, to);
        --cluster_sizes[c];
        ++cluster_sizes[to];

        // Choose one additional vertex to be moved to isolated vertex
        std::vector<NodeID> best_nodes;
        EdgeWeight best_delta = std::numeric_limits<EdgeWeight>::max();
        for ( const NodeID& u : cliques[c] ) {
          bool is_adjacent_to_seed = false;
          EdgeWeight edges_to_from = 0;
          for ( const NodeID& v : graph.neighbors(u) ) {
            if ( c == graph.clique(v) ) {
              ++edges_to_from;
            } else if ( v == seed ) {
              is_adjacent_to_seed = true;
            }
          }

          const NodeID u_degree = graph.degree(u);
          const EdgeWeight from_rating = cluster_sizes[c] - 1 + u_degree - 2 * edges_to_from;
          const EdgeWeight to_rating = 1 + u_degree - 2 * is_adjacent_to_seed;
          const EdgeWeight delta = to_rating - from_rating;
          if ( delta < best_delta ) {
            best_nodes.clear();
            best_nodes.push_back(u);
            best_delta = delta;
          } else if ( delta == best_delta ) {
            best_nodes.push_back(u);
          }
        }

        if ( !best_nodes.empty() ) {
          const NodeID best = best_nodes[
            utils::Randomize::instance().getRandomInt(0, best_nodes.size() - 1)];
          graph.setClique(best, to);
          --cluster_sizes[c];
          ++cluster_sizes[to];
        }
      }
    }
  }

 private:
  CliqueSplit() { }
};

class TestMutation {

 public:
  static void mutate(Graph&, const float) {

  }

 private:
  TestMutation() { }
};

class Mutator {
 public:
  explicit Mutator(const Context& context) :
    _context(context),
    _show_detailed_output(context.general.verbose_output && context.refinement.evo.enable_detailed_output),
    _mutation_selector() {
    activateMutations(_context.refinement.evo.enabled_mutations);
  }

  Mutation mutate(Graph& graph) {
    const Mutation mutation = _mutation_selector.chooseAction(_show_detailed_output);
    const float prob = chooseMutationProbability(mutation);
    if ( mutation == Mutation::LARGE_CLIQUE_ISOLATOR ) {
      if ( _show_detailed_output )
        LOG << "Mutation Action: LARGE_CLIQUE_ISOLATOR ( p =" << prob << ")";
      LargeCliqueIsolator::mutate(graph, _context, prob);
    } else if ( mutation == Mutation::LARGE_CLIQUE_WITH_NEIGHBOR_ISOLATOR ) {
      if ( _show_detailed_output )
        LOG << "Mutation Action: LARGE_CLIQUE_WITH_NEIGHBOR_ISOLATOR ( p =" << prob << ")";
      LargeCliqueWithNeighborIsolator::mutate(graph, _context, prob);
    } else if ( mutation == Mutation::RANDOM_NODE_ISOLATOR ) {
      if ( _show_detailed_output )
        LOG << "Mutation Action: RANDOM_NODE_ISOLATOR ( p =" << prob << ")";
      RandomNodeIsolator::mutate(graph, prob);
    } else if ( mutation == Mutation::RANDOM_NODE_MOVER ) {
      if ( _show_detailed_output )
        LOG << "Mutation Action: RANDOM_NODE_MOVER ( p =" << prob << ")";
      RandomNodeMover::mutate(graph, prob);
    } else if ( mutation == Mutation::CLIQUE_SPLITTER ) {
      if ( _show_detailed_output )
        LOG << "Mutation Action: CLIQUE_SPLITTER ( p =" << prob << ")";
      CliqueSplit::mutate(graph, prob);
    } else if ( mutation == Mutation::TEST_MUTATION ) {
      if ( _show_detailed_output )
        LOG << "Mutation Action: TEST_MUTATION ( p =" << prob << ")";
      TestMutation::mutate(graph, prob);
    }
    return mutation;
  }

  void notifyResult(const Mutation& mutation,
                    const EdgeWeight before_edits,
                    const EdgeWeight after_edits) {
    _mutation_selector.notifyImprovement(
      mutation, std::max( 0, before_edits - after_edits ));
  }

  void activateMutations(const std::string& enabled_mutations) {
    std::vector<Mutation> active_mutations;
    for ( size_t i = 0; i < enabled_mutations.size(); ++i ) {
      if ( enabled_mutations[i] == '1' ) {
        active_mutations.push_back(static_cast<Mutation>(i));
      }
    }
    _mutation_selector.initializeActions(active_mutations);
  }

 private:
  float chooseMutationProbability(const Mutation& mutation) {
    if ( mutation == Mutation::LARGE_CLIQUE_ISOLATOR ) {
      return utils::Randomize::instance().getRandomFloat(
        _context.refinement.evo.min_clique_isolate_prob,
        _context.refinement.evo.max_clique_isolate_prob);
    } else if ( mutation == Mutation::LARGE_CLIQUE_WITH_NEIGHBOR_ISOLATOR ) {
      return utils::Randomize::instance().getRandomFloat(
        _context.refinement.evo.min_neighbor_clique_isolate_prob,
        _context.refinement.evo.max_neighbor_clique_isolate_prob);
    } else if ( mutation == Mutation::RANDOM_NODE_ISOLATOR ) {
      return utils::Randomize::instance().getRandomFloat(
        _context.refinement.evo.min_node_isolation_prob,
        _context.refinement.evo.max_node_isolation_prob);
    } else if ( mutation == Mutation::RANDOM_NODE_MOVER ) {
      return utils::Randomize::instance().getRandomFloat(
        _context.refinement.evo.min_node_move_prob,
        _context.refinement.evo.max_node_move_prob);
    } else if ( mutation == Mutation::CLIQUE_SPLITTER ) {
      return utils::Randomize::instance().getRandomFloat(
        _context.refinement.evo.min_clique_split_mutation_prob,
        _context.refinement.evo.max_clique_split_mutation_prob);
    } else if ( mutation == Mutation::TEST_MUTATION ) {
      return utils::Randomize::instance().getRandomFloat(
        _context.refinement.evo.min_test_mutation_prob,
        _context.refinement.evo.max_test_mutation_prob);
    }
    return 0.0f;
  }

  const Context& _context;
  const bool _show_detailed_output;
  ActionSelector<Mutation> _mutation_selector;
};

}  // namespace cluster_editing
