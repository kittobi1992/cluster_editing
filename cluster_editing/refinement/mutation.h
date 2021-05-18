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
#include "cluster_editing/metrics.h"

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
  static void mutate(Graph& graph, const EdgeWeight current_edits, const float prob) {
    utils::CommonOperations::instance(graph).computeEmptyCliques(graph);
    std::vector<CliqueID>& empty_cliques =
      utils::CommonOperations::instance(graph)._empty_cliques;
    std::vector<NodeID>& cluster_sizes =
      utils::CommonOperations::instance(graph)._cluster_sizes;
    std::vector<NodeID>& nodes =
      utils::CommonOperations::instance(graph)._nodes;
    std::random_shuffle(nodes.begin(), nodes.end());
    EdgeWeight overall_delta = 0;
    const EdgeWeight budget = static_cast<float>(current_edits) * prob;
    for ( const NodeID& u : nodes ) {
      const NodeID u_degree = graph.degree(u);
      const CliqueID from = graph.clique(u);
      if ( cluster_sizes[from] > 1 ) {
        EdgeWeight edges_to_from = 0;
        for ( const NodeID& v : graph.neighbors(u) ) {
          if ( graph.clique(v) == from ) {
            ++edges_to_from;
          }
        }

        const EdgeWeight from_rating = cluster_sizes[from] - 1 + u_degree - 2 * edges_to_from;
        const EdgeWeight isolation_delta = u_degree - from_rating;
        overall_delta += isolation_delta;
        const CliqueID to = empty_cliques.back();
        empty_cliques.pop_back();
        graph.setClique(u, to);
        --cluster_sizes[from];
        ++cluster_sizes[to];
      }

      if ( overall_delta > budget ) {
        break;
      }
    }
  }

 private:
  RandomNodeIsolator() { }
};

class RandomNodeMover {

 public:
  static void mutate(Graph& graph, const EdgeWeight current_edits, const float prob) {
    utils::CommonOperations::instance(graph).computeClusterSizes(graph);
    ds::FixedSizeSparseMap<CliqueID, EdgeWeight>& rating =
      utils::CommonOperations::instance(graph)._rating;
    std::vector<NodeID>& cluster_sizes =
      utils::CommonOperations::instance(graph)._cluster_sizes;
    std::vector<NodeID>& nodes =
      utils::CommonOperations::instance(graph)._nodes;
    std::random_shuffle(nodes.begin(), nodes.end());
    EdgeWeight overall_delta = 0;
    const EdgeWeight budget = static_cast<float>(current_edits) * prob;
    std::vector<CliqueID> target_cliques;
    for ( const NodeID& u : graph.nodes() ) {
      const float p = utils::Randomize::instance().getRandomFloat(0.0, 1.0);
      if ( p <= prob ) {
        rating.clear();
        const NodeID u_degree = graph.degree(u);
        const CliqueID from = graph.clique(u);
        for ( const NodeID& v : graph.neighbors(u) ) {
          const CliqueID to = graph.clique(v);
          ++rating[to];
          if ( from != to ) {
            target_cliques.push_back(to);
          }
        }

        const EdgeWeight from_rating = cluster_sizes[from] - 1 + u_degree - 2 * rating[from];
        if ( !target_cliques.empty() ) {
          std::sort(target_cliques.begin(), target_cliques.end());
          target_cliques.erase( std::unique(
            target_cliques.begin(), target_cliques.end() ), target_cliques.end() );
          const CliqueID to = target_cliques[
            utils::Randomize::instance().getRandomInt(0, target_cliques.size() - 1)];
          const EdgeWeight to_rating = cluster_sizes[to] + u_degree - 2 * rating[to];
          const EdgeWeight move_delta = to_rating - from_rating;
          overall_delta += move_delta;
          graph.setClique(u, to);
          ++cluster_sizes[to];
          --cluster_sizes[from];
          target_cliques.clear();
        }
      }
      if ( overall_delta > budget ) {
        break;
      }
    }
  }

 private:
  RandomNodeMover() { }
};

class TestMutation {

 public:
  static void mutate(Graph& graph, const EdgeWeight current_edits, const float prob) {
    utils::CommonOperations::instance(graph).computeEmptyCliques(graph);
    ds::FixedSizeSparseMap<CliqueID, EdgeWeight>& rating =
      utils::CommonOperations::instance(graph)._rating;
    std::vector<CliqueID>& empty_cliques =
      utils::CommonOperations::instance(graph)._empty_cliques;
    std::vector<NodeID>& cluster_sizes =
      utils::CommonOperations::instance(graph)._cluster_sizes;
    std::vector<NodeID>& nodes =
      utils::CommonOperations::instance(graph)._nodes;
    std::random_shuffle(nodes.begin(), nodes.end());
    EdgeWeight overall_delta = 0;
    const EdgeWeight budget = static_cast<float>(current_edits) * prob;
    std::vector<CliqueID> target_cliques;
    for ( const NodeID& u : graph.nodes() ) {
      const float p = utils::Randomize::instance().getRandomFloat(0.0, 1.0);
      if ( p <= prob ) {
        const bool is_move = utils::Randomize::instance().flipCoin();
        rating.clear();
        const NodeID u_degree = graph.degree(u);
        const CliqueID from = graph.clique(u);
        for ( const NodeID& v : graph.neighbors(u) ) {
          const CliqueID to = graph.clique(v);
          ++rating[to];
          if ( from != to && is_move ) {
            target_cliques.push_back(to);
          }
        }

        const EdgeWeight from_rating = cluster_sizes[from] - 1 + u_degree - 2 * rating[from];
        EdgeWeight delta = 0;
        CliqueID to = INVALID_CLIQUE;
        if ( is_move && !target_cliques.empty() ) {
          std::sort(target_cliques.begin(), target_cliques.end());
          target_cliques.erase( std::unique(
            target_cliques.begin(), target_cliques.end() ), target_cliques.end() );
          const CliqueID to = target_cliques[
            utils::Randomize::instance().getRandomInt(0, target_cliques.size() - 1)];
          const EdgeWeight to_rating = cluster_sizes[to] + u_degree - 2 * rating[to];
          delta = to_rating - from_rating;
        } else if ( !is_move && cluster_sizes[from] > 1 ) {
          delta = u_degree - from_rating;
          to = empty_cliques.back();
          empty_cliques.pop_back();
        }

        if ( to != INVALID_CLIQUE ) {
          overall_delta += delta;
          graph.setClique(u, to);
          ++cluster_sizes[to];
          --cluster_sizes[from];
          target_cliques.clear();
        }
      }
      if ( overall_delta > budget ) {
        break;
      }
    }
  }

 private:
  TestMutation() { }
};

class Mutator {

  struct MutationProbs {
    float clique_isolation_prob;
    float neighbor_clique_isolation_prob;
    float node_isolation_prob;
    float node_move_prob;
    float node_move_or_isolate_prob;
    float last_prob;
  };

 public:
  explicit Mutator(const Context& context) :
    _context(context),
    _show_detailed_output(context.general.verbose_output && context.refinement.evo.enable_detailed_output),
    _probs(),
    _mutation_selector() {
    activateMutations(_context.refinement.evo.enabled_mutations);
    _probs.clique_isolation_prob = _context.refinement.evo.max_clique_isolate_prob;
    _probs.neighbor_clique_isolation_prob = _context.refinement.evo.max_neighbor_clique_isolate_prob;
    _probs.node_isolation_prob = _context.refinement.evo.max_node_isolation_prob;
    _probs.node_move_prob = _context.refinement.evo.max_node_move_prob;
    _probs.node_move_or_isolate_prob = _context.refinement.evo.max_test_mutation_prob;
  }

  Mutation mutate(Graph& graph, const EdgeWeight current_edits) {
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
      RandomNodeIsolator::mutate(graph, current_edits, prob);
    } else if ( mutation == Mutation::RANDOM_NODE_MOVER ) {
      if ( _show_detailed_output )
        LOG << "Mutation Action: RANDOM_NODE_MOVER ( p =" << prob << ")";
      RandomNodeMover::mutate(graph, current_edits, prob);
    } else if ( mutation == Mutation::TEST_MUTATION ) {
      if ( _show_detailed_output )
        LOG << "Mutation Action: TEST_MUTATION ( p =" << prob << ")";
      TestMutation::mutate(graph, current_edits, prob);
    }
    _probs.last_prob = prob;
    return mutation;
  }

  void updateProbs(const Mutation& mutation,
                   const EdgeWeight before_edits,
                   const EdgeWeight after_edits) {
    const EdgeWeight delta = (after_edits - before_edits);
    const bool improvement = delta < 0;
    const float factor = (delta < 0 ? 2.0 : ( delta == 0 ? 1.0 : 0.5 ) );
    auto update_prob = [&](float& prob, const float min_prob, const float max_prob) {
      if ( improvement && prob < _probs.last_prob ) {
        prob = _probs.last_prob;
      } else if ( prob == _probs.last_prob ) {
        prob = std::min(max_prob, std::max( min_prob, prob * factor ));
      }
    };

    if ( mutation == Mutation::LARGE_CLIQUE_ISOLATOR ) {
      update_prob(_probs.clique_isolation_prob,
        _context.refinement.evo.min_clique_isolate_prob,
        _context.refinement.evo.max_clique_isolate_prob);
    } else if ( mutation == Mutation::LARGE_CLIQUE_WITH_NEIGHBOR_ISOLATOR ) {
      update_prob(_probs.neighbor_clique_isolation_prob,
        _context.refinement.evo.min_neighbor_clique_isolate_prob,
        _context.refinement.evo.max_neighbor_clique_isolate_prob);
    } else if ( mutation == Mutation::RANDOM_NODE_ISOLATOR ) {
      update_prob(_probs.node_isolation_prob,
        _context.refinement.evo.min_node_isolation_prob,
        _context.refinement.evo.max_node_isolation_prob);
    } else if ( mutation == Mutation::RANDOM_NODE_MOVER ) {
      update_prob(_probs.node_move_prob,
        _context.refinement.evo.min_node_move_prob,
        _context.refinement.evo.max_node_move_prob);
    } else if ( mutation == Mutation::TEST_MUTATION ) {
      update_prob(_probs.node_move_or_isolate_prob,
        _context.refinement.evo.min_test_mutation_prob,
        _context.refinement.evo.max_test_mutation_prob);
    }

    _mutation_selector.notifyImprovement(mutation, std::max( 0, before_edits - after_edits ));
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
    const bool select_random_probability =
      utils::Randomize::instance().getRandomFloat(0.0, 1.0) <=
      _context.refinement.evo.random_prob_selection_prob;
    if ( mutation == Mutation::LARGE_CLIQUE_ISOLATOR ) {
      if ( !select_random_probability ) {
        return _probs.clique_isolation_prob;
      } else {
        return utils::Randomize::instance().getRandomFloat(
          _context.refinement.evo.min_clique_isolate_prob,
          _context.refinement.evo.max_clique_isolate_prob);
      }
    } else if ( mutation == Mutation::LARGE_CLIQUE_WITH_NEIGHBOR_ISOLATOR ) {
      if ( !select_random_probability ) {
        return _probs.neighbor_clique_isolation_prob;
      } else {
        return utils::Randomize::instance().getRandomFloat(
          _context.refinement.evo.min_neighbor_clique_isolate_prob,
          _context.refinement.evo.max_neighbor_clique_isolate_prob);
      }
      return _probs.neighbor_clique_isolation_prob;
    } else if ( mutation == Mutation::RANDOM_NODE_ISOLATOR ) {
      if ( !select_random_probability ) {
        return _probs.node_isolation_prob;
      } else {
        return utils::Randomize::instance().getRandomFloat(
          _context.refinement.evo.min_node_isolation_prob,
          _context.refinement.evo.max_node_isolation_prob);
      }
    } else if ( mutation == Mutation::RANDOM_NODE_MOVER ) {
      if ( !select_random_probability ) {
        return _probs.node_move_prob;
      } else {
        return utils::Randomize::instance().getRandomFloat(
          _context.refinement.evo.min_node_move_prob,
          _context.refinement.evo.max_node_move_prob);
      }
    } else if ( mutation == Mutation::TEST_MUTATION ) {
      if ( !select_random_probability ) {
        return _probs.node_move_or_isolate_prob;
      } else {
        return utils::Randomize::instance().getRandomFloat(
          _context.refinement.evo.min_test_mutation_prob,
          _context.refinement.evo.max_test_mutation_prob);
      }
    }
    return 0.0f;
  }

  const Context& _context;
  const bool _show_detailed_output;
  MutationProbs _probs;
  ActionSelector<Mutation> _mutation_selector;
};

}  // namespace cluster_editing
