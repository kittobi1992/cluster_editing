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
  static void mutate(Graph& graph, const Context& context, const bool prob) {
    CliqueStats stats(graph);
    for ( const CliqueID& c : graph.nodes() ) {
      if ( stats.cliques[c].size() >= context.refinement.evo.large_clique_threshold ) {
        const float p = utils::Randomize::instance().getRandomFloat(0.0, 1.0);
        if ( p <= prob ) {
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
  static void mutate(Graph& graph, const Context& context, const float prob) {
    CliqueStats stats(graph);
    for ( const CliqueID& c : graph.nodes() ) {
      if ( stats.cliques[c].size() >= context.refinement.evo.large_clique_threshold ) {
        const float p = utils::Randomize::instance().getRandomFloat(0.0, 1.0);
        if ( p <= prob ) {
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

class RandomNodeIsolator {

 public:
  static void mutate(Graph& graph, const float prob) {
    CliqueStats stats(graph);
    for ( const CliqueID& c : graph.nodes() ) {
      if ( stats.cliques[c].size() > 1 ) {
        size_t current_size = stats.cliques[c].size();
        for ( const NodeID& u : stats.cliques[c].nodes ) {
          const float p = utils::Randomize::instance().getRandomFloat(0.0, 1.0);
          if ( p <= prob ) {
            const CliqueID target = stats.empty_cliques.back();
            stats.empty_cliques.pop_back();
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
        for ( const Neighbor& n : graph.neighbors(u) ) {
          const CliqueID to = graph.clique(n.target);
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


enum class Mutation : uint8_t {
  LARGE_CLIQUE_ISOLATOR = 0,
  LARGE_CLIQUE_WITH_NEIGHBOR_ISOLATOR = 1,
  RANDOM_NODE_ISOLATOR = 2,
  RANDOM_NODE_MOVER = 3,
  NUM_MUTATIONS = 4
};


class Mutator {

  struct MutationProbs {
    float clique_isolation_prob;
    float neighbor_clique_isolation_prob;
    float node_isolation_prob;
    float node_move_prob;
    float last_prob;
  };

 public:
  explicit Mutator(const Context& context) :
    _context(context),
    _show_detailed_output(context.general.verbose_output && context.refinement.evo.enable_detailed_output),
    _probs(),
    _active_mutations() {
    activateMutations(_context.refinement.evo.enabled_mutations);
    _probs.clique_isolation_prob = _context.refinement.evo.max_clique_isolate_prob;
    _probs.neighbor_clique_isolation_prob = _context.refinement.evo.max_neighbor_clique_isolate_prob;
    _probs.node_isolation_prob = _context.refinement.evo.max_node_isolation_prob;
    _probs.node_move_prob = _context.refinement.evo.max_node_move_prob;
  }

  Mutation mutate(Graph& graph) {
    const Mutation mutation = chooseMutation();
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
        _context.refinement.evo.max_clique_isolate_prob);
    } else if ( mutation == Mutation::RANDOM_NODE_ISOLATOR ) {
      update_prob(_probs.node_isolation_prob,
        _context.refinement.evo.min_node_isolation_prob,
        _context.refinement.evo.max_node_isolation_prob);
    } else if ( mutation == Mutation::RANDOM_NODE_MOVER ) {
      update_prob(_probs.node_move_prob,
        _context.refinement.evo.min_node_move_prob,
        _context.refinement.evo.max_node_move_prob);
    }
  }

  void activateMutations(const std::string& enabled_mutations) {
    ASSERT(enabled_mutations.size() == static_cast<size_t>(Mutation::NUM_MUTATIONS));
    _active_mutations.clear();
    for ( size_t i = 0; i < enabled_mutations.size(); ++i ) {
      if ( enabled_mutations[i] == '1' ) {
        _active_mutations.push_back(static_cast<Mutation>(i));
      }
    }
  }

 private:
  Mutation chooseMutation() const {
    const int rnd_mutation = utils::Randomize::instance().getRandomInt(
      0, _active_mutations.size() - 1);
    return _active_mutations[rnd_mutation];
  }

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
        return _probs.node_isolation_prob;
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
    }
    return 0.0f;
  }

  const Context& _context;
  const bool _show_detailed_output;
  MutationProbs _probs;
  std::vector<Mutation> _active_mutations;
};

}  // namespace cluster_editing
