#pragma once

#include "cluster_editing/refinement/lp_refiner.h"
#include "cluster_editing/refinement/fm_refiner.h"

#include <unordered_map>

namespace cluster_editing::evolutionary {

using CliqueAssignment = std::vector<CliqueID>;
using Edits = std::vector<Edit>;

struct Solution {
  CliqueAssignment cliques;
  size_t cost;
  size_t together_pairs;  // TODO cache pairs stuck together in this solution
};

class EvolutionaryAlgorithm {
public:
  EvolutionaryAlgorithm(Graph& graph, Context& context) :
      graph(graph),
      shared_context(context),
      lp_refiner(graph, context),
      fm_refiner(graph, context),
      clique_intersection(graph.numNodes(), INVALID_CLIQUE),
      intersection_size(graph.numNodes(), 0),
      c1_size(graph.numNodes(), 0),
      c2_size(graph.numNodes(), 0)
      { }

  void generate_initial_population();

  void evolution_step();

  void mutation();

  void recombine();

  void refine();

  void intersect_cliques(const CliqueAssignment& cliques1, const CliqueAssignment& cliques2);

  Edits convert_clique_assignment_to_edits(const CliqueAssignment& cliques);

  size_t num_different_edits(const Edits& edits1, const Edits& edits2);

  // plain agreements and adjusted rand index
  std::pair<size_t, double> rand_index(const CliqueAssignment& cliques1, const CliqueAssignment& cliques2);

  void print_best_solution() {
    const auto& best_solution = best();
    Edits edits = convert_clique_assignment_to_edits(best_solution.cliques);
    for (const Edit& e : edits) {
      std::cout << e.first << " " << e.second << "\n";
    }
    std::cout << std::flush;
  }

  void apply_best_solution() {
    const auto& best_solution = best();
    for (NodeID u : graph.nodes()) {
      graph.setClique(u, best_solution.cliques[u]);
    }
  }

  void add_current_solution();

  size_t max_pop_size = 10;

  const Solution& worst() const {
    auto comp = [](const Solution& s1, const Solution& s2) { return s1.cost < s2.cost; };
    return *std::max_element(population.begin(), population.end(), comp);
  }

  const Solution& best() const {
    auto comp = [](const Solution& s1, const Solution& s2) { return s1.cost < s2.cost; };
    return *std::min_element(population.begin(), population.end(), comp);
  }

private:
  Graph& graph;
  std::vector<Solution> population;

  Context& shared_context;
  LabelPropagationRefiner lp_refiner;
  FMRefiner fm_refiner;

  std::unordered_map<CliqueID, CliqueID> compactification_map;
  CliqueAssignment clique_intersection;

  std::vector<size_t> intersection_size, c1_size, c2_size;
};

}