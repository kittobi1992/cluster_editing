#include "evolutionary.h"

#include "cluster_editing/utils/randomize.h"
#include "cluster_editing/metrics.h"


namespace cluster_editing::evolutionary {

void EvolutionaryAlgorithm::refine() {
  auto& rnd = utils::Randomize::instance();
  auto& prng = rnd.getGenerator();

  bool use_lp = true;
  bool use_fm = std::bernoulli_distribution(0.05)(prng);    // with 5% probability, run FM
  //use_fm = false;
  if (use_fm) {
    // disable LP with 5% prob, if FM is used
    use_lp = std::bernoulli_distribution(0.95)(prng);
  }

  if (use_lp) {
    auto& lp_opt = shared_context.refinement.lp;
    auto old_lp_opt = lp_opt; // make copy

    lp_opt.maximum_lp_iterations = 10000;
    lp_opt.node_order = NodeOrdering::random_shuffle;
    lp_opt.random_shuffle_each_round = true;
    lp_opt.activate_all_cliques_after_rounds = 5;

    lp_refiner.initialize(graph);
    lp_refiner.refine(graph);

    lp_opt = old_lp_opt;
  }
  if (use_fm) {
    size_t old_fm_iters = shared_context.refinement.maximum_fm_iterations;
    shared_context.refinement.maximum_fm_iterations = 5 + std::geometric_distribution<size_t>(0.3)(prng);

    fm_refiner.initialize(graph);
    fm_refiner.refine(graph);

    shared_context.refinement.maximum_fm_iterations = old_fm_iters;
  }
}


void EvolutionaryAlgorithm::intersect_cliques(const CliqueAssignment& cliques1, const CliqueAssignment& cliques2) {
  CliqueID z = *std::max_element(cliques1.begin(), cliques1.end()) + 1;
  CliqueID fresh = 0;
  compactification_map.clear();
  for (NodeID u : graph.nodes()) {
    CliqueID cu = cliques2[u] * z + cliques1[u];
    auto [it, insertion_happened] = compactification_map.try_emplace(cu, fresh);
    if (insertion_happened) fresh++;
    clique_intersection[u] = it->second;
  }
}


Edits EvolutionaryAlgorithm::convert_clique_assignment_to_edits(const CliqueAssignment& cliques) {
  // group vertices by cluster
  CliqueID num_cliques = *std::max_element(cliques.begin(), cliques.end()) + 1;
  std::vector<uint32_t> clique_start(num_cliques + 2, 0);
  std::vector<NodeID> vertices_by_clique(graph.numNodes());
  for (NodeID u : graph.nodes()) {
    clique_start[cliques[u] + 2]++;
  }
  for (size_t i = 3; i < clique_start.size(); ++i) {
    clique_start[i] += clique_start[i - 1];
  }
  for (NodeID u : graph.nodes()) {
    vertices_by_clique[ clique_start[cliques[u]+1]++ ] = u;
  }

  // enumerate vertex pairs in each cluster. if not adjacent --> insertion.
  // edges to outside clusters --> deletion
  std::vector<bool> neighbors(graph.numNodes(), false);
  std::vector<Edit> edits;
  for (CliqueID c = 0; c < num_cliques; ++c) {
    for (size_t i = clique_start[c]; i < clique_start[c+1]; ++i) {
      const NodeID u = vertices_by_clique[i];

      for (const Neighbor& nb : graph.neighbors(u)) {
        const NodeID v = nb.target;
        if (u < v && c != cliques[v]) {
          edits.emplace_back(u, v); // deletion
        }
        neighbors[v] = true;
      }

      for (size_t j = i + 1; j < clique_start[c+1]; ++j) {
        const NodeID v = vertices_by_clique[j];
        if (!neighbors[v]) {
          edits.emplace_back(u, v);  // insertion
        }
      }

      for (const Neighbor& nb : graph.neighbors(u)) {
        neighbors[nb.target] = false;
      }
    }
  }
  return edits;
}

size_t EvolutionaryAlgorithm::num_different_edits(const Edits& edits1, const Edits& edits2) {
  assert(std::is_sorted(edits1.begin(), edits1.end()));
  assert(std::is_sorted(edits2.begin(), edits2.end()));

  size_t i1 = 0, i2 = 0, n = 0;
  while (i1 < edits1.size() && i2 < edits2.size()) {
    if (edits1[i1] < edits2[i2]) {
      ++n; ++i1;
    } else if (edits1[i1] > edits2[i2]) {
      ++n; ++i2;
    } else {
      ++i1; ++i2;
    }
  }
  n += edits1.size() - i1;
  n += edits2.size() - i2;

#ifndef NDEBUG
  Edits sym_diff;
  std::set_symmetric_difference(edits1.begin(), edits1.end(), edits2.begin(), edits2.end(), std::back_inserter(sym_diff));
  assert(n == sym_diff.size());
#endif

  return n;
}

std::pair<size_t, double> EvolutionaryAlgorithm::rand_index(const CliqueAssignment& cliques1, const CliqueAssignment& cliques2) {
  intersect_cliques(cliques1, cliques2);

  intersection_size.assign(graph.numNodes(), 0);
  c1_size.assign(graph.numNodes(), 0);
  c2_size.assign(graph.numNodes(), 0);

  for (NodeID u : graph.nodes()) {
    c1_size[cliques1[u]]++;
    c2_size[cliques2[u]]++;
    intersection_size[clique_intersection[u]]++;
    // could update sums here
  }

  size_t c1_sum = 0, c2_sum = 0, intersection_sum = 0;
  for (NodeID u : graph.nodes()) {
    c1_sum += c1_size[u]*(c1_size[u]-1)/2;
    c2_sum += c2_size[u]*(c2_size[u]-1)/2;
    intersection_sum += intersection_size[u]*(intersection_size[u]-1)/2;
  }

  size_t n = graph.numNodes();
  size_t pairs = n*(n-1)/2;
  size_t agreements = pairs + 2 * intersection_sum - c1_sum - c2_sum;

  double max_index = 0.5 * (c1_sum + c2_sum);
  double expected_index = static_cast<double>(c1_sum * c2_sum) / static_cast<double>(pairs);
  double ari;
  if ((c1_sum == 0 && c2_sum == 0) || max_index == expected_index) {
    ari = 1.0;
  } else {
    ari = (intersection_sum - expected_index) / (max_index - expected_index);
  }

  return std::make_pair(agreements,ari);
}

void EvolutionaryAlgorithm::generate_initial_population() {
  for (size_t i = 0; i < max_pop_size; ++i) {
    graph.reset();
    refine();
    add_current_solution();
  }
}

void EvolutionaryAlgorithm::evolution_step() {
  static constexpr double recombine_probability = 0.7;
  static constexpr double mutation_probability = 0.15;
  static constexpr double from_scratch_probability = 0.15;
  assert(recombine_probability + mutation_probability + from_scratch_probability == 1.0);
  double toss = utils::Randomize::instance().getRandomFloat(0.0, 1.0);
  if (toss < recombine_probability) {
    LOG << "recombine";
    recombine();
  } else if (toss < recombine_probability + mutation_probability) {
    LOG << "mutate";
    mutation();
  } else {
    LOG << "scratch";
    graph.reset();
  }

  refine();
  add_current_solution();
}

void EvolutionaryAlgorithm::add_current_solution() {
  size_t cost = metrics::edits(graph);

  if (population.size() < max_pop_size || cost < worst().cost) {
    Solution sol;
    sol.cliques.reserve(graph.numNodes());
    for (NodeID u : graph.nodes()) {
      sol.cliques.push_back(graph.clique(u));
    }
    sol.cost = cost;

    if (population.size() < max_pop_size) {
      population.emplace_back(std::move(sol));
    } else {
      size_t replace = population.size();
      size_t min_agreement = graph.numNodes()*graph.numNodes();
      double min_ari = 2.0;
      for (size_t i = 0; i < population.size(); ++i) {
        if (population[i].cost >= sol.cost) {
          auto [agreement, ari] = rand_index(population[i].cliques, sol.cliques);
          if (std::tie(ari, agreement) < std::tie(min_ari, min_agreement)) {
            min_ari = ari; min_agreement = agreement; replace = i;
          }
        }
      }

      if (replace != population.size()) {
        population[replace] = std::move(sol);
      }

    }
  }
}

void EvolutionaryAlgorithm::mutation() {
  size_t sol_pos;
  auto& rnd = utils::Randomize::instance();
  if (population.empty()) {
    graph.reset();
    return;
  } else if (population.size() < 2) {
    sol_pos = 0;
  } else {
    // randomly selected, ro2 tournament
    size_t i1 = rnd.getRandomInt(0, population.size() - 1);
    size_t i2 = rnd.getRandomInt(0, population.size() - 2);
    if (i1 == i2) {
      i2 = population.size() - 1;
    }
    if (population[i1].cost < population[i2].cost
        || (population[i1].cost == population[i2].cost && rnd.flipCoin())) {
      sol_pos = i1;
    } else {
      sol_pos = i2;
    }
  }

  const CliqueAssignment& c = population[sol_pos].cliques;
  for (NodeID u : graph.nodes()) {
    graph.setClique(u, c[u]);
  }

  // TODO store whether solution in population was already compactified?
  CliqueID num_cliques = 0;

  auto compactify = [&] {
    num_cliques = 0;
    compactification_map.clear();
    for (NodeID u : graph.nodes()) {
      auto [it, insertion_happened] = compactification_map.try_emplace(graph.clique(u), num_cliques);
      if (insertion_happened) num_cliques++;
      graph.setClique(u, it->second);
    }
    assert(num_cliques <= graph.numNodes());
  };

  auto& prng = rnd.getGenerator();
  static constexpr size_t num_options = 4;
  std::bernoulli_distribution toss(1.0 / (num_options - 1.0));   // n * (1 / n-1) expected perturbations
  std::array<bool, num_options> options { false };
  bool any = false;
  for (bool& o : options) {
    o = toss(prng);
    any |= o;
  }
  if (!any) {
    options[std::uniform_int_distribution<int>(0, num_options - 1)(prng)] = true;
  }

  compactify();

  if (options[0]) { // split_clique_into_sampled_pieces
    CliqueID clique_to_split = rnd.getRandomInt(0, num_cliques - 1);
    std::vector<NodeID> nodes_in_clique;
    for (NodeID u : graph.nodes()) {
      if (graph.clique(u) == clique_to_split) {
        nodes_in_clique.push_back(u);
      }
    }
    if (nodes_in_clique.size() > 1) {
      size_t ub = std::ceil(std::sqrt(nodes_in_clique.size()));
      size_t num_splits = rnd.getRandomInt(2, ub);
      for (NodeID u : nodes_in_clique) {
        graph.setClique(u, num_cliques + rnd.getRandomInt(0, num_splits - 1));
      }
      num_cliques += num_splits;
    }
  }

  if (options[1]) { // isolate_clique
    CliqueID clique_to_split = rnd.getRandomInt(0, num_cliques - 1);
    for (NodeID u : graph.nodes()) {
      if (graph.clique(u) == clique_to_split) {
        graph.setClique(u, num_cliques++);
      }
    }
  }

  if (options[2]) { // perform random moves
    size_t sqrt_n = std::ceil(std::sqrt(graph.numNodes()));
    size_t num_moves =  std::min(100, rnd.getRandomInt(sqrt_n, 2*sqrt_n));
    for (size_t i = 0; i < num_moves; ++i) {
      NodeID u = rnd.getRandomInt(0, graph.numNodes());   // don't care about duplicates
      CliqueID target = rnd.getRandomInt(0, num_cliques);
      if (target == num_cliques) {
        // isolate
        num_cliques++;
      }
      graph.setClique(u, target);
    }
  }

  if (options[3]) { // isolate_random_nodes
    size_t sqrt_n = std::ceil(std::sqrt(graph.numNodes()));
    size_t num_moves =  std::min(100, rnd.getRandomInt(sqrt_n, 2*sqrt_n));
    for (size_t i = 0; i < num_moves; ++i) {
      NodeID u = rnd.getRandomInt(0, graph.numNodes());
      graph.setClique(u, num_cliques++);
    }
  }

  if (num_cliques > graph.numNodes()) {
    /*
    LOG << "had to compactify after perturbation";
    for (bool option : options) {
      std::cout << std::boolalpha << option << " ";
    }
    std::cout << std::endl;
     */
    compactify();
  }
}

void EvolutionaryAlgorithm::recombine() {
  if (population.size() < 3) {
    graph.reset();
    return;
  } else {
    auto& rnd = utils::Randomize::instance();
    std::array<size_t, 4> starters {};
    for (size_t i = 0; i < starters.size(); ++i) {
      size_t pos = rnd.getRandomInt(0, population.size() - 1 - i);

      // check for duplicates and replace with fall-back values, if so
      for (size_t j = 0; j < i; ++j) {
        if (pos == starters[j]) {
          pos = population.size() - 1 - j;
        }
      }
      starters[i] = pos;
    }

    auto select = [&](size_t p1, size_t p2) {
      size_t cost1 = population[starters[p1]].cost, cost2 = population[starters[p2]].cost;
      if (cost1 < cost2 || (cost1 == cost2 && rnd.flipCoin())) {
        return starters[p1];
      } else {
        return starters[p2];
      }
    };

    const CliqueAssignment& c1 = population[select(0,1)].cliques;
    const CliqueAssignment& c2 = population[select(2,3)].cliques;
    intersect_cliques(c1, c2);
    for (NodeID u : graph.nodes()) {
      graph.setClique(u, clique_intersection[u]);
    }
  }
}

}