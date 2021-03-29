
#include <cluster_editing/exact/instance.h>
#include <optional>

struct StarBoundConfig {
    size_t max_num_unchanged = 5;
    float min_degree_candidate_vs_random_candidate_probability = 0.8;
};

inline StarBoundConfig global_star_bound_config;

int star_bound(const Instance& inst, int limit);

std::optional<Instance> forcedChoicesStarBound(const Instance& inst, int upper_bound, bool verbose);

std::optional<Instance> forcedChoicesSingleMerge(const Instance& inst, int upper_bound, bool verbose);
