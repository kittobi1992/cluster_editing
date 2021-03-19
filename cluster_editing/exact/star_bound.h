
#include <cluster_editing/exact/instance.h>
#include <optional>

int star_bound(const Instance& inst, int limit);

std::optional<Instance> forcedChoicesStarBound(const Instance& inst, int upper_bound, bool verbose);
