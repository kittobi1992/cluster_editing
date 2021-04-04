
#include <optional>

#include <cluster_editing/exact/instance.h>

int kab_bound(const Instance& inst, int limit, bool verbose = false, int time_limit=0);

std::optional<Instance> forcedChoicesKAB(const Instance& inst, int upper_bound, bool verbose=false);

