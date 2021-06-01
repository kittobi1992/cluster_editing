/*******************************************************************************
 * This file is part of KaPoCE.
 *
 * KaPoCE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaPoCE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaPoCE.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include <cluster_editing/exact/instance.h>
#include <optional>

struct StarBoundConfig {
    size_t max_num_unchanged = 5;
    float min_degree_candidate_vs_random_candidate_probability = 0.8;
};

inline StarBoundConfig global_star_bound_config;

int star_bound(const Instance& inst, int limit);

std::optional<Instance> forcedChoicesStarBound(const Instance& inst, int upper_bound, bool verbose, int min_time=0);

std::optional<Instance> forcedChoicesSingleMerge(const Instance& inst, int upper_bound, bool verbose);
