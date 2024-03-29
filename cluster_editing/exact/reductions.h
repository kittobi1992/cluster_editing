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

#include "instance.h"

#include <utility>
#include <optional>

std::vector<int> connectedComponents(const Edges& graph);

std::vector<Instance> constructConnectedComponents(const Instance& graph);

// contracts the edge u,v
Instance merge(const Instance& inst, int u, int v);
void mergeAllINF(Instance& inst);

// compute icp and ice (lower bound for making an edge permanent/forbidden)
std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>> computeICFandICP(const Edges& edges);

// check for unaffordable edge modifications (params by reference because they contain updated values after reducton)
std::optional<Instance> icxReductions(const Instance& inst, int budget);


// set all pairs with dist>=4 to forbidden
std::optional<Instance> distance4Reduction(const Instance& inst);

std::optional<Instance> forcedChoices(const Instance& inst, int upper_bound, bool verbose=false);

std::optional<Instance> simpleTwin(const Instance& inst);

std::optional<Instance> complexTwin(const Instance& inst, bool calc_dp);

std::optional<Instance> mergeCliques(const Instance& inst);

std::optional<Instance> weightedKernel(const Instance& inst);

std::optional<Instance> force_small_components(const Instance& inst);
