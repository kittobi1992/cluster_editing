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

#pragma once

#include <vector>

using Edges = std::vector<std::vector<int>>; // adj matrix
using IDMap = std::vector<std::vector<int>>; // for every node in current graph, a list of original nodes that it represents
using Clust = std::vector<std::vector<int>>; // clusters that were already solved

const int INF = 1e9;

struct Instance {
  Edges orig;
  Edges edges;
  IDMap idmap;

  int spendCost = 0;
  Clust done_clusters;

  Instance() : spendCost(1e9) {}; // unsolvable instance
  Instance(const Edges& _edges, const IDMap& _idmap) : edges(_edges), idmap(_idmap) {};

  Instance(int n) : edges(n, std::vector<int>(n,-1)), idmap(n) {
    for(int i=0; i<n; ++i) idmap[i] = {i};
    orig = edges;
  }
};

int forbiddenEdges(const Instance& inst);

Instance load_exact_instance(int num);
Instance load_exact_instance(); // from stdin

Instance remove_nodes(const Instance& inst, const std::vector<int>& nodes);