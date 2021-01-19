
#pragma once

#include <vector>

using Edges = std::vector<std::vector<int>>; // adj matrix
using IDMap = std::vector<std::vector<int>>; // for every node in current graph, a list of original nodes that it represents

struct Instance {
  Edges edges;
  IDMap idmap;

  int spendCost = 0;

  Instance(const Edges& _edges, const IDMap& _idmap) : edges(_edges), idmap(_idmap) {};

  Instance(int n) : edges(n, std::vector<int>(n,-1)), idmap(n) {
    for(int i=0; i<n; ++i) idmap[i] = {i};
  }
};