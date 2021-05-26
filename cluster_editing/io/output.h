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

#include <chrono>

#include "cluster_editing/definitions.h"
#include "cluster_editing/context/context.h"

namespace cluster_editing::io {

  void printObjectives(const Graph& graph,
                       const std::chrono::duration<double>& elapsed_seconds);
  void printClusterEditingResults(const Graph& hypergraph,
                                  const Context& context,
                                  const std::chrono::duration<double>& elapsed_seconds);
  void printResultLine(const Graph& graph,
                       const Context& context,
                       const std::chrono::duration<double>& elapsed_seconds);

  void printStripe();
  void printBanner();
  void printContext(const Context& context);
  void printInputInfo(const Graph& graph, const Context& context);
  void printGraphInfo(const Graph& graph, const Context& context, const std::string& name);
  void printCliqueInfo(const Graph& graph, const Context& context);
  void printClusteringBanner(const Context& context);
  void printInitialSolutionBanner(const Context& context);
  void printGlobalEvoBanner(const Context& context);
  void printLocalizedEvoBanner(const Context& context);


  void readSolutionFile(Graph& graph, const std::string& filename);
  void writeSolutionFile(const Graph& graph, const std::string& filename);
}  // namespace cluster_editing::io
