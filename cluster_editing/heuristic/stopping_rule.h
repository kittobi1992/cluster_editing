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

#include <cmath>

#include "cluster_editing/definitions.h"
#include "cluster_editing/context/context.h"

namespace cluster_editing {

class FruitlessMovesStoppingRule {

 public:
  FruitlessMovesStoppingRule(const Graph&, const Context& context) :
    _context(context),
    _round_delta(0),
    _best_delta(0),
    _num_fruitless_moves(0) { }


  bool searchShouldStop() {
    return _num_fruitless_moves > _context.refinement.fm.max_fruitless_moves;
  }

  void update(const EdgeWeight delta) {
    _round_delta += delta;
    if ( _round_delta < _best_delta ) {
      _best_delta = _round_delta;
      _num_fruitless_moves = 0;
    } else {
      ++_num_fruitless_moves;
    }
  }

  void reset() {
    _num_fruitless_moves = 0;
  }

 private:
  const Context& _context;
  EdgeWeight _round_delta;
  EdgeWeight _best_delta;
  size_t _num_fruitless_moves;
};

class AdaptiveStoppingRule {
public:
  AdaptiveStoppingRule(const Graph& graph, const Context&) :
    beta(std::log(graph.numNodes())) { }

  bool searchShouldStop() {
    return (numSteps > beta) && (Mk != 0 && numSteps >= ( variance / (Mk*Mk) ) * stopFactor );
  }

  void update(const EdgeWeight delta) {
    ++numSteps;
    if (numSteps == 1) {
      Mk = delta;
      MkPrevious = delta;
      SkPrevious = 0.0;
      variance = 0.0;
    } else {
      Mk = MkPrevious + (delta - MkPrevious) / numSteps;
      Sk = SkPrevious + (delta - MkPrevious) * (delta - Mk);
      variance = Sk / (numSteps - 1.0);

      MkPrevious = Mk;
      SkPrevious = Sk;
    }
  }

  void reset() {
    numSteps = 0;
    variance = 0.0;
  }

private:
  size_t numSteps = 0;
  double variance = 0.0, Mk = 0.0, MkPrevious = 0.0, Sk = 0.0, SkPrevious = 0.0;
  const double alpha = 4.0;   // make parameter if it doesn't work well
  const double stopFactor = (alpha / 2.0) - 0.25;
  double beta;
};

}  // namespace cluster_editing
