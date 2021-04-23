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

#include "cluster_editing/context/context.h"
#include "cluster_editing/utils/randomize.h"

namespace cluster_editing {

template<typename Action>
class ActionSelector {

  static constexpr size_t HISTORY_WINDOW_SIZE = 5;

  struct ActionStats {
    Action action;
    float prob;
    std::vector<EdgeWeight> improvement_history;
  };

 public:
  explicit ActionSelector() :
    _actions() { }

  explicit ActionSelector(const std::vector<Action>& actions) :
    _actions() {
    ASSERT(actions.size() > 0UL);
    initializeActions(actions);
  }

  ActionSelector(const ActionSelector&) = delete;
  ActionSelector(ActionSelector&&) = delete;

  ActionSelector & operator= (const ActionSelector &) = delete;
  ActionSelector & operator= (ActionSelector &&) = delete;

  Action chooseAction(const bool verbose) const {
    ASSERT(_actions.size() > 0);
    if ( verbose ) {
      LOG << "Action Probabilities:";
      for ( size_t i = 0; i < _actions.size(); ++i ) {
        LOG << _actions[i].action << "=" << _actions[i].prob;
      }
    }
    const float p = utils::Randomize::instance().getRandomFloat(0.0, 1.0);
    float prefix_sum = 0.0;
    for ( size_t i = 0; i < _actions.size(); ++i ) {
      if ( p >= prefix_sum && p <= prefix_sum + _actions[i].prob ) {
        return _actions[i].action;
      }
      prefix_sum += _actions[i].prob;
    }
    return _actions[0].action;
  }

  void notifyImprovement(const Action action, const EdgeWeight delta) {
    for ( size_t i = 0; i < _actions.size(); ++i ) {
      if ( action == _actions[i].action ) {
        _actions[i].improvement_history.push_back(delta);
        break;
      }
    }
    updateProbs();
  }

  void initializeActions(const std::vector<Action>& actions) {
    _actions.clear();
    _actions.resize(actions.size());
    for ( size_t i = 0; i < actions.size(); ++i ) {
      _actions[i].action = actions[i];
      _actions[i].prob = 1.0 / static_cast<double>(actions.size());
    }
  }

 private:
  void updateProbs() {
    bool recompute_probs = true;
    std::vector<EdgeWeight> window_delta(_actions.size(), 0);
    EdgeWeight total_delta = 0;

    for ( size_t i = 0; i < _actions.size(); ++i ) {
      const size_t history_size = _actions[i].improvement_history.size();
      if ( history_size >= HISTORY_WINDOW_SIZE ) {
        for ( size_t j = history_size - HISTORY_WINDOW_SIZE; j < history_size; ++j ) {
          window_delta[i] += _actions[i].improvement_history[j];
          total_delta += _actions[i].improvement_history[j];
        }
        if ( window_delta[i] == 0 ) {
          ++window_delta[i];
          ++total_delta;
        }
      } else {
        // not enough entry in history, wait for more stats...
        recompute_probs = false;
        break;
      }
    }

    if ( recompute_probs ) {
      // Compute new probabilities based on improvement history
      for ( size_t i = 0; i < _actions.size(); ++i ) {
        _actions[i].prob = static_cast<double>(window_delta[i]) / total_delta;
      }
    }
  }

  std::vector<ActionStats> _actions;
};
}  // namespace cluster_editing
