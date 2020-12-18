/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2016 Sebastian Schlag <sebastian.schlag@kit.edu>
 * Copyright (C) 2016 Yaroslav Akhremtsev <yaroslav.akhremtsev@kit.edu>
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

#include <algorithm>
#include <cstddef>
#include <random>
#include <type_traits>
#include <utility>
#include <vector>

namespace cluster_editing {
namespace math {
// See: Warren, Henry S. Hacker's Delight (Second Edition), p.61
template <typename T>
static inline T nextPowerOfTwoCeiled(T n) {
  static_assert(std::is_integral<T>::value, "Integer required.");
  static_assert(std::is_unsigned<T>::value, "Incompatible type");
  --n;
  for (size_t k = 1; k != 8 * sizeof(n); k <<= 1) {
    n |= n >> k;
  }
  return n + 1;
}

template <typename T, typename U>
static inline T nearestMultipleOf(T num, U multiple) {
  return (num + multiple - 1) & ~(multiple - 1);
}

template <typename T>
inline double median(const std::vector<T>& vec) {
  double median = 0.0;
  if ((vec.size() % 2) == 0) {
    median = static_cast<double>((vec[vec.size() / 2] + vec[(vec.size() / 2) - 1])) / 2.0;
  } else {
    median = vec[vec.size() / 2];
  }
  return median;
}

// based on: http://mathalope.co.uk/2014/07/18/accelerated-c-solution-to-exercise-3-2/
template <typename T>
inline std::pair<double, double> firstAndThirdQuartile(const std::vector<T>& vec) {
  if (vec.size() > 1) {
    const size_t size_mod_4 = vec.size() % 4;
    const size_t M = vec.size() / 2;
    const size_t ML = M / 2;
    const size_t MU = M + ML;
    double first_quartile = 0.0;
    double third_quartile = 0.0;
    if (size_mod_4 == 0 || size_mod_4 == 1) {
      first_quartile = (vec[ML] + vec[ML - 1]) / 2;
      third_quartile = (vec[MU] + vec[MU - 1]) / 2;
    } else if (size_mod_4 == 2 || size_mod_4 == 3) {
      first_quartile = vec[ML];
      third_quartile = vec[MU];
    }
    return std::make_pair(first_quartile, third_quartile);
  } else {
    return std::make_pair(0.0, 0.0);
  }
}


// see: http://graphics.stanford.edu/~seander/bithacks.html#IntegerLog10
static const uint64_t powers_of_10[] = {
  0,
  10,
  100,
  1000,
  10000,
  100000,
  1000000,
  10000000,
  100000000,
  1000000000,
  10000000000,
  100000000000,
  1000000000000,
  10000000000000,
  100000000000000,
  1000000000000000,
  10000000000000000,
  100000000000000000,
  1000000000000000000,
  10000000000000000000U
};

static inline uint8_t digits(const uint64_t x) {
  uint64_t t = (64 - __builtin_clzll(x | 1)) * 1233 >> 12;
  return static_cast<uint8_t>(t - (x < powers_of_10[t]) + 1);
}

}  // namespace math
}  // namespace cluster_editing
