/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include <atomic>

namespace cluster_editing {

class SpinLock {
public:
  // boilerplate to make it 'copyable'. but we just clear the spinlock. there is never a use case to copy a locked spinlock
  SpinLock() { }
  SpinLock(const SpinLock&) { }
  SpinLock& operator=(const SpinLock&) { spinner.clear(std::memory_order_relaxed); return *this; }

  bool tryLock() {
    return !spinner.test_and_set(std::memory_order_acquire);
  }

  void lock() {
    while (spinner.test_and_set(std::memory_order_acquire)) {
      // spin
      // stack overflow says adding 'cpu_relax' instruction may improve performance
    }
  }

  void unlock() {
    spinner.clear(std::memory_order_release);
  }

private:
  std::atomic_flag spinner = ATOMIC_FLAG_INIT;
};

} // namespace cluster_editing