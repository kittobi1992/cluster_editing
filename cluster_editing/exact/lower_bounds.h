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


struct Triple {
    int u=-1,v=-1,w=-1;
    int cost=0;
    bool valid = false;

    Triple() = default;
    Triple(int a, int b, int c, const Edges& potential);
    // returns cost modification of total packing
    int apply(Edges& potential, bool undo=false);
    inline bool operator<(const Triple& rhs) const { return !rhs.valid || this->cost < rhs.cost; }
};

int packing_local_search_bound(const Instance& inst, int limit);
int meta_lower_bound(const Instance& inst, int limit);

std::vector<Triple> getAPacking(const Instance& inst);
