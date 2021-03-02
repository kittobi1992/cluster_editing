
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
