/*+++++++++++++++++ class LowerBound +++++++++++++++++++ */
#ifndef __lowerbound_h__
#define __lowerbound_h__

// include for IO for testing !!!
#include <sstream>

// includes from stlib
#include <vector>
#include <algorithm>

#include <costsgraph.h>
#include <array.h>
#include <trianglematrix.h>

/* This static class implements three different lower bounds.
   1.) calculate minimal solving costs per conflict triple
       of every edge 
   2.) calculate a set of independent conflict triples
   3.) like 2.) but includes the icf/icp costs of the
       most expensive edge

   1.) + 2.) need therefore just the graph itself. 3.) needs
   in addition the icf/icp costs in a matrix.

   The returned double value is the final lower bound.
*/

class LowerBound {
public:

    LowerBound() = delete;

    ~LowerBound() = delete;

    // structure to save a conflict triple and minimum solving value
    // float is used to avoid memory problems
    struct triple_type {
        unsigned short x;
        unsigned short y;
        unsigned short z;
        float value;

        triple_type(int _x, int _y, int _z, double _value) : x(_x), y(_y), z(_z), value(_value) {
            constexpr int MaxUShort = std::numeric_limits<unsigned short>::max();
            constexpr double MaxFloat = std::numeric_limits<float>::max();
            assert(_x < MaxUShort);
            assert(_y < MaxUShort);
            assert(_z < MaxUShort);
            assert(_value < MaxFloat);
        }

        friend bool operator<(const triple_type& lhs, const triple_type &rhs) {
            return lhs.value > rhs.value;
        }
    };

    // returns a lower bound by checking the minimal solving costs per conflict triple of an edge
    static inline double get(const CostsGraph &graph) {

        int size = graph.getSize();
        int nr_conflict_triples = 0;
        double min_costs_ratio = graph.permanent;

        std::vector<double> edge_list;

        // iterate over all edges
        for (int i = 0; i < size; i++) {
            for (int k = 0; k < i; k++) {

                // count the number of conflict triples this edge is in
                int ik_nr_confl_triples = 0;
                for (int h = 0; h < k; h++) {
                    int not_set = 0;
                    not_set += (graph.getEdge(i, k) <= 0) ? 1 : 0;
                    not_set += (graph.getEdge(h, k) <= 0) ? 1 : 0;
                    not_set += (graph.getEdge(i, h) <= 0) ? 1 : 0;
                    ik_nr_confl_triples += (not_set == 1) ? 1 : 0;
                }
                nr_conflict_triples += ik_nr_confl_triples;
                for (int h = k + 1; h < i; h++) {
                    int not_set = 0;
                    not_set += (graph.getEdge(i, k) <= 0) ? 1 : 0;
                    not_set += (graph.getEdge(h, k) <= 0) ? 1 : 0;
                    not_set += (graph.getEdge(i, h) <= 0) ? 1 : 0;
                    ik_nr_confl_triples += (not_set == 1) ? 1 : 0;
                }
                for (int h = i + 1; h < size; h++) {
                    int not_set = 0;
                    not_set += (graph.getEdge(i, k) <= 0) ? 1 : 0;
                    not_set += (graph.getEdge(h, k) <= 0) ? 1 : 0;
                    not_set += (graph.getEdge(i, h) <= 0) ? 1 : 0;
                    ik_nr_confl_triples += (not_set == 1) ? 1 : 0;
                }

                // calculate edge value divided by #conflict triples
                double ratio = (ik_nr_confl_triples > 0) ? (std::abs(graph.getEdge(i, k)) / ik_nr_confl_triples)
                                                         : graph.permanent;

                // insert edge with the ratio into the edge list as many times as it is occurs in conflict triples
                for (int h = 0; h < ik_nr_confl_triples; h++) {
                    edge_list.push_back(ratio);
                }

                // save minimum ratio (here not longer used...used for simple approximation)
                min_costs_ratio = std::min(min_costs_ratio, ratio);
            }
        }

        // sort edges by ratio
        sort(edge_list.begin(), edge_list.end());

        // sum first most expensive solving costs
        double sum_of_edges = 0;
        for (int i = 0; i < nr_conflict_triples; i++) {
            sum_of_edges += edge_list[i];
        }

        // return final lower bound
        return sum_of_edges;
    }


    // calculates lower bound by finding a set of independent conflict triples
    static inline double get2(const CostsGraph &graph) {
        // get independent conflict triples
        std::vector<Array<double>> final_list = getIndependentKT(graph);

        double bound = 0.0;

        // iterate through list
        for (const auto& triple : final_list) {
            // add minimum solving costs to lower bound
            bound += triple[0];
        }
        return bound;

    }


    // returns a maximal set (not maximum) of indpendent conflict triples, meaning they do not share an edge
    static inline std::vector<Array<double>> getIndependentKT(const CostsGraph &graph) {

        int size = graph.getSize();
        std::vector<triple_type> triple_list;

        // get all conflict triples
        for (int i = 0; i < size; i++) {
            for (int k = 0; k < i; k++) {
                for (int h = 0; h < k; h++) {
                    int not_set = 0;
                    double min;
                    not_set += (graph.getEdge(i, k) <= 0) ? 1 : 0;
                    min = std::abs(graph.getEdge(i, k));
                    not_set += (graph.getEdge(h, k) <= 0) ? 1 : 0;
                    min = std::min(min, std::abs(graph.getEdge(h, k)));
                    not_set += (graph.getEdge(i, h) <= 0) ? 1 : 0;
                    min = std::min(min, std::abs(graph.getEdge(i, h)));

                    // found conflict, save it array, first index=min solving costs
                    if (not_set == 1) {
                        triple_type triple = {i, k, h, min};
                        triple_list.push_back(triple);
                    }
                }
            }
        }

        // sort conflict triples by minimum solving costs
        sort(triple_list.begin(), triple_list.end());

        // create final list with independent conflict triples
        std::vector<Array<double>> final_list;

        // create matrix to reject conflict triples, initialized with false
        TriangleMatrix<bool> used_edges = TriangleMatrix<bool>(size, false);

        // iterate through list
        for (const auto& triple : triple_list) {
            unsigned short x = triple.x;
            unsigned short y = triple.y;
            unsigned short z = triple.z;

            // check if conflict triple was already rejected
            if (!used_edges.pos(x, y) &&
                !used_edges.pos(x, z) &&
                !used_edges.pos(y, z)) {

                // save conflict triple
                Array<double> triple_array = Array<double>(4);
                triple_array[0] = triple.value;
                triple_array[1] = x;
                triple_array[2] = y;
                triple_array[3] = z;
                final_list.push_back(triple_array);

                used_edges.pos(x, y) = true;
                used_edges.pos(x, z) = true;
                used_edges.pos(y, z) = true;
            }
        }

        return final_list;
    }


    // returns a lower bound by including icp/icf and independent conflict triples
    static inline double get3(const CostsGraph &graph, const TriangleMatrix<double> &min_insert_costs,
                              const TriangleMatrix<double> &min_delete_costs) {
        int size = graph.getSize();
        std::vector<triple_type> triple_list;
        Array<double> max_triple = Array<double>(4, 0.0);

        // iterate over all edges
        for (int i = 0; i < size; i++) {
            for (int k = 0; k < i; k++) {

                // check for conflict triples
                for (int h = 0; h < k; h++) {
                    int not_set = 0;
                    double min;
                    not_set += (graph.getEdge(i, k) <= 0) ? 1 : 0;
                    min = std::abs(graph.getEdge(i, k));
                    not_set += (graph.getEdge(h, k) <= 0) ? 1 : 0;
                    min = std::min(min, std::abs(graph.getEdge(h, k)));
                    not_set += (graph.getEdge(i, h) <= 0) ? 1 : 0;
                    min = std::min(min, std::abs(graph.getEdge(i, h)));

                    // found conclict triple
                    if (not_set == 1) {
                        // save conflict triple with the first index corresponding to the minimal solving costs
                        triple_type triple = {i, k, h, min};
                        triple_list.push_back(triple);

                        // check max icp/icf triple
                        double min_induced_costs = graph.permanent;
                        if (graph.getEdge(i, k) <= 0 && min_insert_costs.pos(i, k) < min_induced_costs) {
                            min_induced_costs = min_insert_costs.pos(i, k);
                        } else if (min_delete_costs.pos(i, k) < min_induced_costs) {
                            min_induced_costs = min_delete_costs.pos(i, k);
                        }
                        if (graph.getEdge(h, k) <= 0 && min_insert_costs.pos(h, k) < min_induced_costs) {
                            min_induced_costs = min_insert_costs.pos(h, k);
                        } else if (min_delete_costs.pos(h, k) < min_induced_costs) {
                            min_induced_costs = min_delete_costs.pos(h, k);
                        }
                        if (graph.getEdge(i, h) <= 0 && min_insert_costs.pos(i, h) < min_induced_costs) {
                            min_induced_costs = min_insert_costs.pos(i, h);
                        } else if (min_delete_costs.pos(i, h) < min_induced_costs) {
                            min_induced_costs = min_delete_costs.pos(i, h);
                        }

                        // save triple if icf/icp are maximal
                        if (max_triple[0] < min_induced_costs) {
                            max_triple[0] = min_induced_costs;
                            max_triple[1] = i;
                            max_triple[2] = k;
                            max_triple[3] = h;
                        }
                    }
                }
            }
        }

        // sort conflict triples by minimal solving costs
        sort(triple_list.begin(), triple_list.end());

        // now 2 lower bounds are calculated simultanously
        // the first with used_edges1 like in the previous lower bound function
        // the second with used_edges2 includes the icf/icp costs
        // std::vector<Array<double>> final_list;
        TriangleMatrix<bool> used_edges1 = TriangleMatrix<bool>(size, false);
        TriangleMatrix<bool> used_edges2 = TriangleMatrix<bool>(size, false);
        double bound1 = 0.0;
        double bound2;

        // first set max induced triple, therefore set all edges including one vertex
        // of the conflict triple to false
        //std::cout << "max triple (" << max_triple[1] << "," << max_triple[2] << "," << max_triple[3] << ") with " << max_triple[0] << std::endl;
        bound2 = max_triple[0];

        int max_triple_i = static_cast<int>(max_triple[1]);
        int max_triple_k = static_cast<int>(max_triple[2]);
        int max_triple_h = static_cast<int>(max_triple[3]);

        for (int i = 0; i < max_triple_i; i++) {
            used_edges2.pos(i, max_triple_i) = true;
        }
        for (int i = max_triple_i + 1; i < used_edges2.size(); i++) {
            used_edges2.pos(i, max_triple_i) = true;
        }
        for (int i = 0; i < max_triple_k; i++) {
            used_edges2.pos(i, max_triple_k) = true;
        }
        for (int i = max_triple_k + 1; i < used_edges2.size(); i++) {
            used_edges2.pos(i, max_triple_k) = true;
        }
        for (int i = 0; i < max_triple_h; i++) {
            used_edges2.pos(i, max_triple_h) = true;
        }
        for (int i = max_triple_h + 1; i < used_edges2.size(); i++) {
            used_edges2.pos(i, max_triple_h) = true;
        }


        // now insert to the lower bound all independent conflict triples
        for (const auto& triple : triple_list) {
            unsigned short x = triple.x;
            unsigned short y = triple.y;
            unsigned short z = triple.z;
            // check if conflict triple was already rejected
            if (!used_edges1.pos(x, y) && !used_edges1.pos(x, z) && !used_edges1.pos(y, z)) {
                // save conflict triple
                /* Array<double> triple_array = Array<double>(4);
                triple_array[0] = triple_list[i].value;
                triple_array[1] = x;
                triple_array[2] = y;
                triple_array[3] = z;
                final_list.push_back(triple_array); */
                bound1 += triple.value;
                //std::cout << "to b1 add (" << triple[1] << "," << triple[2] << "," << triple[3] << ") with " << triple[0] << std::endl;

                used_edges1.pos(x, y) = true;
                used_edges1.pos(x, z) = true;
                used_edges1.pos(y, z) = true;
            }
            if (!used_edges2.pos(x, y) && !used_edges2.pos(x, z) && !used_edges2.pos(y, z)) {
                // save conflict triple
                /* Array<double> triple_array = Array<double>(4);
                triple_array[0] = triple_list[i].value;
                triple_array[1] = x;
                triple_array[2] = y;
                triple_array[3] = z;
                final_list.push_back(triple_array); */
                bound2 += triple.value;
                //std::cout << "to b2 add (" << triple[1] << "," << triple[2] << "," << triple[3] << ") with " << triple[0] << std::endl;

                used_edges2.pos(x, y) = true;
                used_edges2.pos(x, z) = true;
                used_edges2.pos(y, z) = true;
            }

        }

        // the greater lower bound is then returned
        return std::max(bound1, bound2);
    }
};


#endif
