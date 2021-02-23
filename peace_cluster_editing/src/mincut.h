/*+++++++++++++++++ class MinCut +++++++++++++++++++ */
#ifndef __mincut_h__
#define __mincut_h__

// include for IO for testing !!!
#include <sstream>

// includes from stlib
#include <vector>
#include <algorithm>
#include <math.h>

#include <costsgraph.h>
#include <array.h>

/* This static class implements an s-t min-cut algorithm
   from Ford & Fulkerson and a graph min-cut algorithm
   from Frank Wagner. Google for clearance.
*/

class MinCut {
public:

    MinCut() = delete;

    ~MinCut() = delete;

    /* ################ Graph Min Cut ############################### */

    // return minimum s-t cut produced by a ford-fulkerson algorithm
    static inline double getMinSTCut(int s, int t, CostsGraph graph) {
        double min_cut = 0.0;

        // check egde itself and delete it
        if (graph.getEdge(s, t) > 0) {
            min_cut += graph.getEdge(s, t);
            graph.setEdge(s, t, 0.0);
        }

        // check all common neighbors first to avoide useless running time and delete them if necessary
        for (int i = 0; i < graph.getSize(); i++) {
            if (i != s && i != t) {
                double edge1 = graph.getEdge(s, i);
                double edge2 = graph.getEdge(t, i);
                if (edge1 > 0 && edge2 > 0) {
                    if (edge1 < edge2) {
                        min_cut += edge1;
                        graph.setEdge(s, i, 0.0);
                    } else {
                        min_cut += edge2;
                        graph.setEdge(t, i, 0.0);
                    }
                }
            }
        }

        // now check for the maximal edge value
        // which is later our minimal edge we take into account for our calculation
        double min_value = 0.0;
        for (int i = 0; i < graph.getSize(); i++) {
            for (int k = 0; k < i; k++) {
                min_value = (graph.getEdge(i, k) > min_value && graph.getEdge(i, k) != graph.permanent) ? graph.getEdge(
                        i, k) : min_value;
            }
        }

        // as long there min_value is large enough...
        while (min_value > 1) {

            // check for augmenting paths, where we only take edges > min_value into account
            // and delete them
            while (augmentingPath(graph, s, t, min_cut, min_value)) {
            }

            // decrease min_value
            min_value = min_value / 2.0;
        }

        // do it again without any lower bound for edges
        while (augmentingPath(graph, s, t, min_cut, 0)) {
        }

        // returned so obtained min-cut
        return min_cut;
    }

    /* ################ Graph Min Cut ############################### */


    // returns graph min-cut as Frank Wagner algorithm suggested
    inline static double getGraphMinCut(CostsGraph clique_graph) {
        int size = clique_graph.getSize();
        double min_cut = clique_graph.permanent;

        // iterate so many time as the graph has vertices
        for (int i = 0; i < size - 1; i++) {
            // get an s-t cut and merge the last two vertices in this function
            double new_cut = getSTCut(clique_graph);

            // get the minimum cut of all so obtained cuts
            min_cut = (new_cut < min_cut) ? new_cut : min_cut;
        }

        return min_cut;
    }


private:

    /* ############ help functions for S-T min-cut algorithm ############### */


    // search for augmenting path from s to to considering only edges greater than min_cut
    // and delete this path
    static inline bool augmentingPath(CostsGraph &graph, int s, int t, double &min_cut, double min_value) {
        // create vector to hold path
        std::vector<int> path = std::vector<int>(0);

        // create vector to save already vistited vertices in the recursion
        std::vector<bool> used_vertices = std::vector<bool>(graph.getSize(), false);

        // search for a path
        if (seekPath(graph, s, t, path, used_vertices, min_value)) {

            // search for the minimal edge in this path
            double min = graph.permanent;
            for (int i = 0; i < path.size() - 1; i++) {
                if (graph.getEdge(path[i], path[i + 1]) < min) {
                    min = graph.getEdge(path[i], path[i + 1]);
                }
            }

            // decrease all edge along the path
            for (int i = 0; i < path.size() - 1; i++) {
                double edge = graph.getEdge(path[i], path[i + 1]);
                graph.setEdge(path[i], path[i + 1], edge - min);
            }

            // increase min_cut
            min_cut += min;
            return true;
        } else {
            return false;
        }
    }


    // recursive function to search an s-t path
    static inline bool
    seekPath(const CostsGraph &graph, int s, int t, std::vector<int> &path, std::vector<bool> &used_vertices,
             double min_value) {
        // abort if s == t
        if (s == t) {
            path.insert(path.end(), s);
            return true;
        }

        // otherwise mark vertex
        used_vertices[s] = true;

        // ceck for all neighbors
        for (int i = 0; i < graph.getSize(); i++) {

            // if edge to neighbor is considered, since min_value < edge-value
            // call recursion
            if (s != i && graph.getEdge(s, i) > min_value && !used_vertices[i]) {
                if (seekPath(graph, i, t, path, used_vertices, min_value)) {
                    // by going up the final path is filled
                    path.insert(path.end(), s);
                    return true;
                }
            }
        }

        return false;
    }

    /* ######################### help function for graph min cut ###################### */

    // get an s-t cut of the remaining graph
    inline static double getSTCut(CostsGraph &clique_graph) {
        // list to save vertices so far in the set
        std::vector<int> vertex_set = std::vector<int>(1, 0);

        double st_cut = 0.0;

        int last_max_vertex = 0;
        int current_max_vertex = 0;

        // now search always for the cheapest way to extend the vertex set
        // add this vertex, until no vertex is left
        while (vertex_set.size() < clique_graph.getSize()) {

            int set_member = 0;
            int max_vertex = 0;
            double max = 0;

            // iterate over all vertices not in vertex_set
            // check which vertex has the strongest connection
            // to the vertices within the set
            for (int i = 0; i < clique_graph.getSize(); i++) {
                if (set_member < vertex_set.size() && vertex_set[set_member] == i) {
                    set_member++;
                } else {
                    // check edge sum of this neighbor to vertex_set
                    double edge_sum = 0.0;
                    for (int vertex : vertex_set) {
                        edge_sum += std::max(0.0, clique_graph.getEdge(i, vertex));
                    }
                    // check wether this is the max vertex
                    if (edge_sum > max) {
                        max_vertex = i;
                        max = edge_sum;
                    }
                }
            }

            // add max vertex to set
            vertex_set.push_back(max_vertex);
            sort(vertex_set.begin(), vertex_set.end());

            current_max_vertex = max_vertex;

            // save if necessary the last cut
            if (vertex_set.size() == clique_graph.getSize()) {
                st_cut = max;
            } else {
                last_max_vertex = max_vertex;
            }
        }

        // if graph is not trivial then merge last two vertices
        if (clique_graph.getSize() > 2) {

            // merge last vertex and the one before
            CostsGraph::double_array_type new_costs(clique_graph.getSize() - 2);

            // swap to be faster
            if (last_max_vertex > current_max_vertex) {
                std::swap(last_max_vertex, current_max_vertex);
            }

            // get new costs column
            for (int i = 0; i < last_max_vertex; i++) {
                new_costs[i] = fabs(clique_graph.getEdge(last_max_vertex, i)) +
                               fabs(clique_graph.getEdge(current_max_vertex, i));
            }
            for (int i = last_max_vertex + 1; i < current_max_vertex; i++) {
                new_costs[i - 1] = fabs(clique_graph.getEdge(last_max_vertex, i)) +
                                   fabs(clique_graph.getEdge(current_max_vertex, i));
            }
            for (int i = current_max_vertex + 1; i < clique_graph.getSize(); i++) {
                new_costs[i - 2] = fabs(clique_graph.getEdge(last_max_vertex, i)) +
                                   fabs(clique_graph.getEdge(current_max_vertex, i));
            }

            // merge vertices
            clique_graph.mergeVertices(last_max_vertex, current_max_vertex, new_costs);
        }

        return st_cut;
    }

};


#endif
