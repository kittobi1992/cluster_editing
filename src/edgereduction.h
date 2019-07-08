/*+++++++++++++++++ class EdgeReduction +++++++++++++++++++ */
#ifndef __edgereduction_h__
#define __edgereduction_h__

// include for IO for testing !!!
#include <sstream>

// includes from stlib
#include <vector>
#include <math.h>

#include <costsgraph.h>
#include <Exceptions/probleminstanceexception.h>
#include <mincut.h>
#include <lowerbound.h>
#include <algorithm>

#include <array.h>
#include <trianglematrix.h>

/* EdgeReduction implements the calculation of icp and icf costs
   (see paper APBC or RECOMB of weighted cluster editing).

  The current icp and icf costs are saved in a triangle matrix for
  every edge. The matrix is initalized by a constructor or by the init
  function.

  Whenever an edge is set to forbidden or permanent function
  setEdgeToValue is called. This leads to the update of these matrices
  in an appropriate way.

*/

class EdgeReduction {
public:
    /* ####  type definitions  #### */

    // vector for edge list indices, when costs will have to be translated into a list number
    typedef std::vector<unsigned short> byte_array_type;

    // matrix for min delete/insert costs
    typedef TriangleMatrix<double> triangle_matrix_type;

    // array of new adjacency values of vertices were merged
    typedef CostsGraph::double_array_type double_array_type;

    // list of pointer to object itself, necessary if objects gets split
    typedef Array<EdgeReduction *> edge_red_list_type;


    /* ####  Constructors  #### */

    // empty constructor
    EdgeReduction() : _graph(NULL), _parameter(0.0) {};

    // initate object with a graph
    EdgeReduction(CostsGraph *graph, double parameter);

    // copy constructor
    EdgeReduction(const EdgeReduction &ER, CostsGraph *graph);

    // copy constructor if elements gets split
    inline EdgeReduction(CostsGraph *graph, double parameter, triangle_matrix_type m_i_c, triangle_matrix_type m_d_c);


    /* ####  Destructor  #### */
    ~EdgeReduction();


    /* ####  main function reduce  #### */

    // make a deep copy of the object
    inline EdgeReduction *clone() const;

    // initate object after it was initalized with empty constructor
    void init(CostsGraph *graph, double parameter, bool merge = false);

    // sets edge to certain value and consider updates for number of (non)-common neighbors
    void setEdgeToValue(int i, int j, double old_value, double parameter, double parameter_before);

    // merge two vertices and all other structures like min insert & delete matrix
    void mergeVertices(int i, int j, double parameter, double_array_type new_costs);

    // return edge to reduce
    bool getReduceableEdges(CostsGraph::edge_list_type &permanent_edges, CostsGraph::edge_list_type &forbidden_edges,
                            double fuzzy_threshold);

    // functions to divide instances
    edge_red_list_type
    divideInstance(CostsGraph::vertex_matrix_type vertex_matrix, CostsGraph::graph_list_type graph_list);

    // sets parameter and therefore changes the problem kernel
    void setParameter(double _new_parameter);

    // return min insert or delete costs
    double getMinDeleteCosts(unsigned short i, unsigned short j) const;

    double getMinInsertCosts(unsigned short i, unsigned short j) const;

    // return min insert and min delete matrix, used for example for lower bound calculation
    triangle_matrix_type getMinInsertMatrix() const;

    triangle_matrix_type getMinDeleteMatrix() const;

    // assign operator
    EdgeReduction &operator=(const EdgeReduction &er);

    /* #### heuristic functions #### */

    // heuristic functions...return the edge with highest icp/icf value
    CostsGraph::edge_list_type getTopEdgesToDelete(int nr_edges);

    CostsGraph::edge_list_type getTopEdgesToInsert(int nr_edges);

    /* ####  help functions - just for testing !!! #### */
    void printLists();


private:
    /* ####  internal types #### */

    // type to save an edge in the heuristic part
    struct edge_type {
        int i;
        int j;
        double value;
    };

    /* ####  additional member variables  #### */

    // costs matrices
    triangle_matrix_type _min_insert_costs;
    triangle_matrix_type _min_delete_costs;

    // pointer to graph
    CostsGraph *_graph;

    // define forbidden & permanent
    double forbidden;
    double permanent;

    // define parameter
    double _parameter;

    /* ####  constructor help functions  #### */

    // calculate costs
    inline void createCostsMatrices();

    // creates to-proceed lists
    inline void calculateCostsForEdge(int i, int k, double &sum_insert, double &sum_delete);

    inline void createCostsLists();

    inline void getInDelCosts(int i, int k, int n, double &sum_insert, double &sum_delete);

    inline void getInDelCosts(double e1, double e2, double &sum_insert, double &sum_delete);

    /* ####  general help functions  #### */

    // converting array for converting costs to lists in O(1)
    inline double minDelEdge(double a, double b);

    inline double minInDelEdge(double a, double b);

    inline double setMinDeleteCosts(unsigned short i, unsigned short j, double value);

    inline double setMinInsertCosts(unsigned short i, unsigned short j, double value);

    // copies edge from one list to another considering the update in the reference matrix
    //void copyEdgeAndUpdateReference(int i, int j, char type, int from, int to);
    // updates 'n' or 'c' for (non)-common neighbors by value add
    inline void updateEdgeReference(int i, int j, char type, double value);

    // update costs matrices after merge operation
    inline void updateEdgeCase1(int i, int j, int a, int b, double new_costs1, double new_costs2, double &up_insert,
                                double &up_delete);

    inline void updateEdgeCase2(int i, int j, int a, int b, double new_costs1, double new_costs2, double &up_insert,
                                double &up_delete);

    inline void updateCostsMatrices(int i, int j, double_array_type new_costs);

    // check for all edges whether we need to insert or delete them
    inline CostsGraph::edge_list_type
    getEdgesToInsert(double min_cut_threshold, triangle_matrix_type &lower_bound_edge);

    inline CostsGraph::edge_list_type getEdgesToDelete(bool strong, triangle_matrix_type &lower_bound_edge);

    // fill lower bound matrix
    triangle_matrix_type createLowerBoundMatrix();

    // compare function for sorting edges...used by heuristic
    static inline bool compare(edge_type t1, edge_type t2) {
        return (t1.value > t2.value);
    }

};


#endif
