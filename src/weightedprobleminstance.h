/*++++++++++++++++ class WeightedProblemInstance +++++++++++++++++++ */

#ifndef __weightedprobleminstance_h__
#define __weightedprobleminstance_h__

// include for IO for testing !!!
#include <sstream>

// includes from stlib
#include <vector>
#include <math.h>
#include <algorithm>

#include <edgereduction.h>
#include <mergereduction.h>
#include <bnmanager.h>
#include <costsgraph.h>
#include <lowerbound.h>
#include <Exceptions/probleminstanceexception.h>
#include <almostclique.h>
#include <criticalclique.h>
#include <cc-kernel.h>

#include <array.h>
#include <trianglematrix.h>

// - implements a weighted problem instance object
// - it implements simple reduction rules as the major reduction cycle
// - it contains a graph and a parameter
// - icf/icp reduction is done by using the edgeredction object
// - merging is done using the mergereduction object
// - almost clique and critical clique are done with the correspondend object

// - the three main functions are reduce, strongReduce and maxReduce
//   with different strong reduction routines
// - these functions will start the inner reduction circle which works somehow like that:
//   - reduction rules are called which return lists of edges which need to be reduced
//   - these edges are passed to a function which sets them to the appropriate value
//   - in each set operation the edgereduction object is updated and maybe the merging rule is called

// - further the objects provides a lower bound for itself, necessary in a search tree
// - it can also return the edge with minimal branching number, necessary in a search tree
// - it is able to split in different objects if the inner graph splits into connected components
// - all changes are saved in a list can be ordered via a function

class WeightedProblemInstance {
public:
    /* ####  Constants  #### */

    // constants for different adjacency status
    // note that any other adjacency values should be below/above the given values
    const double permanent;
    const double forbidden;


    /* ####  type definitions  #### */

    // edge struct to save an edge and its changing costs
    // it also saves the original indices of this edge, since merging could change the indices
    struct edge_type {
        CostsGraph::index_list_type i;
        CostsGraph::index_list_type j;
        double costs;
        char operation;
    };

    // vector to save a list of edge changes
    typedef std::vector<edge_type> edge_list_type;

    // pair type (=edge) for simple reduction rules
    typedef CostsGraph::pair_type pair_type;

    // for new_costs array in case of merging
    typedef CostsGraph::double_array_type double_array_type;

    // list of problem instance objects to divide objects whenever there are different connected components
    typedef std::vector<WeightedProblemInstance *> pi_list_type;

    /* ####  Constructors  #### */
    WeightedProblemInstance(CostsGraph graph, double parameter, bool para_independent = false);

    WeightedProblemInstance(const WeightedProblemInstance &PI);

    WeightedProblemInstance(CostsGraph &graph, double parameter, double start_parameter,
                            const EdgeReduction &edge_reduction,
                            bool para_independent = false);

    /* ####  Destructor  #### */
    ~WeightedProblemInstance() = default;

    /* ####  main function reduce  #### */

    // reduces the problem kernel, by setting edges to permanent/forbidden, since
    // other operations will cause costs bigger then the parameter
    int reduce();

    // stronger reduce which includes more involved reduction rules
    int strongReduce();

    // maximal reduce which includes more involved reduction rules and reduces as long as possible
    int maxReduce();

    // unweighted critical clique merging
    int ccReduce();

    // unweighted critical clique reduction after Guo to achieve a 4k kernel in the unweighted case
    int ccKernelization();

    // use heuristic to solve instance...simply by iterativley setting edges to permanent and forbidden
    // dependend on their icp and icf costs
    int heuristicSolve();

    // set edge to values in the underlying graph and do every other necessary update
    void setEdge(int i, int j, double value);

    // return edge value of the underlying graph
    inline double getEdge(short i, short j) const { return _graph.getEdge(i, j); };

    // returns changes in the graph made by the object
    inline edge_list_type getChangedEdges() const { return _changed_edges; }

    // saves new changes, if it makes sence or not
    // inline void setChangedEdges(edge_list_type new_edges) { _changed_edges = std::move(new_edges); }

    // returns parameter
    inline double getParameter() const { return _parameter; };

    // return start parameter + 1
    // gives the border for costs which exceed the parameter
    inline double getStartParameter() const { return _start_parameter; };

    // returns a copy of the graph
    inline CostsGraph getGraph() const { return _graph; };

    inline double getMinInsertCosts(int i, int j) { return _edge_reduction.getMinInsertCosts(i, j); };

    // return lower bound for instance
    double getLowerBound() const;

    // return a upper bound for this instance...pidr needs to be default=false
    double getUpperBound() const;

    // make a deep copy of the object
    WeightedProblemInstance *clone() const;

    // return edge to branch for, -1 if there exists no edge
    double getEdgeForBranching(int &i, int &j) const;

    // functions to divide and merge instances
    pi_list_type divideInstance(const CostsGraph::vertex_matrix_type &vertex_matrix);

    // sets parameter and therefore changes the problem kernel
    void setParameter(double _new_parameter);

private:
    /* ####  additional member variables  #### */
    EdgeReduction _edge_reduction;

    // reducable parameter for problem instance
    double _parameter;

    // original parameter
    double _start_parameter;

    // option to switch of parameter-dependent red rules
    bool _pid;

    // underlying costs graph
    CostsGraph _graph;

    // contains all changed edges
    edge_list_type _changed_edges;


    /* ####  general help functions  #### */

    // sets edge to certain value and consider updates for number of (non)-common neighbors
    inline void setEdgeToValue(int i, int j, double value);

    // save changed edge in list
    inline void saveChangedEdge(int i, int j, double costs, char operation);

    // merge two vertices and all other structures like tripel lists and min insert & delete lists
    inline void mergeVertices(int i, int j);

    // deletes created clique from linked lists
    inline void mergeVertices(CostsGraph::byte_vector_type &clique);

    inline void deleteClique(CostsGraph::byte_vector_type &clique);

    // inline void deleteClique(int i);

    // set a set of edges to forbidden or permanent respectivley
    inline void setEdgesToForbidden(const CostsGraph::edge_list_type &forbiddden_list);

    inline void setEdgesToPermanent(CostsGraph::edge_list_type permanent_list);

    // try simple reduction rules
    inline int simpleReduction();

    // do the approximation of the critical clique as simple reduction rule
    inline int criticalCliqueApproximation();

    // call the more involved weighted critical clique rule which uses dynamic programming
    inline int criticalCliqueReduction(int max_edge, bool check_for_int = false);

    // use the almost clique rule to detect cliques
    inline int detectCliques();

};

#endif
