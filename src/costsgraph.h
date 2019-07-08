/*+++++++++++++++++ class CostsGraph +++++++++++++++++++ */
#ifndef __costsgraph_h__
#define __costsgraph_h__

// includes for IO
#include <iostream>
#include <sstream>
#include <cassert>

#include <vector>
#include <algorithm>

#include <vertexlists.h>
#include <array.h>
#include <trianglematrix.h>

#include <limits>

#include <Exceptions/graphexception.h>

// - costs graph saves a normal undirected weighted graph
// - used costs matrix is a linearized triangle matrix since its a undirected graph
// - positive costs implies the edge is set, negative or equal to zero are unset
// - its possible to delete vertices but not to add
// - vertices will not really be deleted, intern there are still in the adjacency matrix
//   but externally they are gone. therefore it uses a translation table
//   as the triangle matrix
// - names of vertices are always a list of original indices and a list of names
// - it also returns its connected components and searches for them
// - input and output is managed via function pointers, just turn in the right
//   input or output function an the object will use it

class CostsGraph {
public:
    /* ####  type definitions  #### */

    // vector to save many names per vertex = name list
    typedef std::vector<std::string> vertex_name_type;

    // vector to save indices of vertices = index list
    typedef std::vector<unsigned short> index_list_type;

    // information record of one node, with a list of names and
    // a list of original indices
    struct vertex_info_type {
        vertex_name_type names;
        index_list_type indices;
    };

    // vector to save all the vertex infos
    typedef std::vector<vertex_info_type> vertices_names_type;

    // saves an edge
    struct pair_type {
        unsigned short i;
        unsigned short j;

        pair_type() : i(0), j(0) {}

        pair_type(int _i, int _j) : i(_i), j(_j) {
            constexpr int MaxVal = std::numeric_limits<unsigned short>::max();
            assert(_i < MaxVal);
            assert(_j < MaxVal);
        }
    };

    // saves a list of edges
    typedef std::vector<pair_type> edge_list_type;

    // double vector will used for inititalization
    typedef std::vector<double> double_matrix_row_type;
    typedef std::vector<double_matrix_row_type> double_matrix_type;

    // double array will be used to get weight of edges if two vertices are merged
    typedef Array<double> double_array_type;

    // double triangle matrix
    typedef TriangleMatrix<double> triangle_matrix_type;

    // pos integer matrix and array for saving connected components
    typedef std::vector<unsigned short> vertex_set_type;
    typedef std::vector<vertex_set_type> vertex_matrix_type;

    // unsigned short for vector of indices for cliques
    typedef std::vector<unsigned short> byte_vector_type;

    // unsigned short array of translation table
    typedef Array<unsigned short> byte_array_type;

    // vector of references to graphs to return connected components
    typedef std::vector<CostsGraph *> graph_list_type;

    // function pointer to cost parsing function
    typedef double (*costs_parsing_fct_type)(double value, double threshold);

    // function pointer to double init function, a function filling the internal double matrix, by parsing
    // some values through the given cost function
    typedef void (*double_file_fct_type)(char *fname, CostsGraph &object, costs_parsing_fct_type fct);

    // function pointer for matrix io function
    typedef void (*matrix_file_fct_type)(char *fname, CostsGraph &G);

    // function pointer for matrix init function which needs still cost parsing
    typedef void (*matrix_file_fct_type2)(char *fname, CostsGraph &G, costs_parsing_fct_type fct, double threshold);

    // constants for different adjacency status
    const double permanent;
    const double forbidden;

    /* ####  Constructors  #### */
    // th = threshold, every value smaller then th will be set to set, else to not_set
    explicit CostsGraph(int size = 1, double th = 10e-20);

    CostsGraph(int size, std::string *edge_list, double th = 10e-20);

    CostsGraph(int size, vertex_name_type edge_list, double th = 10e-20);

    CostsGraph(int size, double_matrix_type costs, double th = 10e-20);

    CostsGraph(int size, double_matrix_type costs, vertex_name_type edge_list, double th = 10e-20);

    CostsGraph(int size, double_matrix_type weight_matrix, costs_parsing_fct_type cost_fct, double th = 10e-20);

    CostsGraph(int size, double_matrix_type weight_matrix, costs_parsing_fct_type cost_fct, vertex_name_type edge_list,
               double th = 10e-20);

    // CostsGraph(char *file_name, double_file_fct_type fct, costs_parsing_fct_type cost_fct, double th = 10e-20);

    CostsGraph(const std::string &file_name, double_file_fct_type fct, costs_parsing_fct_type cost_fct, double th = 10e-20);

    // CostsGraph(char *file_name, matrix_file_fct_type2 fct, costs_parsing_fct_type cost_fct, double th = 10e-20);

    CostsGraph(const std::string &file_name, matrix_file_fct_type2 fct, costs_parsing_fct_type cost_fct, double th = 10e-20);

    // CostsGraph(char *file_name, matrix_file_fct_type fct);

    CostsGraph(const std::string &file_name, matrix_file_fct_type fct);

    CostsGraph(CostsGraph const &new_graph);


    /* ####  Destructor  #### */
    ~CostsGraph() = default;


    /* ####  access functions  #### */

    // set edge to given adjacency status
    void setEdge(int i, int j, double value);

    // return costs of edge
    inline double getEdge(int i, int j) const { return _cost_matrix.pos(i, j); };

    // get next neighbor of node, starting at k
    // returns -1 for no neighbor found
    short getNeighbor(int i, int k) const;

    // sets edge weight by using the cost-parsing function
    void setEdgeWeight(int i, int j, double value, costs_parsing_fct_type);

    // returns complete name of vertex divided by |
    std::string getVertexName(int i) const;

    // set first name to value
    void setVertexName(int i, std::string name);

    // return complete vertex info of node
    vertex_info_type getVertexInfo(int i) const;

    // set complete vertex info of a node
    void setVertexInfo(int i, vertex_info_type info);

    // convert the names of all nodes to a string, each node gets a numbering
    std::string vertexNamesToString() const;

    // returns index list of a node
    index_list_type getVertexIndices(int i) const;

    // returns name list of a node
    vertex_name_type getVertexNames(int i) const;

    // returns graph size
    int getSize() const;

    // return index of vertex with this name
    int getIndex(std::string key) const;

    // checks if vertex name is already defined
    bool inGraph(std::string name) const;

    // returns true if any of the names in vector is in the graph defined
    bool inGraph(vertex_name_type names) const;

    // calculates degree of vertex
    int getDegree(int i) const;

    // returns number of set edges
    int getEdgeNumber() const;

    // returns threshold
    double getThreshold() const;

    // returns true of i is member of clique, false if not
    bool inClique(int i) const;

    // returns vector of clique members, empty if i is not in a clique
    byte_vector_type getClique(int i) const;

    /* ####  edit functions  #### */
    // deletes vertex from graph
    void deleteVertex(int index);

    // deletes vertex from graph
    void deleteVertices(byte_vector_type set);

    // merges to vertices and gives them the given weight vector as weights to the other vertices
    void mergeVertices(int index1, int index2, double_array_type costs);

    // deletes clique if vertex 'index' is element of clique
    bool deleteClique(int index);

    // calculate connected components and save them in a matrix
    vertex_matrix_type getConnectedVertices();

    // divide the graph into components, each again a costs graph
    graph_list_type getConnectedComponents();

    // dependend on the given connected components, saved in the vertex matrix
    // the graph will be divided into costs graph again
    graph_list_type getConnectedComponents(vertex_matrix_type vertex_matrix);

    /* ####  IO-functions  #### */
    // creates string from costs matrix
    std::string toString(bool simple) const;

    // do IO operation with given function
    void costsIOOperation(char *fname, matrix_file_fct_type fct);

    void costsIOOperation(std::string fname, matrix_file_fct_type fct);

    /* ####  overloaded operators  #### */
    explicit operator std::string();

    friend std::ostream &operator<<(std::ostream &output, const CostsGraph &g);

    CostsGraph &operator=(const CostsGraph &right);

private:
    /* ####  member variables  #### */

    // names of vertices
    vertices_names_type _info_list;

    // internal costs triangle matrix for the actual graph
    triangle_matrix_type _cost_matrix;

    // threshold if we parse in similarity data
    double _threshold;

    // size of graph
    int _size;


    /* ####  help functions  #### */

    //inits adjacency matrices and help arrays
    void initMatrices();

    void initMatrices(char *fname, double_file_fct_type fct, costs_parsing_fct_type cost_fct);

    void initMatrices(char *fname, matrix_file_fct_type2 fct, costs_parsing_fct_type cost_fct);

    // translate index to internal index
    inline int getIntern(int i) const;

    // conversion functions
    template<typename T>
    inline std::string valueToString(T x) const;

    // prints costs matrix to screen
    std::string matrixToString(bool simple) const;
};

#endif
