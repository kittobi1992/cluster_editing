/*+++++++++++++++++ class GraphSet +++++++++++++++++++ */
#ifndef __graphset_h__
#define __graphset_h__

#include <vector>

#include <list>

#include <vertexlists.h>

#include <sys/types.h>
#include <dirent.h>

#include <costsgraph.h>

// - simple container for costs graph objects
// - saves only pointer to the objects
// - has the ability to read a directory of graphs if a graph-io-function
//   is indicated
// - if the graphset is constructed using a graph, then the graph
//   will be checked for connected components and each of it
//   is saved as a graph within the set

class GraphSet {
public:
    // list of pointers to graph objects
    typedef std::vector<CostsGraph *> graph_list_type;

    // function pointer to cost parsing function
    typedef double (*costs_parsing_fct_type)(double blast_value, double threshold);

    // function pointer to double IO function
    typedef void (*graph_set_parser_fct_type)(char *fname, GraphSet &graph_set, double th, costs_parsing_fct_type fct);

    // function pointer to matrix init function
    typedef CostsGraph::matrix_file_fct_type matrix_file_fct_type;

    // function pointer to weight matrix initializer function using a cost parsing function
    typedef CostsGraph::matrix_file_fct_type2 matrix_file_fct_type2;


    /* ####  Constructor  #### */

    // emtpy constructor
    GraphSet();

    // copy constructor
    GraphSet(const GraphSet &graph_set);

    // gets graph and will save every connected component in the graph as a graph
    explicit GraphSet(const CostsGraph &graph);

    // inits a graph from file, function can be choosen to seperate graph in connected components
    // GraphSet(char *fname, double th, costs_parsing_fct_type fct, graph_set_parser_fct_type parser_fct);

    GraphSet(const std::string &fname, double th, costs_parsing_fct_type fct, graph_set_parser_fct_type parser_fct);

    // reads in a weight graph and using the appropriate costs parser to create costs graph(s) out of it
    // GraphSet(char *file_name, matrix_file_fct_type2 fct, costs_parsing_fct_type cost_fct, double th);

    GraphSet(const std::string &file_name, matrix_file_fct_type2 fct, costs_parsing_fct_type cost_fct, double th);

    // reads a directory full of graphs in
    // GraphSet(char *dir_name, matrix_file_fct_type fct);

    GraphSet(const std::string &dir_name, matrix_file_fct_type fct);

    /* ####  Destructor  #### */
    ~GraphSet();


    /* ####  access functions  #### */

    // returns number of graphs
    int getSetSize();

    // returns reference to graph
    CostsGraph &getGraph(int i);

    // deletes graph from set
    void deleteGraph(int i);

    // add another graph, gets reference and will save pointer to deep copy of graph
    void addGraph(CostsGraph &graph);

    // prints out statistc of the graph set
    void printStat();

    GraphSet &operator=(const GraphSet &graph_set);

private:
    /* ####  member variables  #### */

    // list of graph pointers
    graph_list_type _graph_list;

    // set size
    int _size;
};

#endif
