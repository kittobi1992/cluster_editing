/*+++++++++++++++++ class EdgeFileParser +++++++++++++++++++ */
#ifndef __edgefileparser_h__
#define __edgefileparser_h__

// include for input-output
#include <fstream>
#include <iostream>
#include <sstream>

// includes from stlib
#include <string>
#include <vector>
#include <map>

#include <vertexlists.h>

#include <../costsgraph.h>
#include <Exceptions/graphexception.h>
#include <../graphset.h>

class EdgeFileParser {
public:
    typedef std::vector<std::vector<double> > double_matrix_type;
    typedef std::vector<std::vector<int> > int_matrix_type;
    typedef std::vector<std::vector<short> > short_matrix_type;
    typedef std::vector<std::vector<char> > char_matrix_type;

    struct problem {
        double forward;
        double backward;
        int index1;
        int index2;
    };

    struct cc_element {
        int component;
        int index;
    };

    // pos integer matrix and array for saving connected components
    typedef VertexLists::vertex_set_type vertex_set_type;
    typedef std::vector<vertex_set_type> vertex_matrix_type;

    // function pointer to cost parsing function
    typedef double (*costs_parsing_fct_type)(double blast_value, double threshold);

    /* ####  constructors  #### */
    EdgeFileParser() = default;

    /* ####  destructor  #### */
    ~EdgeFileParser() = default;

    /* ####  access functions  #### */
    static void to5ColumnFile(char *fname, CostsGraph &G, double *matrix);

    static void to5ColumnFile(char *fname, CostsGraph &G, int *matrix);

    static void to5ColumnFile(char *fname, CostsGraph &G, short *matrix);

    static void to5ColumnFile(char *fname, CostsGraph &G, char *matrix);

    static void initFrom5ColumnFile(char *fname, CostsGraph &G, costs_parsing_fct_type fct);

    static void initFrom3ColumnFile(char *fname, CostsGraph &G, costs_parsing_fct_type fct);

    static void initFrom12ColumnFile(char *fname, CostsGraph &G, costs_parsing_fct_type fct);

    static void initGraphSetFrom5ColumnFile(char *fname, GraphSet &graph_set, double th, costs_parsing_fct_type fct);

    static void initGraphSetFrom3ColumnFile(char *fname, GraphSet &graph_set, double th, costs_parsing_fct_type fct);

    static void initGraphSetFrom12ColumnFile(char *fname, GraphSet &graph_set, double th, costs_parsing_fct_type fct);

private:
    /* ####  member variables  #### */
    template<typename Type>
    static inline std::string toString(Type x);

    template<typename Type>
    static void to5ColumnFile(char *fname, CostsGraph &G, Type *matrix);

    template<typename Type>
    static inline Type stringToType(std::string &s);

    template<typename Type>
    static inline Type stringToType(char *s);

    template<typename T>
    static inline void swap(T &i, T &j);

    //inline static double stringToDouble(std::string& s) const;
    static void get5ColumnContent(std::string line, std::string &gene_name1, std::string &gene_name2, double &score);

    static void get3ColumnContent(std::string line, std::string &gene_name1, std::string &gene_name2, double &score);

    static void get12ColumnContent(std::string line, std::string &gene_name1, std::string &gene_name2, double &score);

    static void initFromXColumnFile(char *fname, CostsGraph &G, costs_parsing_fct_type fct, int x);

    //init a graph set
    static vertex_matrix_type getConnectedVertices(char_matrix_type matrix);

    template<typename Type>
    static void extendMatrix(std::vector<std::vector<Type> > &matrix, Type init_value);

    static void
    initGraphSetFromXColumnFile(char *fname, GraphSet &graph_set, double th, costs_parsing_fct_type fct, int x);

};

#endif
