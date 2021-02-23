/*+++++++++++++++++ class MatrixParser +++++++++++++++++++ */
#ifndef __matrixparser_h__
#define __matrixparser_h__

// include for input-output
#include <fstream>
#include <iostream>
#include <sstream>

// includes from stlib
#include <string>
#include <vector>
#include <map>

#include <costsgraph.h>

class MatrixParser {
public:
    /* ####  constructors  #### */
    MatrixParser() = default;

    /* ####  destructor  #### */
    ~MatrixParser() = default;

    // function pointer to cost parsing function
    typedef double (*costs_parsing_fct_type)(double value, double threshold);

    /* ####  access functions  #### */
    static void toMatrixFile(char *fname, CostsGraph &G);

/*
	static void fileToGraph(char* fname, CostsGraph& G, double *matrix);
	static void fileToGraph(char* fname, CostsGraph& G, int *matrix);
	static void fileToGraph(char* fname, CostsGraph& G, short *matrix);
	static void fileToGraph(char* fname, CostsGraph& G, char *matrix);
*/
    static void initFromMatrixFile(char *fname, CostsGraph &G);

    static void initFromWeightMatrixFile(char *fname, CostsGraph &G, costs_parsing_fct_type fct, double threshold);

private:
    /* ####  member variables  #### */
    template<typename Type>
    static inline std::string toString(Type x);

    template<typename Type>
    static inline Type stringToType(std::string &s);

    template<typename Type>
    static inline Type stringToType(char *s);

/*
	template <typename Type>
	static void fromMatrixFile(char* fname, CostsGraph& G, Type *matrix);
*/
};

#endif
