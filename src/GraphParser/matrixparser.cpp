#include <matrixparser.h>
#include <string.h>

void MatrixParser::toMatrixFile(char *fname, CostsGraph &G) {
    int size = G.getSize();

    // create file
    std::ofstream graph_file(fname);
    if (!graph_file) {
        throw GraphException("could not create file ");
    }

    // write header
    graph_file << size << std::endl;
    //for (int i = size - 1; i >= 0; i--) {
    for (int i = 0; i < size; i++) {
        graph_file << G.getVertexName(i) << std::endl;
    }

    // write each line with name at the beginning
    //for (int i = size - 1; i >= 0; i--) {
    for (int i = 0; i < size; i++) {
        //for (int k = 0; k < i; k++) {
        for (int k = i + 1; k < size; k++) {
            graph_file << toString<double>(G.getEdge(k, i));
            graph_file << "\t";
        }
        graph_file << std::endl;
    }

    // close output file
    graph_file.close();
}


template<typename Type>
inline Type MatrixParser::stringToType(std::string &s) {
    std::istringstream i(s);
    Type x;
    i >> x;
    return x;
}

template<typename Type>
inline Type MatrixParser::stringToType(char *s) {
    std::istringstream i(s);
    Type x;
    i >> x;
    return x;
}


template<typename Type>
inline std::string MatrixParser::toString(Type x) {
    std::ostringstream o;
    o << x;
    return o.str();
}


/*
template <typename Type>
void MatrixParser::fromMatrixFile(char* fname, CostsGraph& G, Type *matrix)
{
	int size = G.getSize();

	// open file
	std::ifstream matrix_file(fname);	
	if (!matrix_file) {
		throw GraphException( "could not open file " );
	}
	
	std::string row;
	// get header line
	getline(matrix_file, row);
	char* str = strdup(row.c_str());
	char* element = strtok(str, "\t");
	
	for (int i = 0; i < size; i++ ){
		element = strtok(NULL, "\t");
		if (element == NULL || G.getVertexName(i) != element) {
		     throw GraphException( "matrix file is not matching to graph" );
		}
	}
	
	// read in file line-wise
	for (int i = 0; i < size; i++) {
		// get line from file
		getline(matrix_file,row);
		if (row.empty()) {
		     throw GraphException( "matrix file is not matching to graph" );
		}
		str = strdup(row.c_str());

		// exctract vertex name
		element = strtok(str, "\t");
		if (G.getVertexName(i) != element) {
		     throw GraphException( "matrix file is not matching to graph" );
		}
	
		// get values
		for(int k = 0; k < i; k++) {
		     element = strtok(NULL, "\t");
		     if (element == NULL) {
		          throw GraphException( "matrix file is not matching to graph" );
		     }
		     matrix[matrix_pos(i,k)] = stringToType<Type>(element);
		}
	}
	matrix_file.close();
}
*/
/*
void MatrixParser::toMatrixFile(char* fname, CostsGraph& G, double *matrix)
{
	toMatrixFile<double>(fname, G, matrix);
}

void MatrixParser::toMatrixFile(char* fname, CostsGraph& G, int *matrix)
{
	toMatrixFile<int>(fname, G, matrix);
}

void MatrixParser::toMatrixFile(char* fname, CostsGraph& G, short *matrix)
{
	toMatrixFile<short>(fname, G, matrix);
}

void MatrixParser::toMatrixFile(char* fname, CostsGraph& G, char *matrix)
{
	toMatrixFile<char>(fname, G, matrix);
}
*/
/*
void MatrixParser::fileToGraph(char* fname, CostsGraph& G, double *matrix)
{
	fromMatrixFile<double>(fname, G, matrix);
}

void MatrixParser::fileToGraph(char* fname, CostsGraph& G, int *matrix)
{
	fromMatrixFile<int>(fname, G, matrix);
}

void MatrixParser::fileToGraph(char* fname, CostsGraph& G, short *matrix)
{
	fromMatrixFile<short>(fname, G, matrix);
}

void MatrixParser::fileToGraph(char* fname, CostsGraph& G, char *matrix)
{
	fromMatrixFile<char>(fname, G, matrix);
}
*/
void MatrixParser::initFromMatrixFile(char *fname, CostsGraph &G) {
    // open file
    std::ifstream matrix_file(fname);
    if (!matrix_file) {
        throw GraphException("could not open file ");
    }

    // variable to save line
    std::string row;

    // get first line = get size
    getline(matrix_file, row);
    int size = stringToType<int>(row);

    // create object for names
    CostsGraph::vertex_name_type names = CostsGraph::vertex_name_type(size);
    // get names line by line
    for (int i = 0; i < size; i++) {
        getline(matrix_file, row);
        names[size - i - 1] = row;
    }

    // create graph
    CostsGraph G2 = CostsGraph(size, names);
    G = G2;

    // read in matrix line-wise
    for (int i = size - 1; i > 0; i--) {
        // get line from file
        getline(matrix_file, row);

        if (row.empty()) {
            throw GraphException("matrix file is not matching to graph");
        }
        char *str = strdup(row.c_str());

        // get values
        char *element = strtok(str, "\t");
        for (int k = i - 1; k >= 0; k--) {
            if (element == NULL) {
                throw GraphException("matrix file is not matching to graph");
            }
            G.setEdge(i, k, stringToType<double>(element));
            element = strtok(NULL, "\t");
        }
    }
    matrix_file.close();
}

void MatrixParser::initFromWeightMatrixFile(char *fname, CostsGraph &G, costs_parsing_fct_type fct, double threshold) {
    // open file
    std::ifstream matrix_file(fname);
    if (!matrix_file) {
        throw GraphException("could not open file ");
    }

    // variable to save line
    std::string row;

    // create help vector for first line
    std::vector<double> help = std::vector<double>(0);

    // get first line from file
    getline(matrix_file, row);
    if (row.empty()) {
        throw GraphException("empty matrix file");
    }
    char *str = strdup(row.c_str());

    // get values
    char *element = strtok(str, "\t\n");
    while (element != NULL) {
        help.insert(help.end(), stringToType<double>(element));
        element = strtok(NULL, "\t\n");
    }

    // save size
    int size = help.size();

    // create object for names
    CostsGraph::vertex_name_type names = CostsGraph::vertex_name_type(size);
    // create names
    for (int i = 0; i < size; i++) {
        names[i] = toString<int>(i);
    }

    // create matrix
    std::vector<std::vector<double> > weights = std::vector<std::vector<double> >(size);
    weights[0] = help;


    // read in matrix line-wise
    for (int i = 1; i < size; i++) {
        weights[i] = std::vector<double>(size);
        // get line from file
        getline(matrix_file, row);

        if (row.empty()) {
            throw GraphException("matrix file is not matching to graph");
        }
        char *str = strdup(row.c_str());

        // get values
        char *element = strtok(str, "\t");
        for (int k = 0; k < size; k++) {
            if (element == NULL) {
                throw GraphException("matrix file is not matching to graph");
            }
            weights[i][k] = stringToType<double>(element);
            element = strtok(NULL, "\t");
        }
    }
    matrix_file.close();

    // create graph
    CostsGraph G2 = CostsGraph(size, names);
    G = G2;

    // fill graph using the cost parser
    for (int i = 0; i < size; i++) {
        for (int k = 0; k < i; k++) {
            G.setEdge(i, k, fct((weights[i][k] + weights[k][i]) / 2, threshold));
        }
    }

    //std::cout << "created graph: " << G << std::endl;
}
