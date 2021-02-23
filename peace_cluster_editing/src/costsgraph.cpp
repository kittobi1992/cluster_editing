/*+++++++++++++++++ class CostsGraph +++++++++++++++++++ */
#include <costsgraph.h>
#include <string.h>

/* #################  Constructors  ##################### */

// constructor with graph size and threshold
CostsGraph::CostsGraph(int size, double th) : _size(size), _threshold(th),
                                              _cost_matrix(size, forbidden) {
    // init with blank vertex names
    _info_list = vertices_names_type(_size);
    for (int i = 0; i < _size; i++) {
        vertex_info_type info = {vertex_name_type(1, ""), index_list_type(1, i)};
        _info_list[i] = info;
    }

    initMatrices();
}


// constructor with size, threshold and list of names
CostsGraph::CostsGraph(int size, std::string *edge_list, double th) : _size(size), _threshold(th),
                                                                      _cost_matrix(size, forbidden) {
    // init first vertex name with name from string array
    _info_list = vertices_names_type(_size);
    for (int i = 0; i < _size; i++) {
        vertex_info_type info = {vertex_name_type(1, ""), index_list_type(1, i)};
        _info_list[i] = info;
    }

    initMatrices();
}


// constructor with size, threshold and list of names
CostsGraph::CostsGraph(int size, vertex_name_type edge_list, double th) : _size(size), _threshold(th),
                                                                          _cost_matrix(size, forbidden) {
    // init first vertex name with name from string vector
    _info_list = vertices_names_type(_size);
    for (int i = 0; i < _size; i++) {
        vertex_info_type info = {vertex_name_type(1, edge_list[i]), index_list_type(1, i)};
        _info_list[i] = info;
    }

    initMatrices();
}


// constructor with size, threshold and costs matrix
CostsGraph::CostsGraph(int size, double_matrix_type costs, double th) : _size(size), _threshold(th),
                                                                        _cost_matrix(size, forbidden) {

    // throw error if weight matrix has wrong size
    if (costs.size() != _size) {
        throw GraphException("costs matrix has no valid format");
    }

    // check matrix for triangular shape
    for (int i = 0; i < costs.size(); i++) {
        if (costs[i].size() != i) {
            throw GraphException("costs matrix has no valid format");
        }
    }

    // convert and copy costs matrix to linear cost matrix
    for (int i = 0; i < costs.size(); i++) {
        for (int k = 0; k < i; k++) {
            _cost_matrix.pos(i, k) = costs[i][k];
        }
    }

    // initiate vertex names with emtpy names
    _info_list = vertices_names_type(_size);
    for (int i = 0; i < _size; i++) {
        vertex_info_type info = {vertex_name_type(1, ""), index_list_type(1, i)};
        _info_list[i] = info;
    }
}


// constructor with size, threshold, list of names and cost matrix
CostsGraph::CostsGraph(int size, double_matrix_type costs, vertex_name_type edge_list, double th) : _size(size),
                                                                                                    _threshold(th),
                                                                                                    _cost_matrix(size,
                                                                                                                 forbidden) {
    // throw error if weight matrix has wrong size
    if (costs.size() != _size) {
        throw GraphException("costs matrix has no valid format");
    }

    // check matrix for triangular shape
    for (int i = 0; i < costs.size(); i++) {
        if (costs[i].size() != i) {
            throw GraphException("costs matrix has no valid format");
        }
    }

    // convert and copy costs matrix to linear cost matrix
    for (int i = 0; i < costs.size(); i++) {
        for (int k = 0; k < i; k++) {
            _cost_matrix.pos(i, k) = costs[i][k];
        }
    }

    // initiate names list with given names
    _info_list = vertices_names_type(_size);
    for (int i = 0; i < _size; i++) {
        vertex_info_type info = {vertex_name_type(1, edge_list[i]), index_list_type(1, i)};
        _info_list[i] = info;
    }

}


// constructor with size, threshold, weight_matrix and a cost parsing function, with
// which the weights will be converted to costs
CostsGraph::CostsGraph(int size, double_matrix_type weight_matrix, costs_parsing_fct_type cost_fct, double th) : _size(
        size), _threshold(th), _cost_matrix(size, forbidden) {
    // throw error if weight matrix has wrong size
    if (weight_matrix.size() != _size) {
        throw GraphException("weight matrix has no valid format");
    }

    // check matrix for triangular shape + create cost matrix
    for (int i = 0; i < size; i++) {
        if (weight_matrix[i].size() != i) {
            throw GraphException("weight matrix has no valid format");
        }

        for (int k = 0; k < i; k++) {
            _cost_matrix.pos(i, k) = cost_fct(weight_matrix[i][k], _threshold);
        }
    }

    // initiate names with empty names
    _info_list = vertices_names_type(_size);
    for (int i = 0; i < _size; i++) {
        vertex_info_type info = {vertex_name_type(1, ""), index_list_type(1, i)};
        _info_list[i] = info;
    }

}


// constructor with size, threshold, weight_matrix, vertex names and a cost parsing function, with
// which the weights will be converted to costs
CostsGraph::CostsGraph(int size, double_matrix_type weight_matrix, costs_parsing_fct_type cost_fct,
                       vertex_name_type edge_list, double th) : _size(size), _threshold(th),
                                                                _cost_matrix(size, forbidden) {
    // throw error if weight matrix has wrong size
    if (weight_matrix.size() != _size) {
        throw GraphException("weight matrix has no valid format");
    }

    // check matrix for triangular shape + create cost matrix
    for (int i = 0; i < _size; i++) {
        if (weight_matrix[i].size() != i) {
            throw GraphException("weight matrix has no valid format");
        }

        for (int k = 0; k < i; k++) {
            _cost_matrix.pos(i, k) = cost_fct(weight_matrix[i][k], _threshold);
        }
    }

    // initiate names with given names
    _info_list = vertices_names_type(_size);
    for (int i = 0; i < _size; i++) {
        vertex_info_type info = {vertex_name_type(1, edge_list[i]), index_list_type(1, i)};
        _info_list[i] = info;
    }

}


// constructor with file name, double_file parsing function, cost parsing function and threshold
// the file will be parsed with the double file fct and every read weight value will be converted
// to costs with the costs parsing function using the threshold
/* CostsGraph::CostsGraph(char *file_name, double_file_fct_type fct, costs_parsing_fct_type cost_fct, double th) : _size(
        0), _threshold(th), _cost_matrix() {
    try {
        initMatrices(file_name, fct, cost_fct);
    } catch (GraphException &e) {
        throw GraphException(e);
    }
} */


// constructor with file name, double_file parsing function, cost parsing function and threshold
// the file will be parsed with the double file fct and every read weight value will be converted
// to costs with the costs parsing function using the threshold
CostsGraph::CostsGraph(const std::string &file_name, double_file_fct_type fct, costs_parsing_fct_type cost_fct,
                       double th)
        : _size(0), _threshold(th), _cost_matrix() {
    try {
        char *str = strdup(file_name.c_str());
        initMatrices(str, fct, cost_fct);
    } catch (GraphException &e) {
        throw GraphException(e);
    }
}

// constructor like above but reads in a matrix file which needs still cost parsing
/* CostsGraph::CostsGraph(char *file_name, matrix_file_fct_type2 fct, costs_parsing_fct_type cost_fct, double th) : _size(
        0),
                                                                                                                 _threshold(
                                                                                                                         th),
                                                                                                                 _cost_matrix() {
    if (th < 0) {
        throw GraphException("no negative threshold allowed");
    }

    try {
        initMatrices(file_name, fct, cost_fct);
    } catch (GraphException &e) {
        throw GraphException(e);
    }
} */

// constructor like above but reads in a matrix file which needs still cost parsing
CostsGraph::CostsGraph(const std::string &file_name, matrix_file_fct_type2 fct, costs_parsing_fct_type cost_fct,
                       double th)
        : _size(0),
          _threshold(th), _cost_matrix() {
    try {
        char *str = strdup(file_name.c_str());
        initMatrices(str, fct, cost_fct);
    } catch (GraphException &e) {
        throw GraphException(e);
    }
}


// constructor with file name and a matrix file function, which just reads the raw matrix as cost graph
// the threshold is no longer important...all other infos like names will also be read from file
/* CostsGraph::CostsGraph(char *file_name, matrix_file_fct_type fct) : _size(0), _threshold(1e-20),
                                                                    _cost_matrix() {
    fct(file_name, *this);
} */


// constructor with file name and a matrix file function, which just reads the raw matrix as cost graph
// the threshold is no longer important...all other infos like names will also be read from file
CostsGraph::CostsGraph(const std::string &file_name, matrix_file_fct_type fct) : _size(0), _threshold(1e-20), _cost_matrix() {
    char *str = strdup(file_name.c_str());
    fct(str, *this);
}

// cop constructor
CostsGraph::CostsGraph(CostsGraph const &new_graph) : _size(new_graph._size),
                                                      _cost_matrix(new_graph._cost_matrix) {
    // get additional information
    _info_list = new_graph._info_list;
    _threshold = new_graph._threshold;
}


/* #################  help/private functions  ##################### */


// initializes the adjacency matrix as triangular matrix
void CostsGraph::initMatrices() {}


// init cost matrix with a double file fct and cost parsing function, which
// fills converts the read weights to costs
void CostsGraph::initMatrices(char *fname, double_file_fct_type fct, costs_parsing_fct_type cost_fct) {
    fct(fname, *this, cost_fct);

}

// fills converts the read weights to costs
void CostsGraph::initMatrices(char *fname, matrix_file_fct_type2 fct, costs_parsing_fct_type cost_fct) {
    fct(fname, *this, cost_fct, _threshold);
}


// converts element of any type to string
template<typename T>
inline std::string CostsGraph::valueToString(T x) const {
    std::ostringstream o;
    o << x;
    return o.str();
}


// converts byte-matrix to string by adding the edge names in the first column and row
// simple option for convert positive values to 1, negative to 0 and forbidden to F and permanent to P
std::string CostsGraph::matrixToString(bool simple) const {
    std::string output = "\t";

    // create header = first line of vertex names
    for (int i = 0; i < _size; i++) {
        output += _info_list[getIntern(i)].names[0] + "\t";
    }

    output += "\n";

    // for each line ..
    for (int i = 0; i < _size; i++) {
        // .. add first vertex name ...
        output += _info_list[getIntern(i)].names[0] + "\t";

        // .. then values
        for (int k = 0; k < i; k++) {
            double element = getEdge(i, k);
            if (!simple) {
                output += valueToString<double>(element) + "\t";
            } else {
                if (element > 0) {
                    if (element == permanent) {
                        output += "P\t";
                    } else {
                        output += "1\t";
                    }
                } else {
                    if (element == forbidden) {
                        output += "F\t";
                    } else {
                        output += "0\t";
                    }

                }
            }
        }
        output += "\n";
    }

    return output;
}


/* ################## public functions ################# */

/* ------------------ access functions ----------------- */

// set edge to given adjacency status
void CostsGraph::setEdge(int i, int j, double value) {
    // checks that edge i=j never exists
    if (i == j) {
        throw GraphException(" no self incident edges are allowed ");
    }

    // depended on indices add value to triangular adjacency matrix
    _cost_matrix.pos(i, j) = value;
}


// return internal index of i
inline int CostsGraph::getIntern(int i) const {
    return _cost_matrix.getInternalIndex(i);
}


// get next neighbor of node, starting at k
// returns -1 for no neighbor found
short CostsGraph::getNeighbor(int i, int k) const {
    k = (k == i) ? k + 1 : k;
    while (k < _size && getEdge(i, k) <= 0) {
        k++;
        k = (k == i) ? k + 1 : k;
    }
    if (k >= _size) {
        return -1;
    } else {
        return k;
    }
}


// sets edge weight by using the cost-parsing function to convert a double value to
// an edge weight
void CostsGraph::setEdgeWeight(int i, int j, double value, costs_parsing_fct_type fct) {
    // calls cost-parsing function with given values
    double costs = fct(value, _threshold);
    this->setEdge(i, j, costs);
}


// returns complete name of vertex divided by |
std::string CostsGraph::getVertexName(int i) const {
    if (i < 0 || i > _size || _info_list[i].names.empty()) {
        throw GraphException("vertex not defined");
    }

    // translate to internal index
    i = getIntern(i);

    // add up name to big string
    std::string name = _info_list[i].names[0];
    for (int k = 1; k < _info_list[i].names.size(); k++) {
        name += "|" + _info_list[i].names[k];
    }

    return name;
}


// set first name to value
void CostsGraph::setVertexName(int i, std::string name) {
    if (i < 0 || i > _size) {
        throw GraphException("vertex not defined");
    }

    // translate to internal index
    i = getIntern(i);

    _info_list[i].names[0] = std::move(name);
}


// return complete vertex info of node
CostsGraph::vertex_info_type CostsGraph::getVertexInfo(int i) const {
    if (i < 0 || i > _size) {
        throw GraphException("vertex not defined");
    }

    // translate to internal index
    i = getIntern(i);

    return _info_list[i];
}


// set complete vertex info of a node
void CostsGraph::setVertexInfo(int i, CostsGraph::vertex_info_type info) {
    if (i < 0 || i > _size) {
        throw GraphException("vertex not defined");
    }

    // translate to internal index
    i = getIntern(i);

    _info_list[i] = std::move(info);
}


// convert the names of all nodes to a string, each node gets a numbering
std::string CostsGraph::vertexNamesToString() const {
    std::string output;
    for (int i = 0; i < _size - 1; i++) {
        output += valueToString<int>(i) + "=" + getVertexName(i) + ", ";
    }
    output += valueToString<int>(_size - 1) + "=" + getVertexName(_size - 1);
    return output;
}


// returns index list of a node
CostsGraph::index_list_type CostsGraph::getVertexIndices(int i) const {
    if (i < 0 || i > _size || _info_list[i].indices.empty()) {
        throw GraphException("vertex not defined");
    }

    // translate to internal index
    i = getIntern(i);

    return _info_list[i].indices;
}


// returns name list of a node
CostsGraph::vertex_name_type CostsGraph::getVertexNames(int i) const {
    if (i < 0 || i > _size || _info_list[i].indices.empty()) {
        throw GraphException("vertex not defined");
    }

    // translate to internal index
    i = getIntern(i);

    return _info_list[i].names;
}


// returns graph size
int CostsGraph::getSize() const {
    return _size;
}


// return index of vertex with this name
int CostsGraph::getIndex(const std::string &key) const {
    for (int i = 0; i < _size; i++) {
        for (const auto& name : _info_list[getIntern(i)].names) {
            if (name == key) {
                return i;
            }
        }
    }

    return -1;
}

// checks if vertex name is already defined
bool CostsGraph::inGraph(const std::string &key) const {
    for (int i = 0; i < _size; i++) {
        for (const auto& name : _info_list[getIntern(i)].names) {
            if (name == key) {
                return true;
            }
        }
    }

    return false;
}


// returns true if any of the names in vector is in the graph defined
bool CostsGraph::inGraph(const vertex_name_type &names) const {
    for (const auto& name : names) {
        if (this->inGraph(name)) {
            return true;
        }
    }

    return false;
}


// calculates degree of vertex
int CostsGraph::getDegree(int i) const {
    int sum = 0;

    // count positive adjacencies
    for (int k = 0; k < _size; k++) {
        if (k != i) {
            sum += (getEdge(i, k) > 0) ? 1 : 0;
        }
    }
    return sum;
}


// return number of set edges
int CostsGraph::getEdgeNumber() const {
    int sum = 0;

    for (int i = 0; i < _size; i++) {
        for (int k = 0; k < i; k++) {
            sum += (getEdge(i, k) > 0) ? 1 : 0;
        }
    }

    return sum;
}

// returns threshold
double CostsGraph::getThreshold() const {
    return _threshold;
}


// returns true of i is member of clique, false if not
bool CostsGraph::inClique(int i) const {
    return !(this->getClique(i).empty());
}


// returns vector of clique members, empty if i is not in a clique
CostsGraph::byte_vector_type CostsGraph::getClique(int i) const {
    // create vector of clique members, initiated no members
    byte_vector_type clique;

    // fill clique members with neighbors of i and i itself...sorted
    for (int k = 0; k < i; k++) {
        if (getEdge(i, k) > 0) {
            clique.push_back(k);
        }
    }
    clique.insert(clique.end(), i);
    for (int k = i + 1; k < _size; k++) {
        if (getEdge(i, k) > 0) {
            clique.push_back(k);
        }
    }

    // check members of clique for right neighbors
    for (int k = 0; k < clique.size(); k++) {
        // count neighbors of clique member
        int neighbors = 0;
        for (int h = 0; h < clique[k]; h++) {
            neighbors += (getEdge(clique[k], h) > 0) ? 1 : 0;
        }
        for (int h = clique[k] + 1; h < _size; h++) {
            neighbors += (getEdge(clique[k], h) > 0) ? 1 : 0;
        }
        if (neighbors + 1 != clique.size()) {
            return byte_vector_type(0);
        }

        // check all clique members for edges
        for (int h = 0; h < k; h++) {
            if (getEdge(clique[h], clique[k]) <= 0) {
                return byte_vector_type(0);
            }
        }
        for (int h = k + 1; h < clique.size(); h++) {
            if (getEdge(clique[h], clique[k]) <= 0) {
                return byte_vector_type(0);
            }
        }
    }

    return clique;
}


// deletes vertex from graph, just virtually...size will be decreased
// but not internal matrix
void CostsGraph::deleteVertex(int index) {
    _cost_matrix.deleteIndex(index);

    // decrease size variable of object
    _size--;

}


void CostsGraph::deleteVertices(CostsGraph::byte_vector_type set) {
    sort(set.begin(), set.end());
    for (int i = set.size() - 1; i >= 0; i--) {
        deleteVertex(i);
    }
}

// merges to vertices and gives them the given weight vector as weights to the other vertices
void CostsGraph::mergeVertices(int index1, int index2, const double_array_type &costs) {
    // save internal indices
    int intern_index1 = getIntern(index1);
    int intern_index2 = getIntern(index2);

    // be sure that intern_index2 is the node which disappears
    if (intern_index1 > intern_index2) {
        std::swap(intern_index1, intern_index2);
    }

    // add names of vertex 2 to vertex 1
    _info_list[intern_index1].names.insert(_info_list[intern_index1].names.end(),
                                           _info_list[intern_index2].names.begin(),
                                           _info_list[intern_index2].names.end());
    _info_list[intern_index1].indices.insert(_info_list[intern_index1].indices.end(),
                                             _info_list[intern_index2].indices.begin(),
                                             _info_list[intern_index2].indices.end());

    // be sure that index2 node is the deleted one
    if (index1 > index2) {
        std::swap(index1, index2);
    }

    // delete node index2
    deleteVertex(index2);

    // save new costs to node at position index1
    for (int i = 0; i < index1; i++) {
        setEdge(i, index1, costs[i]);
    }
    for (int i = index1 + 1; i < _size; i++) {
        setEdge(i, index1, costs[i - 1]);
    }
}


// deletes clique if vertex 'index' is element of clique, returns true if
// deletion was successful
bool CostsGraph::deleteClique(int index) {
    // get clique
    byte_vector_type clique = this->getClique(index);
    if (!clique.empty()) {
        // sort clique
        std::sort(clique.begin(), clique.end());

        // delete vertices by starting with biggest index
        for (int k = clique.size() - 1; k >= 0; k++) {
            this->deleteVertex(clique[k]);
        }

        return true;

    } else {
        return false;
    }
}


// calculate connected components and save them in a matrix,
// each row corresponds to the indices of one component
CostsGraph::vertex_matrix_type CostsGraph::getConnectedVertices() const {
    VertexLists vertex_lists = VertexLists(_size, _size);

    // init every list with one vertex
    for (int i = 0; i < _size; i++) {
        vertex_lists.addToList(i, i);
    }

    // check every vertex for neighbors and merge connected components so far
    for (int i = 0; i < _size; i++) {
        for (int k = 0; k < i; k++) {
            if (getEdge(i, k) > 0) {
                unsigned short list_i = vertex_lists.getListNr(i);
                unsigned short list_k = vertex_lists.getListNr(k);
                if (list_i != list_k) {
                    vertex_lists.mergeLists(list_i, list_k);
                }
            }
        }
    }

    vertex_matrix_type output(0);

    // check for non-empty lists and create new graph
    for (int i = 0; i < _size; i++) {
        int sub_graph_size = vertex_lists.getListSize(i);
        if (sub_graph_size != 0) {
            output.push_back(vertex_lists.getList(i));
        }
    }

    return output;
}


// divide the graph into components, each again a costs graph
CostsGraph::graph_list_type CostsGraph::getConnectedComponents() const {
    // get cliques
    vertex_matrix_type vertex_matrix = getConnectedVertices();
    graph_list_type graph_list = graph_list_type(0);

    // check for non-empty lists and create new graph
    for (const auto& component : vertex_matrix) {
        int sub_graph_size = component.size();

        // create new graph with the appropriate size
        auto new_graph = new CostsGraph(sub_graph_size, _threshold);

        // fill graph with names/info and values
        for (int k = 0; k < sub_graph_size; k++) {
            (*new_graph).setVertexInfo(k, this->getVertexInfo(component[k]));
            for (int m = 0; m < k; m++) {
                (*new_graph).setEdge(k, m, this->getEdge(component[k], component[m]));
            }
        }
        // save graph in return list
        graph_list.insert(graph_list.begin(), new_graph);
    }

    return graph_list;
}


// dependend on the given connected components, saved in the vertex matrix
// the graph will be divided into costs graph again
CostsGraph::graph_list_type CostsGraph::getConnectedComponents(const vertex_matrix_type &vertex_matrix) const {
    graph_list_type graph_list;

    // check for non-empty lists and create new graph
    for (const auto& component : vertex_matrix) {
        int sub_graph_size = component.size();

        // create new graph with the appropriate size
        auto new_graph = new CostsGraph(sub_graph_size, _threshold);

        // fill graph with names/info and values
        for (int k = 0; k < sub_graph_size; k++) {
            (*new_graph).setVertexInfo(k, this->getVertexInfo(component[k]));
            for (int m = 0; m < k; m++) {
                (*new_graph).setEdge(k, m, this->getEdge(component[k], component[m]));
            }
        }
        // save graph in return list
        graph_list.push_back(new_graph);
    }

    return graph_list;
}

/* #################### IO-functions ################ */

// creates string from costs matrix
std::string CostsGraph::toString(bool simple) const {
    return matrixToString(simple);
}


// do IO operation with given function
void CostsGraph::costsIOOperation(char *fname, matrix_file_fct_type fct) {
    fct(fname, *this);
}


// do IO operation with given function
void CostsGraph::costsIOOperation(const std::string &fname, matrix_file_fct_type fct) {
    char *str = strdup(fname.c_str());
    fct(str, *this);
}


/* #################### overloaded operators ################## */


// costsGraph to string - cast 
CostsGraph::operator std::string() {
    return toString(false);
}


// costsGraph to stream
std::ostream &operator<<(std::ostream &output, const CostsGraph &g) {
    return output << g.toString(false);
}


// assignment operator
CostsGraph &CostsGraph::operator=(const CostsGraph &right) {
    // copy information
    _info_list = right._info_list;
    _size = right._size;
    _threshold = right._threshold;
    _cost_matrix = right._cost_matrix;

    return *this;
}
