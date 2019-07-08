#include <graphset.h>
#include <string.h>

GraphSet::GraphSet() {
    _graph_list = graph_list_type(0);
    _size = 0;
}

// init graph set with graph, which will be divided in its components
GraphSet::GraphSet(CostsGraph graph) {
    _graph_list = graph.getConnectedComponents();
    _size = _graph_list.size();
}


// inits a graph from file, function can be choosen to seperate graph in connected components
GraphSet::GraphSet(char *fname, double th, costs_parsing_fct_type fct, graph_set_parser_fct_type parser_fct)
        : _size(0) {
    _graph_list = graph_list_type(0);

    parser_fct(fname, *this, th, fct);
}

// inits a graph from file, function can be choosen to seperate graph in connected components
GraphSet::GraphSet(std::string fname, double th, costs_parsing_fct_type fct, graph_set_parser_fct_type parser_fct)
        : _size(0) {
    _graph_list = graph_list_type(0);

    char *str = strdup(fname.c_str());
    parser_fct(str, *this, th, fct);
}


// reads in a weight graph and using the appropriate costs parser to create costs graph(s) out of it
GraphSet::GraphSet(char *file_name, matrix_file_fct_type2 fct, costs_parsing_fct_type cost_fct, double th) : _size(0) {
    CostsGraph G = CostsGraph(file_name, fct, cost_fct, th);
    _graph_list = G.getConnectedComponents();
    _size = _graph_list.size();
}


// reads in a weight graph and using the appropriate costs parser to create costs graph(s) out of it
GraphSet::GraphSet(std::string file_name, matrix_file_fct_type2 fct, costs_parsing_fct_type cost_fct, double th)
        : _size(0) {
    CostsGraph G = CostsGraph(file_name, fct, cost_fct, th);
    _graph_list = G.getConnectedComponents();
    _size = _graph_list.size();
}

// inits graph with a set of graphs from a directory
GraphSet::GraphSet(char *dir_name, matrix_file_fct_type fct) : _size(0) {
    _graph_list = graph_list_type(0);

    DIR *dir;
    struct dirent *dirzeiger;

    // open directory
    if ((dir = opendir(dir_name)) != NULL) {
        // save directory pointer
        dirzeiger = readdir(dir);

        // iterate through directory
        while ((dirzeiger = readdir(dir))) {
            // get file name
            std::string fname = std::string(dir_name) + (*dirzeiger).d_name;
            int point_pos = fname.find_last_of('.');

            // check whether its a cm file
            if (fname.substr(point_pos + 1, point_pos + 3) == "cm") {
                CostsGraph *new_graph = new CostsGraph(fname, fct);
                _graph_list.insert(_graph_list.end(), new_graph);
                _size++;
            }
        }
    }

    closedir(dir);

}


// inits graph with a set of graphs from a directory
GraphSet::GraphSet(std::string dir_name, matrix_file_fct_type fct) : _size(0) {
    _graph_list = graph_list_type(0);

    char *str = strdup(dir_name.c_str());

    DIR *dir;
    struct dirent *dirzeiger;

    // open directory
    if ((dir = opendir(str)) != NULL) {

        // save directory pointer
        dirzeiger = readdir(dir);

        // iterate through directory
        while ((dirzeiger = readdir(dir))) {
            // get file name

            std::string fname = dir_name + (*dirzeiger).d_name;

            int point_pos = fname.find_last_of('.');

            // check whether its a cm file
            if (fname.substr(point_pos + 1, point_pos + 3) == "cm") {
                CostsGraph *new_graph = new CostsGraph(fname, fct);
                _graph_list.insert(_graph_list.end(), new_graph);
                _size++;
            }
        }
    }

    closedir(dir);
}


// copy constructor
GraphSet::GraphSet(GraphSet const &graph_set) {
    _size = graph_set._size;

    _graph_list = graph_list_type(_size);

    for (int i = 0; i < _size; i++) {
        _graph_list[i] = new CostsGraph(*(graph_set._graph_list[i]));
    }
}


GraphSet::~ GraphSet() {
    for (int i = 0; i < _size; i++) {
        delete (_graph_list[i]);
    }
}


// return set size
int GraphSet::getSetSize() {
    return _size;
}


// return reference to graph
CostsGraph &GraphSet::getGraph(int i) {
    if (i < 0 || i >= _size) {
        throw GraphException("invalid graph index");
    }

    return *_graph_list[i];
}


// delete graph from list
void GraphSet::deleteGraph(int i) {
    if (i < 0 || i >= _size) {
        throw GraphException("invalid graph index");
    }

    delete (_graph_list[i]);
    _graph_list.erase(_graph_list.begin() + i);
}


// add graph by making a deep copy of it and set pointer to it
void GraphSet::addGraph(CostsGraph &graph) {
    CostsGraph *new_graph = new CostsGraph(graph);
    _graph_list.insert(_graph_list.end(), new_graph);
    _size++;
}


// overloaded assignment operator = makes a deep copy
GraphSet &GraphSet::operator=(const GraphSet &graph_set) {
    if (this != &graph_set) {
        // delete old graphs
        for (int i = 0; i < _size; i++) {
            delete (_graph_list[i]);
        }
        _size = graph_set._size;

        // copy graphs to this object
        _graph_list = graph_list_type(_size);
        for (int i = 0; i < _size; i++) {
            _graph_list[i] = new CostsGraph(*(graph_set._graph_list[i]));
        }
    }
    return *this;
}


// print out statistic of graph set
void GraphSet::printStat() {
    int maximum = 0;
    int max = 1000;
    std::cout << "Number of connected components: " << _size << std::endl;
    int stat[max];
    for (int i = 0; i < max; i++) {
        stat[i] = 0;
    }

    for (int i = 0; i < _size; i++) {
        int size = (*_graph_list[i]).getSize() - 1;
        if (size < max) stat[size]++;
        if (size > maximum) maximum = size;
    }

    for (int i = 0; i < max; i++) {
        if (stat[i] != 0) {
            std::cout << "Number of connected componentes of size " << i + 1 << " : " << stat[i] << std::endl;
        }
    }
    std::cout << std::endl;
}
