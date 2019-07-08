#include <edgefileparser.h>
#include <string.h>


template<typename Type>
void EdgeFileParser::to5ColumnFile(char *fname, CostsGraph &G, Type *matrix) {
    int size = G.getSize();

    // create file
    std::ofstream graph_file(fname);
    if (!graph_file) {
        throw GraphException("could not create file ");
    }

    // iterate through edges
    for (int i = 0; i < size; i++) {
        for (int k = 0; k < i; k++) {
            if (G.getEdge(i, k) != G.forbidden) {
                graph_file << G.getVertexName(i) << "\t " << "\t" << G.getVertexName(k) << "\t " << "\t";
                graph_file << toString<double>(matrix[i / 2 * (i + 1) + k]);
                graph_file << "\n";
            }
        }
    }

    // close output file
    graph_file.close();
}


template<typename Type>
inline Type EdgeFileParser::stringToType(std::string &s) {
    std::istringstream i(s);
    Type x;
    i >> x;
    return x;
}

template<typename Type>
inline Type EdgeFileParser::stringToType(char *s) {
    std::istringstream i(s);
    Type x;
    i >> x;
    return x;
}


template<typename Type>
inline std::string EdgeFileParser::toString(Type x) {
    std::ostringstream o;
    o << x;
    return o.str();
}


void
EdgeFileParser::get5ColumnContent(std::string line, std::string &gene_name_1, std::string &gene_name_2, double &score) {
    char *str = strdup(line.c_str());
    gene_name_1 = strtok(str, "\t");
    std::string gene_pos_1 = strtok(NULL, "\t");
    gene_name_2 = strtok(NULL, "\t");
    std::string gene_pos_2 = strtok(NULL, "\t");

    std::string blast_score;
    blast_score = strtok(NULL, "\n\t\r");

    score = stringToType<double>(blast_score);
}


void
EdgeFileParser::get3ColumnContent(std::string line, std::string &gene_name_1, std::string &gene_name_2, double &score) {
    std::cout << "line in fkt: " << line << std::endl;

    char *str = strdup(line.c_str());

    gene_name_1 = strtok(str, "\t");

    gene_name_2 = strtok(NULL, "\t");

    std::string blast_score;
    blast_score = strtok(NULL, "\n\t\r");

    score = stringToType<double>(blast_score);
}


void EdgeFileParser::get12ColumnContent(std::string line, std::string &gene_name_1, std::string &gene_name_2,
                                        double &score) {
    char *str = strdup(line.c_str());
    gene_name_1 = strtok(str, "\t");
    gene_name_2 = strtok(NULL, "\t");
    // get 8 other columns
    for (int i = 0; i < 8; i++) {
        std::string blast_score = strtok(NULL, "\t");
    }

    std::string blast_score;
    blast_score = strtok(NULL, "\n\t\r");

    score = stringToType<double>(blast_score);
}

template<typename T>
void EdgeFileParser::swap(T &i, T &j) {
    T help = i;
    i = j;
    j = help;
}

EdgeFileParser::vertex_matrix_type EdgeFileParser::getConnectedVertices(char_matrix_type matrix) {
    int size = matrix.size();

    VertexLists vertex_lists = VertexLists(size, size);

    // init every list with one vertex
    for (int i = 0; i < size; i++) {
        vertex_lists.addToList(i, i);
    }

    // check every vertex for neighbors and merge connected components so far
    for (int i = 0; i < size; i++) {
        for (int k = 0; k < i; k++) {
            if (matrix[i][k] == '+') {
                int list_i = vertex_lists.getListNr(i);
                int list_k = vertex_lists.getListNr(k);
                if (list_i != list_k) {
                    vertex_lists.mergeLists(list_i, list_k);
                }
            }
        }
    }

    vertex_matrix_type output = vertex_matrix_type(0);

    // check for non-empty lists and create new graph
    for (int i = 0; i < size; i++) {
        int sub_graph_size = vertex_lists.getListSize(i);
        if (sub_graph_size != 0) {
            output.insert(output.end(), vertex_lists.getList(i));
        }
    }

    return output;
}

template<typename Type>
void EdgeFileParser::extendMatrix(std::vector<std::vector<Type> > &matrix, Type init_value) {
    int size = matrix.size();
    // new line in weight matrix will be added
    std::vector<Type> newv = std::vector<Type>(size, init_value);
    matrix.insert(matrix.end(), newv);
}


void
EdgeFileParser::initGraphSetFromXColumnFile(char *fname, GraphSet &graph_set, double th, costs_parsing_fct_type fct,
                                            int x) {
    // create associated map for vertix names, name + its index in the main graph
    std::map<std::string, int> vertices = std::map<std::string, int>();
    // create reverse as vector to get name by vertex index
    std::vector<std::string> names = std::vector<std::string>();

    // saves to name of problem vertex pair their scores and their indices
    std::map<std::string, problem> problems;

    // open file
    std::ifstream graph_file(fname);
    if (!graph_file) {
        throw GraphException("could not open file ");
    }

    char_matrix_type heuristic = char_matrix_type(0);

    // read in file line-wise
    std::string row;
    while (getline(graph_file, row)) {
        if (strtok(strdup(row.c_str()), "\t") == NULL) {
            continue;
        }

        // get line content
        std::string gene_name_1;
        std::string gene_name_2;
        double scr;

        switch (x) {
            case 3 :
                EdgeFileParser::get3ColumnContent(row, gene_name_1, gene_name_2, scr);
                break;
            case 5 :
                EdgeFileParser::get5ColumnContent(row, gene_name_1, gene_name_2, scr);
                break;
            case 12 :
                EdgeFileParser::get12ColumnContent(row, gene_name_1, gene_name_2, scr);
                break;
        }

        double costs = fct(scr, th);

        // to read in the data from Bielefeld
        if (gene_name_1 == gene_name_2) continue;

        // get positions from vertix map and reduce by 1
        int pos_gene_1 = vertices[gene_name_1] - 1;
        int pos_gene_2 = vertices[gene_name_2] - 1;

        if (pos_gene_1 != -1 && pos_gene_2 != -1) {
            if (pos_gene_1 < pos_gene_2) {
                swap<int>(pos_gene_1, pos_gene_2);
            }

            if (heuristic[pos_gene_1][pos_gene_2] != ' ') {
                if ((costs > 0 && heuristic[pos_gene_1][pos_gene_2] != '+') ||
                    (costs <= 0 && heuristic[pos_gene_1][pos_gene_2] != '-')) {
                    heuristic[pos_gene_1][pos_gene_2] = '0';
                    // save problematic pair in problem list
                    if (gene_name_1.compare(gene_name_2) < 0) {
                        swap<std::string>(gene_name_1, gene_name_2);
                        //swap<int>(pos_gene_1, pos_gene_2);
                    }
                    if (problems.find((gene_name_1 + gene_name_2)) == problems.end()) {
                        problem p = {0.0, 0.0, pos_gene_1, pos_gene_2};
                        problems[(gene_name_1 + gene_name_2)] = p;
                    }
                }
            } else {
                heuristic[pos_gene_1][pos_gene_2] = (costs > 0) ? '+' : '-';
            }
        } else if (pos_gene_1 != -1 && pos_gene_2 == -1) {
            extendMatrix<char>(heuristic, ' ');
            vertices[gene_name_2] = heuristic.size();
            names.insert(names.end(), gene_name_2);
            heuristic[heuristic.size() - 1][pos_gene_1] = (costs > 0) ? '+' : '-';
        } else if (pos_gene_1 == -1 && pos_gene_2 != -1) {
            extendMatrix<char>(heuristic, ' ');
            vertices[gene_name_1] = heuristic.size();
            names.insert(names.end(), gene_name_1);
            heuristic[heuristic.size() - 1][pos_gene_2] = (costs > 0) ? '+' : '-';
        } else {
            extendMatrix<char>(heuristic, ' ');
            extendMatrix<char>(heuristic, ' ');
            vertices[gene_name_1] = heuristic.size() - 1;
            names.insert(names.end(), gene_name_1);
            vertices[gene_name_2] = heuristic.size();
            names.insert(names.end(), gene_name_2);
            heuristic[heuristic.size() - 1][heuristic.size() - 2] = (costs > 0) ? '+' : '-';
        }
    }

    graph_file.close();

    // ++++++++++++ update problems

    // open file
    std::ifstream graph_file2(fname);
    if (!graph_file2) {
        throw GraphException("could not open file ");
    }

    while (getline(graph_file2, row)) {
        if (strtok(strdup(row.c_str()), "\t") == NULL) {
            continue;
        }

        // get line content
        std::string gene_name_1;
        std::string gene_name_2;
        double scr;

        switch (x) {
            case 3 :
                EdgeFileParser::get3ColumnContent(row, gene_name_1, gene_name_2, scr);
                break;
            case 5 :
                EdgeFileParser::get5ColumnContent(row, gene_name_1, gene_name_2, scr);
                break;
            case 12 :
                EdgeFileParser::get12ColumnContent(row, gene_name_1, gene_name_2, scr);
                break;
        }

        double costs = fct(scr, th);

        bool did_swap = false;
        if (gene_name_1.compare(gene_name_2) < 0) {
            swap<std::string>(gene_name_1, gene_name_2);
            did_swap = true;
        }

        if (problems.find((gene_name_1 + gene_name_2)) != problems.end()) {
            problem p = problems[(gene_name_1 + gene_name_2)];
            // save maximum
            if (did_swap) {
                p.backward = (p.backward > costs) ? p.backward : costs;
            } else {
                p.forward = (p.forward > costs) ? p.forward : costs;
            }
            problems[(gene_name_1 + gene_name_2)] = p;
            int i1 = p.index1;
            int i2 = p.index2;
            /*if (i1 < i2) {
                 swap<int>(i1,i2);
            }*/
            if (p.forward > p.backward) {
                heuristic[i1][i2] = (p.backward > 0) ? '+' : '-';
            } else {
                heuristic[i1][i2] = (p.forward > 0) ? '+' : '-';
            }

        }
    }

    graph_file2.close();

    // calculate connected components
    vertex_matrix_type cc = getConnectedVertices(heuristic);

    // create second graph set for backward scores
    GraphSet graph_set2;

    // create graph sets for forward and backward direction
    // + create vector with index in the main graph as key with
    // the information of their component and internal index
    // in their component

    // its index and its component number
    std::vector<cc_element> component = std::vector<cc_element>(names.size());

    for (int i = 0; i < cc.size(); i++) {
        CostsGraph graph = CostsGraph(cc[i].size(), th);

        for (int k = 0; k < cc[i].size(); k++) {
            cc_element e = {i, k};
            component[cc[i][k]] = e;
            graph.setVertexName(k, names[cc[i][k]]);
        }

        graph_set.addGraph(graph);
        graph_set2.addGraph(graph);
    }

    double forbidden = ((graph_set).getGraph(0)).forbidden;
    double permanent = ((graph_set).getGraph(0)).permanent;

    // ++++++++++++ create graph set

    // open file
    std::ifstream graph_file3(fname);
    if (!graph_file3) {
        throw GraphException("could not open file ");
    }

    while (getline(graph_file3, row)) {
        if (strtok(strdup(row.c_str()), "\t") == NULL) {
            continue;
        }

        // get line content
        std::string gene_name_1;
        std::string gene_name_2;
        double scr;

        switch (x) {
            case 3 :
                EdgeFileParser::get3ColumnContent(row, gene_name_1, gene_name_2, scr);
                break;
            case 5 :
                EdgeFileParser::get5ColumnContent(row, gene_name_1, gene_name_2, scr);
                break;
            case 12 :
                EdgeFileParser::get12ColumnContent(row, gene_name_1, gene_name_2, scr);
                break;
        }

        double costs = fct(scr, th);

        // to read in the data from Bielefeld
        if (gene_name_1 == gene_name_2) continue;

        // get positions from vertix map and reduce by 1
        int pos_gene_1 = vertices[gene_name_1] - 1;
        int pos_gene_2 = vertices[gene_name_2] - 1;
        cc_element cc_info1 = component[pos_gene_1];
        cc_element cc_info2 = component[pos_gene_2];

        if (cc_info1.component == cc_info2.component) {

            if (pos_gene_1 < pos_gene_2) {
                // if yes, get already existing score
                double old_costs = (graph_set.getGraph(cc_info1.component)).getEdge(cc_info1.index, cc_info2.index);

                // if score is already defined take average of new and old score
                if (old_costs != forbidden) {
                    (graph_set.getGraph(cc_info1.component)).setEdge(cc_info1.index, cc_info2.index,
                                                                     (costs > old_costs) ? costs : old_costs);
                } else {
                    (graph_set.getGraph(cc_info1.component)).setEdge(cc_info1.index, cc_info2.index, costs);
                }
            } else {
                // if yes, get already existing score
                double old_costs = (graph_set2.getGraph(cc_info1.component)).getEdge(cc_info1.index, cc_info2.index);

                // if score is already defined take average of new and old score
                if (old_costs != forbidden) {
                    (graph_set2.getGraph(cc_info1.component)).setEdge(cc_info1.index, cc_info2.index,
                                                                      (costs > old_costs) ? costs : old_costs);
                } else {
                    (graph_set2.getGraph(cc_info1.component)).setEdge(cc_info1.index, cc_info2.index, costs);
                }

            }
        }

    }

    graph_file3.close();

    // get min of backward and forward scores
    for (int i = 0; i < graph_set.getSetSize(); i++) {
        CostsGraph &graph = graph_set.getGraph(i);
        CostsGraph &graph2 = graph_set2.getGraph(i);
        for (int k = 0; k < (graph).getSize(); k++) {
            for (int h = 0; h < k; h++) {
                if ((graph2).getEdge(k, h) != forbidden && (graph2).getEdge(k, h) < (graph).getEdge(k, h)) {
                    (graph).setEdge(k, h, (graph2).getEdge(k, h));
                }
                if (graph.getEdge(k, h) == forbidden) {
                    graph.setEdge(k, h, forbidden + 1);
                } else if (graph.getEdge(k, h) == permanent) {
                    graph.setEdge(k, h, permanent - 1);
                };
            }
        }

    }

}


void EdgeFileParser::initFromXColumnFile(char *fname, CostsGraph &G, costs_parsing_fct_type fct, int x) {
    double th = G.getThreshold();
    double forbidden = G.forbidden;


    double_matrix_type costs_matrix = double_matrix_type(0);
    double_matrix_type costs_matrix2 = double_matrix_type(0);

    // open file
    std::ifstream graph_file(fname);
    if (!graph_file) {
        throw GraphException("could not open file ");
    }

    // create associated map for vertix names
    std::map<std::string, int> *vertices = new std::map<std::string, int>();

    CostsGraph::vertex_name_type names = CostsGraph::vertex_name_type(0);

    // read in file line-wise
    std::string row;
    while (getline(graph_file, row)) {

        if (strtok(strdup(row.c_str()), "\t") == NULL) {
            continue;
        }

        // get line content
        std::string gene_name_1;
        std::string gene_name_2;
        double scr;

        switch (x) {
            case 3 :
                EdgeFileParser::get3ColumnContent(row, gene_name_1, gene_name_2, scr);
                break;
            case 5 :
                EdgeFileParser::get5ColumnContent(row, gene_name_1, gene_name_2, scr);
                break;
            case 12 :
                EdgeFileParser::get12ColumnContent(row, gene_name_1, gene_name_2, scr);
                break;
        }

        double costs = fct(scr, th);

        // to read in the data from Bielefeld
        if (gene_name_1 == gene_name_2) continue;

        // get positions from vertix map and reduce by 1
        int pos_gene_1 = (*vertices)[gene_name_1] - 1;
        int pos_gene_2 = (*vertices)[gene_name_2] - 1;

        // check if vertices is already defined
        if (pos_gene_1 != -1 && pos_gene_2 != -1) {

            if (pos_gene_1 < pos_gene_2) {
                // if yes, get already existing score
                double old_costs = costs_matrix[pos_gene_1][pos_gene_2];

                // if score is already defined take average of new and old score
                if (old_costs != forbidden) {
                    costs_matrix[pos_gene_1][pos_gene_2] = (costs > old_costs) ? costs : old_costs;
                } else {
                    costs_matrix[pos_gene_1][pos_gene_2] = costs;
                }
            } else {
                // if yes, get already existing score
                double old_costs = costs_matrix2[pos_gene_1][pos_gene_2];

                // if score is already defined take average of new and old score
                if (old_costs != forbidden) {
                    costs_matrix2[pos_gene_1][pos_gene_2] = (costs > old_costs) ? costs : old_costs;
                } else {
                    costs_matrix2[pos_gene_1][pos_gene_2] = costs;
                }

            }

        } else if (pos_gene_1 == -1 && pos_gene_2 != -1) {

            // if vertex 1 is not defined add vertex
            extendMatrix<double>(costs_matrix, forbidden);
            extendMatrix<double>(costs_matrix2, forbidden);
            names.insert(names.end(), gene_name_1);

            // save position in associated map
            (*vertices)[gene_name_1] = G.getSize();

            // set edge weight from file
            costs_matrix[costs_matrix.size() - 1][pos_gene_2] = forbidden;
            costs_matrix2[costs_matrix2.size() - 1][pos_gene_2] = costs;

        } else if (pos_gene_1 != -1 && pos_gene_2 == -1) {

            // if vertex 2 is not defined add vertex
            extendMatrix<double>(costs_matrix, forbidden);
            extendMatrix<double>(costs_matrix2, forbidden);
            names.insert(names.end(), gene_name_2);

            // save position in associated map
            (*vertices)[gene_name_2] = G.getSize();

            // set edge weight from file
            costs_matrix[costs_matrix.size() - 1][pos_gene_1] = costs;
            costs_matrix2[costs_matrix2.size() - 1][pos_gene_1] = forbidden;

        } else {

            // if neither of them is defined, add both
            extendMatrix<double>(costs_matrix, forbidden);
            extendMatrix<double>(costs_matrix2, forbidden);
            extendMatrix<double>(costs_matrix, forbidden);
            extendMatrix<double>(costs_matrix2, forbidden);

            names.insert(names.end(), gene_name_1);
            names.insert(names.end(), gene_name_2);

            // save positions in associated map
            (*vertices)[gene_name_1] = G.getSize() - 1;
            (*vertices)[gene_name_2] = G.getSize();

            // set edge weight from file
            costs_matrix[costs_matrix.size() - 1][costs_matrix.size() - 2] = costs;
            costs_matrix2[costs_matrix2.size() - 1][costs_matrix2.size() - 2] = forbidden;

        };
    }

    for (int i = 0; i < costs_matrix.size(); i++) {
        for (int k = 0; k < i; k++) {
            if (costs_matrix2[i][k] != forbidden && costs_matrix2[i][k] < costs_matrix[i][k]) {
                costs_matrix[i][k] = costs_matrix2[i][k];
            }
        }
    }

    CostsGraph new_graph(costs_matrix.size(), costs_matrix, names, th);

    G = new_graph;

    // close input file and delete vertix map
    graph_file.close();
    delete (vertices);

}


void EdgeFileParser::to5ColumnFile(char *fname, CostsGraph &G, double *matrix) {
    to5ColumnFile<double>(fname, G, matrix);
}

void EdgeFileParser::to5ColumnFile(char *fname, CostsGraph &G, int *matrix) {
    to5ColumnFile<int>(fname, G, matrix);
}

void EdgeFileParser::to5ColumnFile(char *fname, CostsGraph &G, short *matrix) {
    to5ColumnFile<short>(fname, G, matrix);
}

void EdgeFileParser::to5ColumnFile(char *fname, CostsGraph &G, char *matrix) {
    to5ColumnFile<char>(fname, G, matrix);
}

void EdgeFileParser::initFrom5ColumnFile(char *fname, CostsGraph &G, costs_parsing_fct_type fct) {
    initFromXColumnFile(fname, G, fct, 5);
}

void EdgeFileParser::initFrom3ColumnFile(char *fname, CostsGraph &G, costs_parsing_fct_type fct) {
    initFromXColumnFile(fname, G, fct, 3);
}

void EdgeFileParser::initFrom12ColumnFile(char *fname, CostsGraph &G, costs_parsing_fct_type fct) {
    initFromXColumnFile(fname, G, fct, 12);
}

void
EdgeFileParser::initGraphSetFrom3ColumnFile(char *fname, GraphSet &graph_set, double th, costs_parsing_fct_type fct) {
    std::cout << "call init from 3 column file " << std::endl;
    initGraphSetFromXColumnFile(fname, graph_set, th, fct, 3);
}

void
EdgeFileParser::initGraphSetFrom5ColumnFile(char *fname, GraphSet &graph_set, double th, costs_parsing_fct_type fct) {
    initGraphSetFromXColumnFile(fname, graph_set, th, fct, 5);
}

void
EdgeFileParser::initGraphSetFrom12ColumnFile(char *fname, GraphSet &graph_set, double th, costs_parsing_fct_type fct) {
    initGraphSetFromXColumnFile(fname, graph_set, th, fct, 12);
}
