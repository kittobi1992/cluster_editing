#include <fstream>
#include <iostream>
#include <sstream>

#include <matrixparser.h>
#include <costsgraph.h>
#include <math.h>
#include <iomanip>
#include <string.h>

void fillParameters(int argc, char **argv, std::string &fname, char *&output_fname, bool &names);

void readOutputFile(char *output_name, CostsGraph::vertex_matrix_type &cc, CostsGraph &G);

int getComp(int i, CostsGraph::vertex_matrix_type &cc);

std::string toString(int x);

int main(int argc, char **argv) {
    // initialize parameters with standard settings
    std::string fname = "";
    char *output_fname = "";
    bool names = true;
    CostsGraph::matrix_file_fct_type matrix_init_fct = MatrixParser::initFromMatrixFile;

    // fill parameters with values from shell
    fillParameters(argc, argv, fname, output_fname, names);

    try {
        // create costs graph by reading from file
        CostsGraph G = CostsGraph(fname, matrix_init_fct);

        CostsGraph::vertex_matrix_type cc = CostsGraph::vertex_matrix_type(0);

        // create file
        std::ofstream graph_file(output_fname);
        if (!graph_file) {
            throw GraphException("could not create file ");
        }


        // get max weight
        double max = 0.0;
        double min = 0.0;
        for (int i = 0; i < G.getSize(); i++) {
            for (int k = 0; k < i; k++) {
                max = (G.getEdge(i, k) > max) ? G.getEdge(i, k) : max;
                min = (G.getEdge(i, k) < min) ? G.getEdge(i, k) : min;
            }
        }
        bool n = names;


        // write graph file from graph
        graph_file << "Creator \"yFiles\"" << std::endl;
        graph_file << "Version \"2.4.2.2\"" << std::endl;
        graph_file << "graph" << std::endl;
        graph_file << "[" << std::endl;
        graph_file << "        hierarchic      1" << std::endl;
        graph_file << "        label   \"\"" << std::endl;
        graph_file << "        directed        0" << std::endl;
        // write nodes
        for (int i = 0; i < G.getSize(); i++) {
            graph_file << "\tnode" << std::endl;
            graph_file << "\t[" << std::endl;
            graph_file << "\t\tid\t" << i << std::endl;
            if (names) {
                graph_file << "\t\tlabel\t\"" << G.getVertexName(i) << "\"" << std::endl;
            } else {
                graph_file << "\t\tlabel\t\"" << i << "\"" << std::endl;
            }
            graph_file << "\t]" << std::endl;
        }

        for (int i = 0; i < G.getSize(); i++) {
            for (int k = 0; k < i; k++) {
                if (G.getEdge(i, k) > 0) {
                    char buffer[2];
                    int n = 200 - ceil(200 * (G.getEdge(i, k) / max));
                    //std::string s = printf("%x", n);
                    sprintf(buffer, "%2.0x", n);
                    for (int j = 0; j < 2; j++) {
                        if (buffer[j] == ' ') {
                            buffer[j] = '0';
                        }
                    }
                    graph_file << "\tedge" << std::endl;
                    graph_file << "\t[" << std::endl;
                    graph_file << "\t\tsource " << i << std::endl;
                    graph_file << "\t\ttarget " << k << std::endl;
                    graph_file << "\t\tgraphics" << std::endl;
                    graph_file << "\t\t[" << std::endl;
                    graph_file << "\t\t\tfill\t\"#" << buffer << buffer << buffer << "\"" << std::endl;
                    graph_file << "\t\t]" << std::endl;
                    graph_file << "\t]" << std::endl;
                } else {
                    /*if (G.getEdge(i,k)/min < 0.1) {
                         graph_file << "\tedge" << std::endl;
                         graph_file << "\t[" << std::endl;
                         graph_file << "\t\tsource " << i << std::endl;
                         graph_file << "\t\ttarget " << k << std::endl;
                         graph_file << "\t\tgraphics" << std::endl;
                         graph_file << "\t\t[" << std::endl;
                         graph_file << "\t\t\tstyle\t\"dashed\"" << std::endl;
                         graph_file << "\t\t\tfill\t\"#d7d7d7\"" << std::endl;
                         graph_file << "\t\t]" << std::endl;
                         graph_file << "\t]" << std::endl;
                    }*/
                }
            }
        }
        graph_file << "]" << std::endl;
        graph_file.close();

    } catch (GraphException e) {
        std::cout << "Graph Exception: " << e.getMessage() << std::endl;
    }

    return 0;
}

void readOutputFile(char *output_fname, CostsGraph::vertex_matrix_type &cc, CostsGraph &G) {
    // open output file, which is the second file
    std::ifstream output_file(output_fname);
    if (!output_file) {
        throw GraphException("could not open file ");
    }

    // save if "Component" was already once found
    bool first_occurence_of_component = false;

    // container for each row
    std::string row;

    //skip call line, which is the first line in the output file
    getline(output_file, row);

    // create matrix to save connected component indices
    cc = CostsGraph::vertex_matrix_type(0);

    // counter for number of component
    int comp = 0;

    // iterate through lines
    while (getline(output_file, row)) {
        // if line says not "Component" skip it
        if (strstr(strdup(row.c_str()), "Component") == NULL) {
            continue;
        } else if (!first_occurence_of_component) {
            // if this is the first occurence of "Component" skip it
            first_occurence_of_component = true;
            continue;
        }

        // create container for connected component
        cc.insert(cc.end(), CostsGraph::vertex_set_type(0));

        // get position of colon
        char *str = strdup(row.c_str());
        char *p_to_colon;
        p_to_colon = strstr(str, ":");
        if (p_to_colon == NULL) {
            continue;
        } else {
            // read elements after position of colon
            char *p_to_elem;
            p_to_elem = strtok(p_to_colon + 2, " ");

            // iterate through elements
            while (p_to_elem != NULL) {
                cc[comp].insert(cc[comp].end(), G.getIndex(p_to_elem));
                p_to_elem = strtok(NULL, " ");
            }
        }
        comp++;
    }
}


int getComp(int i, CostsGraph::vertex_matrix_type &cc) {
    for (int j = 0; j < cc.size(); j++) {
        for (int k = 0; k < cc[j].size(); k++) {
            if (cc[j][k] == i) return j;
        }
    }
}


std::string toString(int x) {
    std::ostringstream o;
    o << x;
    return o.str();
}


void fillParameters(int argc, char **argv, std::string &fname, char *&output_fname, bool &names) {
    int g = 1;
    while (argc > g) {
        std::string flag(argv[g]);
        if (flag == "--help") {
            std::cout << "GraphViz File Creator" << std::endl;
            std::cout << "---------------------" << std::endl << std::endl;
            std::cout << "Usage: graphvizcreator [cm-file] [output-file] [0=number,1=names]" << std::endl << std::endl;
            std::cout << "Creates a file for the graph visualization\n software yEd from the given\n *.cm file. "
                      << std::endl << std::endl;
            exit(1);
        } else if (g == 1) {
            fname = argv[g];
            g++;
        } else if (g == 2) {
            output_fname = argv[g];
            g++;
        } else if (g == 3) {
            std::istringstream i(argv[g]);
            int help;
            i >> help;
            if (help == 0) {
                names = false;
            } else {
                names = true;
            }
            return;
        }

    }
}

