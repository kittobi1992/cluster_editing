#include <iostream>
#include <sstream>

#include <searchtree.h>
#include <probleminstanceeditoperations.h>
#include <probleminstanceeditcosts.h>
#include <edgefileparser.h>
#include <blastparser.h>
#include <graphset.h>

int main(int argc, char **argv) {
    if (argc != 4) {
        if (argc < 2) {
            std::cerr << "No input file!" << std::endl;
        } else {
            if (argc < 3) {
                std::cerr << "No threshold given!" << std::endl;
            } else {
                std::cerr << "Too many parameters!" << std::endl;
            }
        }
        exit(1);
    }

    std::string fname = argv[1];
    std::istringstream i(argv[2]);
    double th;
    i >> th;
    std::istringstream k(argv[3]);
    int pi_type;
    k >> pi_type;
    if (pi_type != 1 && pi_type != 2) {
        std::cerr << "Unknown problem instance type!" << std::endl;
        exit(1);
    }

    std::cout << std::endl << "test with " << fname << " file" << std::endl;
    try {
        std::cout << "read graph file..." << std::endl;
        CostsGraph H(fname, EdgeFileParser::initFrom5ColumnFile, BlastParser::exponentToCosts, th);
        std::cout << "graph size: " << H.getSize() << std::endl;
        std::cout << "create graph set by dividing connected components" << std::endl;
        GraphSet graph_set(H);
        std::cout << "--> in sum " << graph_set.getSetSize() << " connected components" << std::endl;

        int min_parameter = 0;
        long solutions_nr = 0;

        for (int i = 0; i < graph_set.getSetSize(); i++) {
            CostsGraph &G = graph_set.getGraph(i);
            std::cout << "reduction with graph " << i << " of size " << G.getSize() << " ";

            ProblemInstance *PI;
            int done = 0;
            int par = 1;
            int par_add = 1;
            if (pi_type == 2) {
                par_add = 1;
            }
            while (done == 0) {
                try {
                    std::cout << "try with parameter = " << par << std::endl;
                    par += par_add;

                    if (pi_type == 1) {
                        PI = (ProblemInstance_EditOperations(G, par - par_add)).clone();
                    } else {
                        PI = (ProblemInstance_EditCosts(G, par - par_add)).clone();
                    }

                    std::cout << G << std::endl;
                    int count = (*PI).reduce();
                    std::cout << " with parameter " << par - par_add << " caused " << count
                              << " changes. Parameter was reduced by " << par - (*PI).getParameter() - par_add
                              << std::endl;

                    // SEARCH TREE !!!
                    std::cout << "create SearchTree..." << std::endl;
                    SearchTree ST(PI);
                    std::cout << "search..." << std::endl;
                    ST.search();
                    SearchTree::solutions_type solutions = ST.getSolutions();
                    std::cout << "found " << solutions.size() << " solutions!" << std::endl;
                    if (solutions.size() > 0) {
                        std::cout << "solutions: " << std::endl;
                        for (int k = 0; k < solutions.size(); k++) {
                            std::cout << (*solutions[k]).getGraph() << std::endl;
                            std::cout << "changed: ";
                            ProblemInstance::edge_list_type changed = (*solutions[k]).getChangedEdges();
                            for (int i = 0; i < changed.size(); i++) {
                                std::cout << "(" << changed[i].i << "," << changed[i].j << ") costs: "
                                          << changed[i].costs << " ";
                            }
                            std::cout << std::endl;
                            std::cout << "with costs"
                                      << (*solutions[k]).getStartParameter() - (*solutions[k]).getParameter() - 1
                                      << std::endl;
                        }
                        solutions_nr = solutions_nr + solutions.size();
                        min_parameter += par - par_add;
                        done = 1;
                    }
                } catch (ProblemInstanceException e) {
                    std::cout << "PK-Exception: " << e.getMessage() << std::endl;
                } catch (GraphException e) {
                    std::cout << "Graph Exception: " << e.getMessage() << std::endl;
                }
            }
            delete (PI);
            std::cout << "-----------------------" << std::endl;
        }
        std::cout << "################" << std::endl;
        std::cout << "orignal graph size: " << H.getSize() << std::endl;
        std::cout << "statistic of created graph set: " << std::endl;
        graph_set.printStat();
        std::cout << "REDUCTION: " << std::endl;
        std::cout << "minimal parameter of problem kernel reductions: " << min_parameter << std::endl;
        std::cout << "number of solutions (added): " << solutions_nr << std::endl;

    } catch (EdgeListsException e) {
        std::cout << "EdgeLists-Exception: " << e.getMessage() << std::endl;
    } catch (VertexListsException e) {
        std::cout << "VertexLists-Exception: " << e.getMessage() << std::endl;
    } catch (ProblemInstanceException e) {
        std::cout << "PK-Exception: " << e.getMessage() << std::endl;
    } catch (GraphException e) {
        std::cout << "Graph Exception: " << e.getMessage() << std::endl;
    }

    return 0;
}
