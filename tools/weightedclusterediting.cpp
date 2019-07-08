#include <fstream>
#include <iostream>
#include <sstream>

#include <string.h>

#include <searchtreeweighted.h>
#include <weightedprobleminstance.h>
#include <edgefileparser.h>
#include <matrixparser.h>

#include <stopwatch.h>

#include <math.h>

#include <blastparser.h>
#include <doubleparser.h>
#include <graphset.h>

void fillParameters(int argc, char **argv, std::string &fname, char *&output_fname,
                    GraphSet::costs_parsing_fct_type &parsing_fct, GraphSet::graph_set_parser_fct_type &read_graph_fct,
                    double &threshold, double &parameter, int &parameter_steps, short &pi_type, int &max_component,
                    bool &merging, short &mode, bool &just_max, bool &split, bool &loginfo, char *&log_fname);

double solveForConnectedComponent(CostsGraph G, int nr, double &min_parameter, long &solutions_nr, double parameter,
                                  int parameter_steps, std::ofstream &log_file, short pi_type, bool merging,
                                  bool just_max, bool split, double &costs, std::ofstream &output_file, bool &loginfo);

void writeVertexList(CostsGraph::index_list_type &vertex_list, std::ofstream &log_file, bool &loginfo);

void writeObviousComponent(int nr, CostsGraph &G, std::ofstream &log_file, bool &loginfo);

void writeHeaderOfSolution(int nr, CostsGraph &graph, std::ofstream &log_file, bool &loginfo);

void writeCalculatedSolutions(int nr, SearchTreeWeighted::solutions_type solutions, std::ofstream &log_file,
                              std::ofstream &output_file, bool &loginfo);

int main(int argc, char **argv) {
    // initialize parameters with standard settings
    std::string fname = "";
    char *output_fname = "";
    char *log_fname = "";
    GraphSet::costs_parsing_fct_type parsing_fct = parsing_fct = DoubleParser::valueToCosts;
    GraphSet::graph_set_parser_fct_type read_graph_fct = EdgeFileParser::initGraphSetFrom5ColumnFile;
    GraphSet::matrix_file_fct_type matrix_file_fct = MatrixParser::initFromMatrixFile;
    double threshold = 1e-40;
    double parameter = 0;
    int parameter_steps = 0;
    short pi_type = 2;
    bool merging = true;
    bool just_max = true;
    bool split = true;
    bool loginfo = false;
    int max_component = 32767;

    short mode = 1;

    double used_time = 0;

    // fill parameters with values from shell
    fillParameters(argc, argv, fname, output_fname, parsing_fct, read_graph_fct, threshold, parameter, parameter_steps,
                   pi_type, max_component, merging, mode, just_max, split, loginfo, log_fname);

    if (strcmp("", log_fname) != 0) { //checking if log_fname is empty to decide logging or not
        loginfo = true;
    }/*else{
		log_fname="delete_me.temp";
	}	*/

    // open output file for writing
    std::ofstream output_file(output_fname);
    if (!output_file) {
        std::cout << "Can not open output file!" << std::endl;
        exit(0);
    }


    std::ofstream log_file(log_fname);

    if (loginfo) {
        std::cout << "logging is active" << std::endl;
        // open logging file for writing
        //std::ofstream log_file(log_fname);
        if (!log_file) {
            std::cout << "Can not create logging file" << std::endl;
            exit(0);
        }
    } else {
        std::cout << "logging is off, nothing is logged" << std::endl;

    }

    // write all parameters to file
    for (int i = 0; i < argc; i++) {
        if (loginfo)log_file << argv[i] << " ";
    }
    if (loginfo)log_file << std::endl << std::endl;

    try {
        // create main costs graph by reading from file
        std::cout << "Read graph file..." << std::endl;

        GraphSet graph_set;

        if (mode == 1) {
            graph_set = GraphSet(fname, threshold, parsing_fct, read_graph_fct);
        } else if (mode == 3) {
            graph_set = GraphSet(fname, matrix_file_fct);
        } else if (mode == 2) {
            graph_set = GraphSet(CostsGraph(fname, matrix_file_fct));
        }

        graph_set.printStat();
        // get rid of main graph
        //delete main_graph;

        double min_parameter = 0;
        long solutions_nr = 0;

        // solve problem for
        for (int i = 0; i < graph_set.getSetSize(); i++) {
            CostsGraph &G = graph_set.getGraph(i);

            if (G.getSize() <= max_component) {

                double x = min_parameter;
                double costs = 0;
                //double this_time =0;
                double this_time = solveForConnectedComponent(G, i, min_parameter, solutions_nr, parameter,
                                                              parameter_steps, log_file, pi_type, merging, just_max,
                                                              split, costs, output_file, loginfo);
                used_time += this_time;
                std::cout << "STATISTICS" << std::endl;
                std::cout << "nodes: " << G.getSize() << std::endl;
                std::cout << "edges: " << G.getEdgeNumber() << std::endl;
                std::cout << "time: " << this_time << std::endl;
                std::cout << "costs: " << costs << std::endl;
                std::cout << "-----------------------" << std::endl;
                // test...include heuristic
                //WeightedProblemInstance PI = WeightedProblemInstance(G, 1E5);
                //std::cout << "start heuristic..." << std::endl;
                //PI.heuristicSolve();
                //std::cout << "heuristic costs = " << (1E5 - PI.getParameter()) << std::endl;
                //SearchTreeWeighted ST(PI, 0);
                //ST.search(just_max, split);
                //SearchTreeWeighted::solutions_type solutions = ST.getSolutions();
                //writeCalculatedSolutions(0, solutions, output_file);
            } else {
                output_file << "Component " << i << std::endl;
                //std::cout << "STEFAN IST COOL" << std::endl;
                // output_file << "> " << max_component << std::endl << std::endl;
                std::cout << "STATISTICS" << std::endl;
                std::cout << "nodes: " << G.getSize() << std::endl;
                std::cout << "edges: " << G.getEdgeNumber() << std::endl;
            }
        }

        output_file << "Sum of costs: " << min_parameter << std::endl;

        std::cout << "################" << std::endl;
        //std::cout << "orignal graph size: " << H.getSize() << std::endl;
        std::cout << "statistic of created graph set: " << std::endl;
        graph_set.printStat();
        std::cout << "REDUCTION: " << std::endl;
        std::cout << "minimal cost for cluster editing: " << min_parameter << std::endl;
        std::cout << "number of solutions (added): " << solutions_nr << std::endl;
        std::cout << "time needed: " << used_time << std::endl;

    } catch (VertexListsException &e) {
        std::cout << "VertexLists-Exception: " << e.getMessage() << std::endl;
    } catch (ProblemInstanceException &e) {
        if (parameter != 0) {
            std::cout << "ProblemInstance-Exception :" << std::endl;
            std::cout << "=> Given problem with parameter " << parameter << " insolvable." << std::endl;
        }
    } catch (GraphException &e) {
        std::cout << "Graph Exception: " << e.getMessage() << std::endl;
    }
    if (loginfo)log_file.close();
    /*if(!loginfo){
        output_file << std::endl;
        output_file << std::endl;
           std::ifstream f;  // Datei-Handle
          string s;
          f.open(log_fname, std::ios::in); // oeffne Datei aus Parameter
          while (!f.eof())          // Solange noch Daten vorliegen
          {
              getline(f, s);        // Lese eine Zeile
            output_file << s << std::endl;
          }
    }*/
    output_file.close();


    return 0;
}


void writeObviousComponent(int nr, CostsGraph &G, std::ofstream &log_file, bool &loginfo) {
    writeHeaderOfSolution(nr, G, log_file, loginfo);

    if (loginfo)log_file << "#Solutions: 1" << std::endl;
    if (loginfo)log_file << "Solution 1: " << std::endl;
    if (loginfo)log_file << "Changes: none" << std::endl << "Total cost: 0" << std::endl;
    if (loginfo)log_file << "Component 1: " << G.vertexNamesToString() << std::endl;
    if (loginfo)log_file << std::endl;
}

void writeHeaderOfSolution(int nr, CostsGraph &graph, std::ofstream &log_file, bool &loginfo) {
    if (loginfo)log_file << "Component " << nr << std::endl;
    if (loginfo)log_file << "Size " << graph.getSize() << std::endl;
    if (loginfo)log_file << "Vertices: " << graph.vertexNamesToString() << std::endl;
}

void writeVertexList(CostsGraph::index_list_type &vertex_list, std::ofstream &log_file, bool &loginfo) {
    if (!vertex_list.empty()) {
        if (loginfo)log_file << vertex_list[0];
        for (int i = 1; i < vertex_list.size(); i++) {
            if (loginfo)log_file << "," << vertex_list[i];
        }
    }
}

void writeCalculatedSolutions(int nr, SearchTreeWeighted::solutions_type solutions, std::ofstream &log_file,
                              std::ofstream &output_file, bool &loginfo) {

    if (loginfo)log_file << "#Solutions: " << solutions.size() << std::endl;
    for (int k = 0; k < solutions.size(); k++) {
        if (loginfo)log_file << "Solution " << k << ": " << std::endl;
        if (loginfo)log_file << "Changes: ";
        WeightedProblemInstance::edge_list_type changed = (solutions[k].changes);
        if (changed.empty()) {
            if (loginfo)log_file << "none";
        } else {
            switch (changed[0].operation) {
                case 'a' :
                    if (loginfo)log_file << "add";
                    break;
                case 'd' :
                    if (loginfo)log_file << "del";
                    break;
                case 'm' :
                    if (loginfo)log_file << "merge";
                    break;
                default:
                    break;
            }
            if (loginfo)log_file << "(";
            writeVertexList(changed[0].i, log_file, loginfo);
            if (loginfo)log_file << ";";
            writeVertexList(changed[0].j, log_file, loginfo);
            if (loginfo)log_file << ") costs: " << changed[0].costs;

            for (int j = 1; j < changed.size(); j++) {
                switch (changed[j].operation) {
                    case 'a' :
                        if (loginfo)log_file << ", add";
                        break;
                    case 'd' :
                        if (loginfo)log_file << ", del";
                        break;
                    case 'm' :
                        if (loginfo)log_file << ", merge";
                        break;
                    default:
                        break;
                }
                if (loginfo)log_file << "(";
                writeVertexList(changed[j].i, log_file, loginfo);
                if (loginfo)log_file << ";";
                writeVertexList(changed[j].j, log_file, loginfo);
                if (loginfo)log_file << ") costs: " << changed[j].costs;
            }
        }
        if (loginfo)log_file << std::endl;
        if (loginfo)log_file << "Total costs: " << solutions[k].start_parameter - solutions[k].parameter << std::endl;
        for (int h = 0; h < (solutions[k].components).getSetSize(); h++) {
            output_file << "Component " << (h + 1) << ": ";
            CostsGraph solved_graph = (solutions[k].components).getGraph(h);
            for (int j = 0; j < solved_graph.getSize(); j++) {
                CostsGraph::vertex_name_type names = solved_graph.getVertexNames(j);
                for (const auto& name : names) {
                    output_file << name << " ";
                }
            }
            output_file << std::endl;
        }
    }
    if (loginfo)log_file << std::endl;
}

double solveForConnectedComponent(CostsGraph G, int nr, double &min_parameter, long &solutions_nr, double parameter,
                                  int parameter_steps, std::ofstream &log_file, short pi_type, bool merging,
                                  bool just_max, bool split, double &costs, std::ofstream &output_file, bool &loginfo) {
    int done = 0;
    double par = parameter;
    double sum = 0;
    double par_steps = 1;

    double used_time = 0;

    /*if (pi_type == 1) {

        par_steps = (parameter_steps > 0) ? parameter_steps : (G.getSize()/ 10 + 1);

    } else if (pi_type == 2) {
        for (int i=0; i < G.getSize(); i++) {
             for(int k=0; k < i;k++) {
                  double e = G.getEdge(i,k);
                  if (static_cast<int>(e/G.forbidden+0.001) != 1  && static_cast<int>(e/G.permanent+0.001) != 1) {
                       sum += std::abs(e);
                  }
             }
        }

        sum = sum / (G.getSize()*G.getSize());

        par_steps = (parameter_steps > 0) ? parameter_steps : (((G.getSize() / 5 * sum) + 1));
    }

    par_steps = par_steps*3;*/

    WeightedProblemInstance *PI = nullptr;

    writeHeaderOfSolution(nr, G, log_file, loginfo);

    if (G.getSize() <= 2) {
        done = 1;
        if (loginfo)log_file << std::endl;
        std::cout << "Component " << nr << "(size=" << G.getSize() << "): obvious solution" << std::endl;
        //writeObviousComponent(nr, G, output_file);

    }

    Stopwatch timer;
    timer.start();

    std::cout << "Component " << nr << ". Computation starts..." << std::endl;
    /*get lower bound*/
    PI = new WeightedProblemInstance(G, 1E12, false);
    par = PI->getLowerBound();
    double upper_bound = PI->getUpperBound();
    std::cout << " Lower Bound for Problem Instance: " << par << std::endl;
    std::cout << " Upper Bound for Problem Instance: " << upper_bound << std::endl;
    delete (PI);
    PI = nullptr;

    par_steps = (upper_bound - par) / 23 + 0.01;

    while (done == 0) {
        try {
            //if (parameter != 0 && min_parameter + par > parameter) throw ProblemInstanceException(">");

            std::cout << "Component " << nr << "(size=" << G.getSize() << "): try with costs = " << par << std::endl;

            par += par_steps;

            PI = new WeightedProblemInstance(G, par - par_steps);

            int count = PI->maxReduce();
            std::cout << " Reduction with maximum costs " << par - par_steps << " caused " << count
                      << " edge changes with costs " << par - PI->getParameter() - par_steps << " new graph size is "
                      << (PI->getGraph()).getSize() << std::endl;

            //std::cout << "this is a test...write reduced graph to test.cm" << std::endl;
            //(PI->getGraph()).costsIOOperation("test.cm",MatrixParser::toMatrixFile);

            // SEARCH TREE !!!
            std::cout << " Create SearchTree" << std::endl;
            SearchTreeWeighted ST(*PI, 0);
            std::cout << " Search for solution(s)..." << std::endl;
            ST.search(just_max, split);
            SearchTreeWeighted::solutions_type solutions = ST.getSolutions();
            std::cout << "Found " << solutions.size() << " solution(s)!" << std::endl << std::endl;
            if (!solutions.empty()) {

                used_time = timer.elapsed();

                writeCalculatedSolutions(nr, solutions, log_file, output_file, loginfo);

                solutions_nr = solutions_nr + solutions.size();

                double min = 100000E100;
                for (const auto &solution : solutions) {
                    double solution_costs = solution.start_parameter - solution.parameter;
                    min = std::min(solution_costs, min);
                }

                costs = min;
                min_parameter += min;
                done = 1;

                // delete solutions
                /*for (int i = 0; i < solutions.size(); i++) {
                     for (int k = 0; k < (solutions[i].components).size(); k++) {
                          delete( (solutions[i].components)[k] );
                     }
                }*/
            }
            delete (PI);
            PI = nullptr;
        } catch (ProblemInstanceException &e) {
            std::cout << "Problem Instance Exception: " << e.getMessage() << std::endl;
            delete (PI);
            PI = nullptr;
        } catch (GraphException &e) {
            std::cout << "Graph Exception: " << e.getMessage() << std::endl;
            delete (PI);
            PI = nullptr;
        }
    }

    return used_time;

}


template<typename Type>
Type convertStringTo(const std::string &s) {
    std::istringstream i(s);
    Type x;
    i >> x;
    return x;
}

void fillParameters(int argc, char **argv, std::string &fname, char *&output_fname,
                    GraphSet::costs_parsing_fct_type &parsing_fct, GraphSet::graph_set_parser_fct_type &read_graph_fct,
                    double &threshold, double &parameter, int &parameter_steps, short &pi_type, int &max_component,
                    bool &merging, short &mode, bool &just_max, bool &split, bool &loginfo, char *&log_fname) {
    int g = 1;
    int files = 0;
    while (argc > g) {
        if (argv[g][0] == '-') {
            std::string flag(argv[g]);
            flag = flag.erase(0, 1);

            if (flag == "-threshold" || flag == "T") {
                std::istringstream i(argv[g + 1]);
                i >> threshold;
            } else if (flag == "-costparser" || flag == "C") {
                std::string argument(argv[g + 1]);
                if (argument == "BlastParser.exponentToCosts") {
                    parsing_fct = BlastParser::exponentToCosts;
                } else if (argument == "BlastParser.simpleCutOff") {
                    parsing_fct = BlastParser::simpleCutOff;
                } else if (argument == "DoubleParser.simpleCutOff") {
                    parsing_fct = DoubleParser::simpleCutOff;
                } else if (argument == "DoubleParser.valueToCosts") {
                    parsing_fct = DoubleParser::valueToCosts;
                } else {
                    std::cout << "Unknown cost parser: " << argv[g + 1] << std::endl;
                    exit(0);
                }
            } else if (flag == "-graphparser" || flag == "G") {
                std::string argument(argv[g + 1]);
                if (argument == "EdgeFileParser.initFrom5ColumnFile") {
                    read_graph_fct = EdgeFileParser::initGraphSetFrom5ColumnFile;
                } else if (argument == "EdgeFileParser.initFrom3ColumnFile") {
                    read_graph_fct = EdgeFileParser::initGraphSetFrom3ColumnFile;
                } else if (argument == "EdgeFileParser.initFrom12ColumnFile") {
                    read_graph_fct = EdgeFileParser::initGraphSetFrom12ColumnFile;
                } else {
                    std::cout << "Unknown graph parser: " << argv[g + 1] << std::endl;
                    exit(0);
                }
            } else if (flag == "-parameter" || flag == "P") {
                std::istringstream i(argv[g + 1]);
                i >> parameter;
                if (parameter < 0) {
                    std::cout << parameter << " is not a valid parameter!" << std::endl;
                    exit(0);
                }
            } else if (flag == "-parameter_steps" || flag == "S") {
                parameter_steps = convertStringTo<int>(argv[g + 1]);
                //std::istringstream i(argv[g+1]);
                //i >> parameter_steps;
                if (parameter_steps < 1) {
                    std::cout << parameter_steps << " is not a valid step size!" << std::endl;
                    exit(0);
                }
            } else if (flag == "-unweighted" || flag == "U") {
                pi_type = 1;
                g--;
            } else if (flag == "-max_component" || flag == "M") {
                std::istringstream i(argv[g + 1]);
                i >> max_component;
                if (max_component < 1) {
                    std::cout << "Maximum component should be bigger then zero" << std::endl;
                    exit(0);
                }
            } else if (flag == "-all_solutions" || flag == "A") {
                just_max = false;
                g--;
            } else if (flag == "-no_split" || flag == "s") {
                split = false;
                g--;
//		     } else if (flag == "-no_log" || flag == "L") {
//		          loginfo = false;
//		          g--;
            } else if (flag == "-mode" || flag == "X") {
                std::istringstream i(argv[g + 1]);
                i >> mode;
                if (mode < 1 || mode > 3) {
                    std::cout << "Unknown mode" << std::endl;
                    exit(0);
                }
            } else if (flag == "-help") {
                std::cout << "Cluster Editing Tool " << std::endl;
                std::cout << "--------------------" << std::endl;
                std::cout << "Usage: [options] <input_file> <output_file>" << std::endl << std::endl;

                std::cout << "Options:" << std::endl << std::endl;
                std::cout
                        << " --mode OR -X <1|2|3> \tIndicates input mode\n\t\t\t\t1 for edge files\n\t\t\t\t2 for matrix files\n\t\t\t\t3 for parsing a directory with connected \n\t\t\t\tcomponents in matrix format\n\t\t\t\tDefault: 1"
                        << std::endl;
                std::cout << " --threshold OR -T <REAL_VALUE>	Indicates threshold for edges" << std::endl
                          << "\t\t\t\tDefault: 1e-40" << std::endl;
                std::cout
                        << " --parameter OR -P <INTEGER_VALUE>\tIndicates parameter, 0 for \n\t\t\t\tsearch for parameter by increasement"
                        << std::endl << "\t\t\t\tDefault: 0 (search for parameter by \n\t\t\t\tincreasement)"
                        << std::endl;
                std::cout
                        << " --parameter_steps OR -S <INTEGER_VALUE>\tIndicates step size for \n\t\t\t\tparameter increasement \n\t\t\t\t 0 for automatic calculation of step size"
                        << std::endl << "\t\t\t\tDefault: 0" << std::endl;
                std::cout << " --unweighted OR -U \tUse on unweighted graph" << std::endl
                          << "\t\t\t\tDefault: weighted graph" << std::endl;
                std::cout
                        << " --max_component OR -M <INTEGER_VALUE> \tIndicates maximum size \n\t\t\t\t of components considered in the calculation \n\t\t\t\tDefault: 32767"
                        << std::endl;
                std::cout
                        << " --costparser OR -C <COST_PARSER>\tIndicates which cost parser will be \n\t\t\t\tused for transforming edge weights into costs"
                        << std::endl
                        << "\t\t\t\tOptions: BlastParser.exponentToCosts, \n\t\t\t\t BlastParser.simpleCutOff"
                        << std::endl << "\t\t\t\tDefault: BlastParser.exponentToCosts" << std::endl;
                std::cout
                        << " --graphparser OR -G <GRAPH_PARSER>\tIndicates which graph parser will be \n\t\t\t\tused to read an edge file"
                        << std::endl
                        << "\t\t\t\tOptions: EdgeFileParser.initFrom5ColumnFile,\n\t\t\t\t EdgeFileParser.initFrom3ColumnFile, \n\t\t\t\t EdgeFileParser.initFrom12ColumnFile"
                        << std::endl << "\t\t\t\tDefault: EdgeFileParser.initFrom5ColumnFile" << std::endl;
                std::cout
                        << " --all_solutions OR -A\tReturning all possible solutions\n\t\t\t\t(Note: not compatible with option -s) \n\t\t\t\tDefault: returning only the best solution"
                        << std::endl;
                std::cout << "DISABLE SPEED-UP OPTIONS" << std::endl;
                std::cout
                        << " --no_split OR -s\tDisable splitting of search tree\n\t\t\t\tfor connected components\n\t\t\t\tNote: Returns just optimal solutions\n\t\t\t\tDefault: Splitting enabled"
                        << std::endl;
                exit(0);
            } else {
                std::cout << "Unknown option! Use help to get more information!" << std::endl;
                exit(0);
            }
            g += 2;
        } else {
            if (files == 0) {
                fname = argv[g];
                g++;
                files++;
            } else if (files == 1) {
                output_fname = argv[g];
                g++;
                files++;
            } else if (files == 2) {
                log_fname = argv[g];
                g++;
                files++;
            } else {
                std::cout << "Wrong flag start symbol (should be '-')!" << std::endl;
                exit(0);
            }
        }
    }

}
