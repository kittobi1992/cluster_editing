#include <fstream>
#include <iostream>
#include <sstream>

#include <sys/time.h>
#include <math.h>
#include <vector>

#include <searchtreeweighted.h>
#include <weightedprobleminstance.h>
#include <graphexception.h>
#include <costsgraph.h>
#include <matrixparser.h>
#include <stopwatch.h>

using std::vector;
using std::rand;
using std::srand;
using std::cout;
using std::endl;
using std::ofstream;
using std::istringstream;
using std::string;


void setPath(int Vred, int krest, string &path, bool &help);

void getBounds(CostsGraph &G, double &lb, double &ub);

int solve(CostsGraph &G, short pi_type, bool merging, bool just_max, bool split, double &lb, double &ub);

void createRandomGraph(CostsGraph &CG, int Graph_size);

void createTesterGraph(CostsGraph &CG, int graph_size, int c);

string toString(int x);

void
fillParameters(int argc, char **argv, int &Graph_size, double &c, double &c_max, int &nr, int &number, int &size_max,
               int &runs);


int main(int argc, char **argv) {
    // initialize parameters with standard settings
    string cm_fname = "";
    string path = "/vol/clusterdata/cluster_editing/unweighted_data/reduced_samples_V_50to1000_k_0.5to10V/";
    string path_help;
    bool merging = true;
    bool just_max = true;
    bool split = true;
    int Graph_size = 0;
    double c = -1;
    double c_max = 0;
    int nr = 0;
    int number = 10;
    int size_max = 0;
    int Vred;
    int Vorg;
    int krest;
    int korg;
    int mod; // number of edge modifications during reduction
    // Number of the negative existig edges (falsch-negativ)
    int sum = 0;
    int runs = 1;
    bool help = true;
    double ub; //upper_bound
    double lb; //lower_bound

    // fill parameters with values from shell
    fillParameters(argc, argv, Graph_size, c, c_max, nr, number, size_max, runs);

    Graph_size = (Graph_size / 10) * 10;
    size_max = (size_max / 10) * 10;

    srand(time(NULL));
    int help2 = 0;
    if ((c_max > 0) || (size_max > 0) || (number > 10)) {
        for (int i = Graph_size; i < size_max + 1; i += 10) {
            for (int j = c * i; j <= c_max * i; j += (i / 2)) {
                std::cout << "number=" << number << std::endl;
                for (int k = 1; k < number + 1; k++) {
                    try {
                        CostsGraph CG;
                        createTesterGraph(CG, i /*Graph_size*/, (j));

                        Vorg = CG.getSize();
                        korg = j - (pow(j, 2) / ((pow(Vorg, 2) - Vorg) / 2) * 0, 8);

                        // use bounds of non-reduced instance as correction for korg estimation
                        getBounds(CG, lb, ub);
                        if (korg < lb || korg > ub) {
                            // use this correction if estimation is corrupt
                            korg = (lb + ub) / 2;
                        }

                        cout << "-------------------------------------------------------------------" << endl;
                        // Reduzieren
                        // ub_after and lb_after are saved only temporary to correct the residual costs
                        double lb_after;
                        double ub_after;
                        mod = solve(CG, 2, merging, just_max, split, lb_after, ub_after);
                        GraphSet graph_set = GraphSet(CG);
                        for (int l = 0; l < graph_set.getSetSize(); l++) {
                            CostsGraph &G = graph_set.getGraph(l);
                            //CostsGraph G = graph_set.getGraph(l);

                            Vred = G.getSize();
                            cout << "Vred: " << Vred << ", size CG: " << CG.getSize() << endl;
                            krest = korg - mod;
                            if (krest < 0) {
                                krest = 0;
                            }
                            if ((krest < lb_after) || (krest > ub_after)) {
                                krest = (lb_after + ub_after) / 2;
                            }
                            // create graph file name
                            // to be clear: lb and ub are the bounds before reduction starts
                            // lb_after and ub_after are used to correct the residual costs estimation
                            // in the same manner we use lb and ub as correction for the inital costs before reduction starts
                            string cm_fname = "graph_";
                            cm_fname = cm_fname.append("Vred_").append(toString(Vred)).append("_krest_").append(
                                    toString(krest)).append("_c_").append(toString(j)).append("_Vorg_").append(
                                    toString(Vorg)).append("_korg_").append(toString(korg)).append("_lb_").append(
                                    toString((int) lb)).append("_ub_").append(toString((int) ub)).append("_nr_").append(
                                    toString(nr)).append(".cm");
                            if (Vred < 3) {
                                help = false;
                            } else {
                                path_help = path;
                                setPath(Vred, krest, path_help, help);
                            }
                            if (help) {
                                path_help = path_help.append(cm_fname);
                                cout << path_help << endl;
                                G.costsIOOperation(path_help.c_str(), MatrixParser::toMatrixFile);
                                help2++;
                                cout << "Vred: " << Vred << ", krest: " << krest << ", c: " << (j) << ", Vorg: " << Vorg
                                     << ", korg: " << korg << ", nr: " << nr << ", #Komponenten: " << help2 << endl;
                                nr++;
                                cout << "-------------------------------------------------------------------" << endl;
                                //continue;
                            }
                            help = true;
                        }
                        help2 = 0;

                    } catch (GraphException e) {
                        cout << "Graph Exception: " << e.getMessage() << endl;
                    }
                    nr++;
                }
            }
            if ((i > 100) && (i % 50 != 0)) {
                i = i + 40 - (i % 50);
            } else if (i >= 100) {
                i += 40;
            }
        }
    } else {
        int help_int = 0;
        for (int i = nr; i < runs + nr; i++) {

            try {
                CostsGraph CG;

                createTesterGraph(CG, Graph_size, c * Graph_size);

                Vorg = CG.getSize();
                korg = c * Graph_size - (pow(c, 2) / ((pow(Vorg, 2) - Vorg) / 2) * 0, 8);

                // use bounds of non-reduced instance as correction for korg estimation
                getBounds(CG, lb, ub);
                if (korg < lb || korg > ub) {
                    // use this correction if estimation is corrupt
                    korg = (lb + ub) / 2;
                }

                cout << "-------------------------------------------------------------------" << endl;
                // Reduzieren
                double lb_after;
                double ub_after;
                mod = solve(CG, 2, merging, just_max, split, lb_after, ub_after);

                GraphSet graph_set = GraphSet(CG);

                for (int l = 0; l < graph_set.getSetSize(); l++) {
                    //CostsGraph& G = graph_set.getGraph(l);
                    CostsGraph G = graph_set.getGraph(l);

                    Vred = G.getSize();
                    cout << "Vred: " << Vred << ", size CG: " << CG.getSize() << endl;
                    krest = korg - mod;
                    if (krest < 0) {
                        krest = 0;
                    }
                    if ((krest < lb_after) || (krest > ub_after)) {
                        krest = (lb_after + ub_after) / 2;
                    }

                    // create graph file name
                    string cm_fname = "graph_";
                    cm_fname = cm_fname.append("Vred_").append(toString(Vred)).append("_krest_").append(
                            toString(krest)).append("_c_").append(toString(c * Graph_size)).append("_Vorg_").append(
                            toString(Vorg)).append("_korg_").append(toString(korg)).append("_lb_").append(
                            toString((int) lb)).append("_ub_").append(toString((int) ub)).append("_nr_").append(
                            toString(nr)).append(".cm");
                    if (Vred < 3) {
                        help = false;
                    } else {
                        path_help = path;
                        setPath(Vred, krest, path_help, help);
                    }
                    if (help) {
                        path_help = path_help.append(cm_fname);
                        cout << path_help << endl;
                        G.costsIOOperation(path_help.c_str(), MatrixParser::toMatrixFile);
                        help2++;
                        cout << "B: Vred: " << Vred << ", krest: " << krest << ", c: " << (c * Graph_size) << ", Vorg: "
                             << Vorg << ", korg: " << korg << ", nr: " << help_int << ", #Komponenten: " << help2
                             << endl;
                        cout << "-------------------------------------------------------------------" << endl;
                        help_int++;
                        //continue;
                    }
                    help_int++;
                    help = true;
                }
                help2 = 0;
            } catch (GraphException e) {
                cout << "Graph Exception: " << e.getMessage() << endl;
            }
            help_int++;
        }
    }
    cout << "Job done!" << endl;
    return 0;
}


void setPath(int Vred, int krest, string &path_help, bool &help) {
    path_help = path_help.append("random_graphs_reduced_");
    switch (Vred / 50) {
        case 0:
            path_help = path_help.append("V_1-49_");
            break;
        case 1:
            path_help == path_help.append("V_50-99_");
            break;
        case 2:
            path_help = path_help.append("V_100-149_");
            break;
        case 3:
            path_help == path_help.append("V_150-199_");
            break;
        case 4:
            path_help = path_help.append("V_200-249_");
            break;
        case 5:
            path_help == path_help.append("V_250-299_");
            break;
        case 6:
            path_help = path_help.append("V_300-349_");
            break;
        case 7:
            path_help == path_help.append("V_350-399_");
            break;
        case 8:
            path_help = path_help.append("V_400-449_");
            break;
        case 9:
            path_help == path_help.append("V_450-499_");
            break;
        case 10:
            path_help == path_help.append("V_500-549_");
            break;
        case 11:
            path_help == path_help.append("V_550-599_");
            break;
        case 12:
            path_help == path_help.append("V_600-649_");
            break;
        case 13:
            path_help == path_help.append("V_650-699_");
            break;
        case 14:
            path_help == path_help.append("V_700-749_");
            break;
        case 15:
            path_help == path_help.append("V_750-799_");
            break;
        case 16:
            path_help == path_help.append("V_800-849_");
            break;
        case 17:
            path_help == path_help.append("V_850-899_");
            break;
        case 18:
            path_help == path_help.append("V_900-949_");
            break;
        case 19:
            path_help == path_help.append("V_950-999_");
            break;
            //default: help=false;
        default:
            path_help == path_help.append("V_1000_");
    }

    switch ((int) (4.0 * krest / Vred)) {
        case 0:
        case 1:
        case 2:
            path_help = path_help.append("k_0.5V/");
            break;
        case 3:
        case 4:
            path_help == path_help.append("k_1V/");
            break;
        case 5:
        case 6:
            path_help = path_help.append("k_1.5V/");
            break;
        case 7:
        case 8:
            path_help == path_help.append("k_2V/");
            break;
        case 9:
        case 10:
            path_help = path_help.append("k_2.5V/");
            break;
        case 11:
        case 12:
            path_help == path_help.append("k_3V/");
            break;
        case 13:
        case 14:
            path_help = path_help.append("k_3.5V/");
            break;
        case 15:
        case 16:
            path_help == path_help.append("k_4V/");
            break;
        case 17:
        case 18:
            path_help = path_help.append("k_4.5V/");
            break;
        case 19:
        case 20:
            path_help == path_help.append("k_5V/");
            break;
        case 21:
        case 22:
            path_help == path_help.append("k_5.5V/");
            break;
        case 23:
        case 24:
            path_help == path_help.append("k_6V/");
            break;
        case 25:
        case 26:
            path_help == path_help.append("k_6.5V/");
            break;
        case 27:
        case 28:
            path_help == path_help.append("k_7V/");
            break;
        case 29:
        case 30:
            path_help == path_help.append("k_7.5V/");
            break;
        case 31:
        case 32:
            path_help == path_help.append("k_8V/");
            break;
        case 33:
        case 34:
            path_help == path_help.append("k_8.5V/");
            break;
        case 35:
        case 36:
            path_help == path_help.append("k_9V/");
            break;
        case 37:
        case 38:
            path_help == path_help.append("k_9.5V/");
            break;
        case 39:
        case 40:
            path_help == path_help.append("k_10V/");
            break;
            //default: help=false;
        default:
            path_help == path_help.append("k_10V/");
    }
}

void getBounds(CostsGraph &G, double &lb, double &ub) {
    WeightedProblemInstance PI_B = WeightedProblemInstance(G, 1E12);
    ub = PI_B.getUpperBound();
    lb = PI_B.getLowerBound();
}

int solve(CostsGraph &G, short pi_type, bool merging, bool just_max, bool split, double &lb, double &ub) {
    int done = 0;
    double par_steps = 1;
    double par = 0;

    if (G.getSize() <= 2) {
        done = 1;
    }
    // start timer
    double used_time = 0;
    Stopwatch timer;

    // do pid reduction with upper bound
    std::cout << "start pid reduction with upper bound..." << std::endl;
    used_time = 0;
    timer.start();
    WeightedProblemInstance PI_B = WeightedProblemInstance(G, 1E12);
    double upper_bound = PI_B.getUpperBound();
    WeightedProblemInstance PI_R = WeightedProblemInstance(G, upper_bound);
    int count = PI_R.maxReduce();
    used_time = timer.elapsed();
    ub = PI_R.getUpperBound();
    lb = PI_R.getLowerBound();
    double red_costs = upper_bound - PI_R.getParameter();
    std::cout << "pid reduction with upper bound in " << used_time << " sec." << std::endl;
    std::cout << "  reduced from " << G.getSize() << " to " << (PI_R.getGraph()).getSize() << std::endl;
    std::cout << "  modification costs are " << red_costs << std::endl;
    G = PI_R.getGraph();

    // das hier brauchst du nicht mehr...macht ja nix sinnvolles??!?!!!
    /*WeightedProblemInstance* PI = NULL;

    PI = new WeightedProblemInstance(G, 1E12, false);
    par = PI->getLowerBound();
    upper_bound = PI->getUpperBound();
    delete(PI);
    PI=NULL;*/
    return (red_costs);
}

void createRandomGraph(CostsGraph &CG, int Graph_size) {
    ofstream cr_file;
    //Number of nodes
    int node;
    //Number of clusters
    int cluster;
    int i = 0;
    int help;
    double h;
    int index1;
    int index2;
    //size of the clusters
    int size;
    bool test;
    double x;
    double y;
    /* generate a random number for # nodes: */
    if (Graph_size == 0) {
        node = rand() % 100 + 10;
    } else {
        node = Graph_size;
    }
    CostsGraph::double_matrix_type Graph;
    vector<vector<int> > ClusterList;
    /* initialize random seed: */
    help = node;
    // find the clusters
    for (i; i < node; i++) {
        size = rand() % help + 1;
        ClusterList.push_back(vector<int>(size, 0));
        for (int j = 0; j < size; j++) {
            ClusterList[i][j] = help - j - 1;
        }
        help -= size;
        if (help == 0) {
            break;
        }
        // if less then 5 nodes are left, these will be the last cluster
        if (help < 5) {
            i++;
            size = help;
            ClusterList.push_back(vector<int>(size, 0));
            for (int j = 0; j < size; j++) {
                ClusterList[i][j] = help - j - 1;
            }
            help -= size;
            break;
        }
    }
    cluster = i + 1;
    int l = 0;
    // fill the graph with vaules
    for (int j = 0; j < node; j++) {
        Graph.push_back(CostsGraph::double_matrix_row_type());
        for (int k = 0; k < j; k++) {
            bool help = true;
            l = 0;
            while (help) {
                index1 = ClusterList[l][0];
                index2 = ClusterList[l][ClusterList[l].size() - 1];
                if ((j >= index2) && (j <= index1)) {
                    if (k == 0) {
                        cr_file << l << " ";
                    }
                    if ((k >= index2)) {
                        h = 1.0;
                    } else {
                        h = -1.0;
                    }
                    help = false;
                }
                l++;
            }
            Graph[j].push_back(h);
        }
    }
    CostsGraph::vertex_name_type EdgeNameList;
    for (int i = 0; i < Graph.size(); i++) {
        EdgeNameList.push_back(toString(i));
    }
    CG = CostsGraph(Graph.size(), Graph, EdgeNameList);

}

void createTesterGraph(CostsGraph &CG, int graph_size, int c) {
    // c number of the edges to change
    int size;
    int node1;
    int node2;
    double help;
    int i = 0;

    createRandomGraph(CG, graph_size);

    size = CG.getSize();
    if (c < 0) {
        c = size / 2;
    }

    vector<vector<bool> > changed_edges = vector<vector<bool> >(size, vector<bool>(size, false));
    while (i < c) {
        node1 = rand() % (size - 1);
        node2 = rand() % (size - 1);
        if (node1 == node2) {
            continue;
        }
        if (node1 < node2) {
            std::swap(node1, node2);
        }
        if (changed_edges[node1][node2]) {
            continue;
        }
        help = CG.getEdge(node1, node2);
        CG.setEdge(node1, node2, -help);
        i++;
    }
}


string toString(int x) {
    std::ostringstream o;
    o << x;
    return o.str();
}

void
fillParameters(int argc, char **argv, int &Graph_size, double &c, double &c_max, int &nr, int &number, int &size_max,
               int &runs) {

    int g = 1;
    int files = 0;
    while (argc > g) {
        if (argv[g][0] == '-') {
            string flag(argv[g]);
            flag = flag.erase(0, 1);
            if (flag == "-size") {
                istringstream i(argv[g + 1]);
                i >> Graph_size;
            } else if (flag == "-c") {
                istringstream i(argv[g + 1]);
                i >> c;
            } else if (flag == "-c_max") {
                istringstream i(argv[g + 1]);
                i >> c_max;
            } else if (flag == "-nr") {
                istringstream i(argv[g + 1]);
                i >> nr;
            } else if (flag == "-number") {
                istringstream i(argv[g + 1]);
                i >> number;
            } else if (flag == "-size_max") {
                istringstream i(argv[g + 1]);
                i >> size_max;
            } else if (flag == "-runs") {
                istringstream i(argv[g + 1]);
                i >> runs;

            } else if (flag == "-help") {
                cout << "-----------------------------------" << endl;
                cout << "| Graphproduzer |" << endl;
                cout << "-----------------------------------" << endl;
                cout
                        << "You can create a random unweigthed CostGraph and write a *.cm file. Filename: graph_x1_kred_x2_c_x3_Vorg_x4_korg_x5_nr_x6.cm"
                        << endl;
                cout << "--------------------" << endl;
                cout << "Options:" << endl << endl;
                cout << " --c (number of edge modifications). C is multiplyed with |V|. Default: |V|/2" << endl;
                cout
                        << " --c_max (number of edge modifications). If you want to create a set of graphs, with different c. The step width is 1."
                        << endl;
                cout
                        << " --size (|V|) if you want to set the number of nodes in the Graph to a specific value. Default: random number of nodes, 10 to 100."
                        << endl;
                cout
                        << " --size_max (|V|) if you want to create a set of graphs. The step width is 10, how long size is unter 100, and 50 over 100."
                        << endl;
                cout << " --number spezifes the number of graphs with a spezific size and a spezific c. Default: 10"
                     << endl;
                cout
                        << " --nr to specify the beginning number of a graph. Do not spezifie runs in this case, it will be ignored."
                        << endl;
                cout << " --runs to specify the number of runs. Must be unsigned int, else default. Default: 1" << endl;

                exit(0);
            } else {
                cout << "Unknown option! Use help to get more information!" << endl;
                exit(0);
            }
            g += 2;
        } else {
            cout << "Wrong flag start symbol (should be '-')!" << endl;
            exit(0);
        }
    }
}

