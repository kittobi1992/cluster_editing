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

#include <../tools/stopwatch.h>

using std::vector;
using std::rand;
using std::srand;
using std::cout;
using std::endl;
using std::ofstream;
using std::istringstream;
using std::string;
using std::swap;
using std::abs;



void solveForConnectedComponent(CostsGraph &G, short pi_type, bool merging, bool just_max, bool split, double &costs);
void createRandomGraph1(CostsGraph &CG, double mue_ex, double sigma_ex, double mue_in, double sigma_in, bool save_cr, bool unweigthed, int Graph_size, double &k, string cr_fname);
void createRandomGraph2(CostsGraph &CG, bool save_cr, bool unweigthed, string cr_fname, int Graph_size);	
double ranf();
double generateGaussianDistributedNumber(double mue, double sigma);
void createTesterGraph(CostsGraph &CG, int mue_ex, int sigma_ex, int mue_in, int sigma_in, bool save_cr, bool unweigthed, string cr_fname, double &k, int Graph_size, short  RandomGraph, int &k2);
string toString(int x);
void fillParameters(int argc, char **argv, int &mue_ex, int &sigma_ex, int &mue_in, int &sigma_in, int &Graph_size, bool &unweighted, bool &test, short &RandomGraph, int &k2, int &runs, string &output_path);



int main(int argc, char **argv) {
	// initialize parameters with standard settings
	string cm_fname = "";
	string cr_fname = "";
	string output_path = "";
	bool merging = true;
	bool just_max = true;
	bool split = true;
	int mue_ex = -21;
	int sigma_ex = 20;
	int mue_in = 21;
	int sigma_in = 20;
	bool unweigthed = false;
	bool test = true;
	short RandomGraph = 1;
	int Graph_size = 0;
	int k2 = -1;
	// Number of the negative existig edges (falsch-negativ) 
	int sum = 0;
	int runs = 1;
	// fill parameters with values from shell
	fillParameters(argc, argv, mue_ex, sigma_ex, mue_in, sigma_in, Graph_size, unweigthed, test, RandomGraph, k2, runs, output_path);
	double k = 0.0;
	srand(time(NULL));
	for (int i=0; i < runs; i++){	
		try {
			CostsGraph CG;
			// create graph file name
			cm_fname=output_path;
			cm_fname=cm_fname.append("graph_").append(toString(i)).append(".cm");
			// create result file name
			cr_fname=output_path;
			cr_fname=cr_fname.append("graph_").append(toString(i)).append(".cr");
			cout << "-----------------------------------" << endl;
			createTesterGraph(CG, mue_ex, sigma_ex, mue_in, sigma_in, !test, unweigthed, cr_fname, k, Graph_size, RandomGraph, k2);
			if (!test) {
				std::cout << "try to write " << cm_fname << std::endl;
				CG.costsIOOperation(cm_fname.c_str(),MatrixParser::toMatrixFile); 
				std::cout << "wrote.." << std::endl;
			} else {
				std::cout << "start solving..." << std::endl;
				double costs = 0;
				CG.costsIOOperation("error.cm",MatrixParser::toMatrixFile);
				solveForConnectedComponent(CG, 2, merging, just_max, split, costs);
				if (costs-0.00000001 > k){
					cout << "k: "<< costs << endl;
					cout << "ERROR!!!" << endl;
	        			cout << "costs:  " << costs << " > "<<k<<endl;
					CG.costsIOOperation("error.cm",MatrixParser::toMatrixFile); 
					exit(1);
				} else {
					cout << "k: "<< costs <<endl;
					cout << "Done!!!" << endl;
				}
				cout << "-----------------------------------" << endl;
			}
		} catch ( GraphException e ) {
			cout << "Graph Exception: " << e.getMessage() << endl;
		}
	}
	return 0;
}


void solveForConnectedComponent(CostsGraph &G, short pi_type, bool merging, bool just_max, bool split, double &costs){
	int done = 0;
	//double sum = 0;
	double par_steps = 1;
	double par = 0;
	/*if (pi_type == 2) {
		for (int i=0; i < G.getSize(); i++) {
			for(int k=0; k < i;k++) {
		        	double e = G.getEdge(i,k);
				if (static_cast<int>(e/G.forbidden+0.001) != 1  && static_cast<int>(e/G.permanent+0.001) != 1) {
					sum += abs(e);
				}
			}
		}
		sum = sum / (G.getSize()*G.getSize());
		par_steps = (((G.getSize() / 5 * sum) + 1));
	}
	par_steps = par_steps*3;*/
	
	if (G.getSize() <= 2) {
		done = 1;
	}

/*	// first do cc reduction
	WeightedProblemInstance WPI = WeightedProblemInstance(G, 1E12, true);*/
	// start timer
	double used_time = 0;
	Stopwatch timer;
/*	timer.start();
	std::cout << "start cc reduction..." << std::endl;
	int count = WPI.ccKernelization();
	used_time = timer.elapsed();
	std::cout << "cc reduction in " << used_time << " sec." << std::endl;
	std::cout << "  reduced from " << G.getSize() << " to " << (WPI.getGraph()).getSize() << std::endl;
	std::cout << "  modification costs are " << (1E12-WPI.getParameter())<< std::endl;
	//costs = (1E12-WPI.getParameter());
	//G = WPI.getGraph();*/

	// do pid reduction with upper bound
	std::cout << "start pid reduction with upper bound..." << std::endl;
	used_time = 0;
	timer.start();
	WeightedProblemInstance PI_B = WeightedProblemInstance(G, 1E12);
	double upper_bound = PI_B.getUpperBound();
	WeightedProblemInstance PI_R = WeightedProblemInstance(G, upper_bound);
	int count = PI_R.maxReduce();
	used_time = timer.elapsed();
	double red_costs = upper_bound - PI_R.getParameter();
	std::cout << "pid reduction with upper bound in " << used_time << " sec." << std::endl;
	std::cout << "  reduced from " << G.getSize() << " to " << (PI_R.getGraph()).getSize() << std::endl;
	std::cout << "  modification costs are " << red_costs << std::endl;
	costs = red_costs;
	G = PI_R.getGraph();

	/*// do pid reduction
	WeightedProblemInstance WPI2 = WeightedProblemInstance(G, 1E12, true);
	used_time = 0;
	timer.start();
	std::cout << "start pid reduction..." << std::endl;
	count = WPI2.maxReduce();
	used_time = timer.elapsed();
	std::cout << "pid reduction in " << used_time << " sec." << std::endl;
	std::cout << "  reduced from " << G.getSize() << " to " << (WPI2.getGraph()).getSize() << std::endl;
	std::cout << "  modification costs are " << (1E12-WPI2.getParameter())<< std::endl;
	costs = (1E12-WPI2.getParameter());
	G = WPI2.getGraph();*/

	WeightedProblemInstance* PI = NULL;
	
	/*get lower bound*/
	PI = new WeightedProblemInstance(G, 1E12, false);
	par = PI->getLowerBound();
	upper_bound = PI->getUpperBound();
	delete(PI);
	PI=NULL;
	
	par_steps = (upper_bound - par) / 2;

	// do real programm
	used_time = 0;
	timer.start();
	std::cout << "start fpt algorithm..." << std::endl;
	while( done == 0 ) {
		try {
			par += par_steps;
			std::cout << "try with " << par-par_steps << std::endl;
			PI = new WeightedProblemInstance(G, par - par_steps);
			
			count = PI->maxReduce();
			// SEARCH TREE !!!
			std::cout << "create search tree...size here=" << (PI->getGraph()).getSize() << std::endl;
			SearchTreeWeighted ST(*PI, 0);
			ST.search(just_max, split);
			SearchTreeWeighted::solutions_type solutions = ST.getSolutions();
			if (solutions.size() > 0) {
				double min = 100000E100;
				for (int i = 0; i < solutions.size(); i++) {
					double costs = solutions[i].start_parameter - solutions[i].parameter;
					min = (costs < min) ? costs : min; 
				}
				used_time = timer.elapsed();
				std::cout << "fpt algorithm in " << used_time << " sec." << std::endl;
				std::cout << "  modification costs are " << min<< std::endl;

				costs += min;
				done = 1;
			}
			delete(PI);
			PI = NULL;
		} 
		catch ( ProblemInstanceException e ) {
			//cout << "Problem Instance Exception: " << e.getMessage() << endl;
			delete(PI);
			PI = NULL;
		} 
		catch ( GraphException e ) {
			cout << "Graph Exception: " << e.getMessage() << endl;
			delete(PI);
			PI = NULL;
		}
	}
}


void createRandomGraph1 (CostsGraph &CG, double mue_ex, double sigma_ex, double mue_in, double sigma_in, bool save_cr, bool unweigthed, int Graph_size, double &k, string cr_fname){
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
	if (Graph_size==0){
		node = rand() % 100 + 10;
	}
	else{
		node=Graph_size;
	}
        cout << "Anzahl der Konten: "<<node<<endl;
	CostsGraph::double_matrix_type Graph;
	//vector <vector <double> > G = vector <vector <double> > (node, vector < double > (node,0.0));
	vector <vector<int> > ClusterList;
	/* initialize random seed: */
	help=node;
	// find the clusters
	for (i; i < node; i++){
		size = rand() % help+1;
		ClusterList.push_back(vector <int> (size,0));
		for (int j=0; j < size; j++){
			ClusterList[i][j]=help-j-1;
		}
		help-=size;
		//cout << "Cluster" << i << " : " << size << endl;
		if (help==0){
			break;
		}
		// if less then 5 nodes are left, these will be the last cluster
		if (help < 10){
			i++;
			size = help;
			ClusterList.push_back(vector <int> (size,0));
			for (int j=0; j < size;j++){
				ClusterList[i][j]=help-j-1;
			}
			help-=size;
			//cout << "Cluster" << i << " : " << size << endl;
			break;
		} 
	}
	cluster=i+1;
        cout << "Anzahl der Cluster: "<< cluster << endl;
	/*cout << "Ausgabe der Cluster: " << endl;
	for (int j=0; j < ClusterList.size(); j++){
		cout << "Cluster"<<j<<": ";
		for (int k=0; k < ClusterList[j].size(); k++){
			cout<<ClusterList[j][k]<< " ";
		}
		cout <<endl;
	}*/
	if(save_cr){
		cr_file.open (cr_fname.c_str());
		if (!cr_file) {
			cout << "Can not open output file!" << endl;
			exit(0);
		}
	}
	int l=0;
	// fill the graph with vaules
	for (int j = 0; j < node; j++) {
		Graph.push_back(CostsGraph::double_matrix_row_type());
		for (int m = 0; m < j; m++) {
			//cout << "j: " << j << ", m: " << m << endl; 
			bool help=true;
			l=0;
			while (help){
				//cout << "Cluster"<<l<<": ";
				index1 = ClusterList[l][0];
				index2 = ClusterList[l][ClusterList[l].size()-1];
				//cout << "Index1: " << index1 << ", Index2: " << index2 << endl;
				if ((j >= index2) && (j <= index1)){
					if (m == 0){
						cr_file << l << " ";
					}
					if ((m >= index2)){
						h = generateGaussianDistributedNumber(mue_in, sigma_in);
						if(unweigthed){
							if (h < 0){
								h=-1.0;
								k+=1;
					
							}
							else{
								h=1.0;
							}
						}
						else if (h < 0){
							k+=-h;
						}
					}
					else{
						h = generateGaussianDistributedNumber(mue_ex, sigma_ex);
						if(unweigthed){
							h=0.0;
							if (h > 0){
								h=1.0;
								k+=1;
					
							}
							else{
								h=-1.0;
							}
						}
						else if (h > 0){
							k+=h;
						}
					}
				help=false;
				}
			l++;
			}
			Graph[j].push_back(h);
			//G[j][k]=h;
			//G[k][j]=h;
		}
	}
	CostsGraph::vertex_name_type EdgeNameList;
	for (int i=0; i < Graph.size(); i++){
		EdgeNameList.push_back(toString(i));
	}
	CG = CostsGraph(Graph.size(), Graph, EdgeNameList);
	
	if (save_cr){
		cr_file.close();
	}
}


void createRandomGraph2 (CostsGraph &CG, bool save_cr, bool unweigthed, string cr_fname, int Graph_size){
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
	if (Graph_size==0){
		node = rand() % 100 + 10;
	}
	else{
		node=Graph_size;
	}
        cout << "Anzahl der Konten: "<<node<<endl;
	CostsGraph::double_matrix_type Graph;
	//vector <vector <double> > G = vector <vector <double> > (node, vector < double > (node,0.0));
	vector <vector <int> > ClusterList;
	/* initialize random seed: */
	help=node;
	// find the clusters
	for (i; i < node; i++){
		size = rand() % help+1;
		ClusterList.push_back(vector <int> (size,0));
		for (int j=0; j < size; j++){
			ClusterList[i][j]=help-j-1;
		}
		help-=size;
		//cout << "Cluster" << i << " : " << size << endl;
		if (help==0){
			break;
		}
		// if less then 5 nodes are left, these will be the last cluster
		if (help < 5){
			i++;
			size = help;
			ClusterList.push_back(vector <int> (size,0));
			for (int j=0; j < size;j++){
				ClusterList[i][j]=help-j-1;
			}
			help-=size;
			//cout << "Cluster" << i << " : " << size << endl;
			break;
		} 
	}
	cluster=i+1;
        cout << "Anzahl der Cluster: "<< cluster << endl;
	/*cout << "Ausgabe der Cluster: " << endl;
	for (int j=0; j < ClusterList.size(); j++){
		cout << "Cluster"<<j<<": ";
		for (int k=0; k < ClusterList[j].size(); k++){
			cout<<ClusterList[j][k]<< " ";
		}
		cout <<endl;
	}*/
	if(save_cr){
		cr_file.open (cr_fname.c_str());
		if (!cr_file) {
			cout << "try to open: " << cr_fname << std::endl;
			cout << "Can not open output file!" << endl;
			exit(0);
		}
	}
	int l=0;
	// fill the graph with vaules
	for (int j = 0; j < node; j++) {
		Graph.push_back(CostsGraph::double_matrix_row_type());
		for (int k = 0; k < j; k++) {
			//cout << "j: " << j << ", k: " << k << endl; 
			bool help=true;
			l=0;
			while (help){
				//cout << "Cluster"<<l<<": ";
				index1 = ClusterList[l][0];
				index2 = ClusterList[l][ClusterList[l].size()-1];
				//cout << "Index1: " << index1 << ", Index2: " << index2 << endl;
				if ((j >= index2) && (j <= index1)){
					if (k == 0){
						cr_file << l << " ";
					}
					if ((k >= index2)){
						if(unweigthed){
							h=1.0;
						}
						else{
							h = ((double) rand() / (double) RAND_MAX) * 100.0;
							if(h < 0){
								h=-h;
							}
						}
					}
					else{
						if(unweigthed){
							h=-1.0;
						}
						else{
							h = ((double) rand() / (double) RAND_MAX) * 100.0;
							if(h > 0){
								h=-h;
							}
						}
					}
				help=false;
				}
			l++;
			}
			Graph[j].push_back(h);
			//G[j][k]=h;
			//G[k][j]=h;
		}
	}
	CostsGraph::vertex_name_type EdgeNameList;
	for (int i=0; i < Graph.size(); i++){
		EdgeNameList.push_back(toString(i));
	}
	CG = CostsGraph(Graph.size(), Graph, EdgeNameList);
	
	if (save_cr){
		cr_file.close();
	}
	//std::cout << CG << std::endl;
}

/*Random Numbers (uniform distributed)*/
double ranf(){
	//srand ( time(NULL) );
	return (double) rand() / (double) RAND_MAX;
}


/*
Box-Müller Transformation
http://en.wikipedia.org/wiki/Box_Muller
http://www.taygeta.com/random/gaussian.html
*/
double generateGaussianDistributedNumber(double m, double s){
	double x1, x2, w, y1;
	static double y2;
	static int use_last = 0;

	if (use_last){
		y1 = y2;
		use_last = 0;
	}
	else{
		do {
			x1 = 2.0 * ranf() - 1.0;
			x2 = 2.0 * ranf() - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );
		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}
	return( m + y1 * s );
}


void createTesterGraph(CostsGraph &CG, int mue_ex, int sigma_ex, int mue_in, int sigma_in, bool save_cr, bool unweigthed, string cr_fname, double &k, int graph_size, short  RandomGraph, int &k2){
	// k2 number of the edges to change
	k = 0.0; // costs of the changed edges
	int size;
	int node1;
	int node2;
	double help;
	int i = 0;
	if (RandomGraph==1){
		createRandomGraph1 (CG, mue_ex, sigma_ex, mue_in, sigma_in, save_cr, unweigthed, graph_size, k, cr_fname);
	}
	else{
		if (RandomGraph==2){
			createRandomGraph2 (CG, save_cr, unweigthed, cr_fname, graph_size);
		}
		else{
			cout << "Error RandomGraph has to be 1 or 2" << endl;
			exit(0);
		}
	
		size = CG.getSize();
		if (k2 < 0){
			k2 = size/2;
		}
		
		vector <vector <bool> > changed_edges = vector <vector <bool> >(size, vector<bool>(size, false));
		while (i<k2){
			node1 = rand() % (size-1);
			node2 = rand() % (size-1);
			if (node1 == node2){
				continue;
			}
			if (node1 < node2){
				swap(node1,node2);
			}
			if (changed_edges[node1][node2]) {
				continue;
			}
			help = CG.getEdge(node1,node2);
			CG.setEdge(node1,node2,-help);
			k+=abs(help);
			i++;
		}

	}
}


string toString(int x){	
	std::ostringstream o;
	o << x;
	return o.str();
}

void fillParameters(int argc, char **argv, int &mue_ex, int &sigma_ex, int &mue_in, int &sigma_in, int &Graph_size, bool &unweighted, bool &test, short &RandomGraph, int &k2, int &runs, string &output_path){
	int g = 1;
	int files = 0;
	while( argc > g ) {
		if (argv[g][0] == '-') {
			string flag(argv[g]);
			flag = flag.erase(0,1);

			if (flag == "-output_path" || flag == "O") {
				string argument(argv[g+1]);
				output_path=argument;
			}
			else if (flag == "-mu_ex" ) {
				istringstream i(argv[g+1]);
				i >> mue_ex;
			}
			else if (flag == "-sigma_ex" ) {
				istringstream i(argv[g+1]);
				i >> sigma_ex;
			}
			else if (flag == "-mu_in" ) {
				istringstream i(argv[g+1]);
				i >> mue_in; 
			}
			else if (flag == "-sigma_in" ) {
				istringstream i(argv[g+1]);
				i >> sigma_in;
			}
			else if (flag == "-size" ) {
				istringstream i(argv[g+1]);
				i >> Graph_size;
			}
			else if (flag == "-unweighted" || flag == "U") {
				unweighted = true;
				g--;
			}
			else if (flag == "-sample" || flag == "S") {
				test = false;
				g--;
			}
			else if (flag == "-RandomGraph" ) {
				istringstream i(argv[g+1]);
				i >> RandomGraph;

			}
			else if (flag == "-k2" ) {
				istringstream i(argv[g+1]);
				i >> k2;
			}
			else if (flag == "-runs" ) {
				istringstream i(argv[g+1]);
				i >> runs;

			}
			else if (flag == "-help") {
				cout << "-----------------------------------" << endl;
				cout << "| Tester Tool for Cluster Editing |" << endl;
				cout << "-----------------------------------" << endl;
				cout << "You can create a random CostGraph and write a *.cm and *.cr file, or test the Cluster Editing Tool with a random Graph." << endl;
				cout << "--------------------" << endl;
				cout << "Options:" << endl << endl;
				cout << " --RandomGraph [1|2] (Default: 1)" << endl;
				cout << "\t 1 = similarities between mates and beween nomates are Gaussian distributed." << endl;
				cout << "\t\t --mu_ex (mue-value for non-existing edges. Default: -21)" << endl;
				cout << "\t\t --sigma_ex (sigma-value for non-existing edges. Default: 20)" << endl;
				cout << "\t\t --mu_in (mue-value for edges. Default: 21)" << endl;
				cout << "\t\t --sigma_in (sigma-value for edges. Default: 20)" << endl;
				cout << "\t 2 = random graph with randomly (uniformlly) switched edges (-s(uv))" << endl;
				cout << "\t\t --k2 (number of edge modifications). Default: |V|/2" << endl;
				cout << " --sample or -S if you create random sample graphs. Default: false" << endl;
				cout << " --output_path or -O needs to be indicated if --test option is not used" << endl;
				cout << " --size if you want to set the number of nodes in the Graph to a specific value. Default: random number of nodes, 10 to 100." << endl;
				cout << " --runs to specify the number of runs. You shoud not specify a filename, if you want to create more than on Random CostsGraph. \n Must be unsigned int, else default. Default: 1" << endl;
				cout << " --unweighted OR -U If you want to work with an unweighted graph. Default: weighted graph" << endl;
				
				exit(0);
			} else {
				cout << "Unknown option! Use help to get more information!" << endl;
				exit(0);
			}
			g += 2;
		}
		 else {
		 	cout << "Wrong flag start symbol (should be '-')!" << endl;
		 	exit(0);
		}
	}
}

