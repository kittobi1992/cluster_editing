#include <fstream>
#include <iostream>
#include <sstream>

#include <weightedprobleminstance.h>
#include <edgefileparser.h>
#include <matrixparser.h>

#include <stopwatch.h>

#include <math.h>

#include <blastparser.h>
#include <doubleparser.h>
#include <graphset.h>

void fillParameters(int argc, char **argv, std::string &fname, char* &output_fname);

int main(int argc, char **argv)
{
	// initialize parameters with standard settings
	std::string fname = "";
	char* output_fname = "";
	GraphSet::matrix_file_fct_type matrix_file_fct = MatrixParser::initFromMatrixFile;
	CostsGraph::matrix_file_fct_type write_graph_fct = MatrixParser::toMatrixFile;
	
	// fill parameters with values from shell
	fillParameters(argc, argv, fname, output_fname);

	std::cout << std::endl << "reduce " << fname  << std::endl;
	try {
		     // create main costs graph by reading from file
		     std::cout << "read graph file..." << std::endl;

		     CostsGraph G = CostsGraph(fname, matrix_file_fct);
		     int old_size = G.getSize();

		     // get rid of main graph
		     //delete main_graph;

		     double min_parameter = 0;
		     long solutions_nr = 0;
	
		     // solve problem for
		     WeightedProblemInstance PI = WeightedProblemInstance(G, 1E12, true);
		     
		     // start timer
		     double used_time = 0;
		     Stopwatch timer;
		     timer.start();

		     // reduce
		     int count = PI.ccKernelization();

		     // stop timer
		     used_time = timer.elapsed();

		     std::cout << "  Reduction in " << used_time << " sec caused " << count << " edge changes. Graph reduced from " << old_size << " to " << (PI.getGraph()).getSize() << std::endl;
		     std::cout << "  Modification costs are " << (1E12-PI.getParameter())<< std::endl;

		     GraphSet graph_set = GraphSet(PI.getGraph());
		     int max_size = 0;
		     for (int i = 0; i < graph_set.getSetSize(); i++ ) {
		         CostsGraph& G2 = graph_set.getGraph(i);
		         std::ostringstream o;
		         o << (i+1);
		         std::string file = output_fname + o.str() + ".cm";
		         //std::cout << "write " << file << std::endl;
		         G2.costsIOOperation(file,write_graph_fct);
		         max_size = (max_size<G2.getSize()) ? G2.getSize() : max_size;
		     }

		     std::cout << "  Instance splitted in " << graph_set.getSetSize() << " components. " << "Greates component of size: " << max_size << std::endl;

	} catch ( VertexListsException e ) {
		std::cout << "VertexLists-Exception: " << e.getMessage() << std::endl;
	} catch ( ProblemInstanceException e ) {
	} catch ( GraphException e ) {
		std::cout << "Graph Exception: " << e.getMessage() << std::endl;
	}

	return 0;
}


void fillParameters(int argc, char **argv, std::string &fname, char* &output_fname)
{
	int g = 1;
	int files = 0;
	while( argc > g ) {
		if (argv[g][0] == '-') {
		     std::string flag(argv[g]);
		     flag = flag.erase(0,1);
		     if (flag == "-help") {
		          std::cout << "Critical Clique Parameter Independent Reducer " << std::endl;
		          std::cout << "--------------------" << std::endl;
		          std::cout << "Usage: <input_file> <output_file>" << std::endl << std::endl;
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
		     } else {
		          std::cout << "Wrong flag start symbol (should be '-')!" << std::endl;
		          exit(0);
		     }
		}
	}

}
