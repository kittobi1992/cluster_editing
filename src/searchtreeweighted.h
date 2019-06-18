/*+++++++++++++++++ class SearchTreeWeighted +++++++++++++++++++ */
#ifndef __searchtreeweighted_h__
#define __searchtreeweighted_h__

#include <vector>

#include <weightedprobleminstance.h>
#include <graphset.h>

// - object for search solutions for in the problem instance object
// - includes problem instance object as root and will run a search tree on it
// - solutions of the search tree are saved with changes, costs and real components as costs graphs
// - is able to create sub search tree objects if instance can be divided in smaller instances

class SearchTreeWeighted {
public:
	/* ####  Type definitions #### */

	// record to save a solution
	struct component_solution_type {
		// adjacency matrices = costs graphs
		GraphSet components;
		// changes made during the search
		WeightedProblemInstance::edge_list_type changes;

		// = end parameter
		double parameter;
		// = starting point
		double start_parameter;
	};

	// vector for saving all solutions
	typedef std::vector< component_solution_type > solutions_type;


	/* ####  Constants  #### */
	
	// constants for different adjacency status
	const double permanent;
	const double forbidden;


	/* ####  Constructor  #### */

	SearchTreeWeighted(WeightedProblemInstance &rootInstance, int d);


	/* ####  Destructor  #### */
	~SearchTreeWeighted();


	/* ####  main function #### */

	// searching for all / just best solution in the search tree
	void search(bool just_max, bool split);
	
	// returns solutions
	solutions_type getSolutions();

private:

	/* ####  members  #### */

	// root of the search tree
	WeightedProblemInstance _root;

	// list of solutions
	solutions_type _solutions;

	// options
	bool _split;
	bool _just_max;

	int _depth;

	// keep lower bound
	double _lower_bound;

	/* ####  methods  #### */

	// simple recursive function, each corresponds to a node in the search tree
	void recSearch(WeightedProblemInstance &node, int depth);

	// test whether the search tree and the instance shouls split
	bool splitWorthIt(CostsGraph::vertex_matrix_type vertex_matrix);

	// creates sub search tree objects and will run their search to get solutions
	// for sub instances of the given problem instance
	void splitSearch(WeightedProblemInstance &new_root, CostsGraph::vertex_matrix_type vertex_matrix, int depth);
};

#endif
