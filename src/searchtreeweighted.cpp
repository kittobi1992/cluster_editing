#include <searchtreeweighted.h>

/* ########## Constructor ########### */

// initate this object with a root of the search tree
SearchTreeWeighted::SearchTreeWeighted(WeightedProblemInstance &root, int d) : _solutions(0), _root(root), permanent(root.permanent), forbidden(root.forbidden), _split(true), _depth(d)
{
	// get an initial lower bound for the root instance
	_lower_bound = root.getLowerBound();
}

/* ########## Destructor ########### */

SearchTreeWeighted::~SearchTreeWeighted()
{
}

/* ########## Methods ########### */


// return whether it is useful to split the search tree
bool SearchTreeWeighted::splitWorthIt(CostsGraph::vertex_matrix_type vertex_matrix)
{
	
	// if there is  more than one connected componentn..
	if (vertex_matrix.size() > 1) {
		
		// count number of not trivial instances and singletons
		int not_simple = 0;
		int singletons = 0;
		for (int i = 0; i < vertex_matrix.size(); i++) {
		     if (vertex_matrix[i].size() > 5) {
		          not_simple++;
		     } else {
		          singletons++;
		     }
		}
		
		// decide whether we should split the search tree or not
		return (not_simple > 1 || singletons > 5) ? true : false;
	} else {
		return false;
	}
}


// public search methods which start from the root and stops if a solution is found
void SearchTreeWeighted::search(bool just_max, bool split)
{
	_just_max = just_max;
	_split = split;

	// get connected components of instance graphs
	CostsGraph::vertex_matrix_type vertex_matrix = (_root.getGraph()).getConnectedVertices();

	// check whether to split already or not
	if (_split && splitWorthIt(vertex_matrix)) {
		// start many search trees
		splitSearch(_root, vertex_matrix, _depth);
	} else {
		// start search in this search tree
		recSearch(_root, _depth);
	}
}


// returned obtained solution saved in the global member variable
SearchTreeWeighted::solutions_type SearchTreeWeighted::getSolutions()
{
	return _solutions;
}


// start split tree search
void SearchTreeWeighted::splitSearch(WeightedProblemInstance &new_root, CostsGraph::vertex_matrix_type vertex_matrix, int depth)
{
	// divide Instance
	WeightedProblemInstance::pi_list_type pi_list = new_root.divideInstance(vertex_matrix);

	// create solutionscontainer to solutions in subtrees
	solutions_type solutions_container = solutions_type(vertex_matrix.size());

	// save initial parameter
	double parameter = (pi_list[0])->getParameter();

	// save changes already made
	WeightedProblemInstance::edge_list_type sum_of_changes =  new_root.getChangedEdges();

	// search for solutions in subtrees
	for (int i = 0; i < vertex_matrix.size(); i++) {
	
		// set parameter for subtree to rest of given parameter
		try {
		     (pi_list[i])->setParameter(parameter);

		     // create subtree and search for solutions
		     SearchTreeWeighted subtree(*pi_list[i], depth+1);
		     subtree.search(_just_max, _split);
		     solutions_type solutions = subtree.getSolutions();

		     // if any subtree is not solvable the whole tree is not solvable-->exception
		     if (solutions.size() < 1) {
		          throw ProblemInstanceException("branch not solvable");
		     }

		     // find best solution and save it in solution container
		     // get rid of rest
		     solutions_container[i] = solutions[0];
		     for (int k = 1; k < solutions.size(); k++) {
		          if ((solutions_container[i]).parameter < (solutions[k]).parameter) {
		               solutions_container[i] = solutions[k];
		          }
		     }

		     // add changes to sum_of_changes
		     sum_of_changes.insert(sum_of_changes.end(), solutions_container[i].changes.begin(), solutions_container[i].changes.end());
	
		     // get rest of parameter
		     parameter = (solutions_container[i]).parameter;

		     delete(pi_list[i]);
		     pi_list[i] = NULL;

		} catch (ProblemInstanceException pe) {
		     // delete all other created problem instances
		     for (int k = i; k < vertex_matrix.size(); k++) {
		          delete(pi_list[k]);
		     }
		     throw ProblemInstanceException("branch not solvable");
		}

	}

	// merge solution of subtrees to one solution, saved in the first solution
	for (int i = 1; i < vertex_matrix.size(); i++) {
		for (int k = 0; k < solutions_container[i].components.getSetSize(); k++) {
			solutions_container[0].components.addGraph( solutions_container[i].components.getGraph(k) );
		}
	}

	// consequently the end parameter is the end parameter of the last subtree
	solutions_container[0].parameter = solutions_container[vertex_matrix.size() - 1].parameter;

	// save the vector of changed edges in every component in the first component
	solutions_container[0].changes = sum_of_changes;

	// if not just maximum solution is requested or no solution so far...just add solution
	if (!_just_max || _solutions.size() == 0) {
		_solutions.insert(_solutions.end(), solutions_container[0]);

	// otherwise test if solution is better then existing one, and delete the rejected solution
	} else if (_solutions[0].parameter < solutions_container[0].parameter) {
		_solutions[0] = solutions_container[0];
	} 

}


// normal recursive search for a problem instance
// this function will call itself recursivley for every branch
void SearchTreeWeighted::recSearch(WeightedProblemInstance &node, int depth)
{
	if (node.getParameter() >= 0) {

		// if no solution can be obtained abort here
		if (_just_max && _solutions.size() > 0 && node.getParameter() < _solutions[0].parameter) {
		     return;
		}

		// if parameter is lower then lower bound check for a new one and check if instance is still solvable
		if (node.getParameter() < _lower_bound) {
		     _lower_bound = node.getLowerBound();
		     if (node.getParameter() < _lower_bound) {
		         return;
		     }
		}


		// define and get branching edge
		int i = 0;
		int j = 0;
		double branch_nr = node.getEdgeForBranching(i, j);

		// if no edge exists to branch for you found a solution
		if (branch_nr == -1) {

		     // if not just maximum solution is requested or no solution so far...just add solution
		     if (!_just_max || _solutions.size() == 0 ) {
			  
			  // create a new solution
                          component_solution_type solution = { GraphSet(node.getGraph()), node.getChangedEdges(), node.getParameter(), node.getStartParameter()};

			 // add solution
 		          _solutions.insert(_solutions.end(), solution);

		     // otherwise just add solution if it is better then the one before & reject the old solution
		     } else if (_solutions[0].parameter < node.getParameter()) {
		          // creatre solution
		          component_solution_type solution = { GraphSet(node.getGraph()), node.getChangedEdges(), node.getParameter(), node.getStartParameter()};

		          // replace solution
		          _solutions[0] = solution;
		     }
		     // now stop recursion here
		     return;
		}
		
				
		// otherwise try to find a solution first

		// decide if branch should go first to the left or to the right
		double first = permanent;
		double second = forbidden;
		if (_lower_bound + node.getMinInsertCosts(i,j) > node.getParameter()) {
		     first = forbidden;
		     second = permanent;
		}

		// start recursive calls
		// case B1 (set (u,v) to forbidden, which decreases parameter
		WeightedProblemInstance *PI;
		try {
		     PI = node.clone();
		     PI->setEdge(i, j, first);

		     // maybe check if a stronger reduction is helpful
		     if (depth % 7 == 6) {
		          PI->strongReduce();
		     } else {
		          PI->reduce();
		     };

		     // test if graph is divided and start recursion
		     CostsGraph::vertex_matrix_type vertex_matrix = (PI->getGraph()).getConnectedVertices();
		     if (_split && splitWorthIt(vertex_matrix)) {
		         splitSearch(*PI, vertex_matrix, depth);
		     } else {
		         recSearch(*PI, depth+1);
		     }

		} catch ( ProblemInstanceException e ) {
		}
		delete PI;
		
		// case B2 (set (u,v) = (j,k) to permanent, (v,w) = (i,j) to forbidden
		// and (u,w) = (i,k) to forbidden, which decreases parameter
		try {
		     PI = node.clone();
		     PI->setEdge(i, j, second);

		     // check whether a strong reduction is helpful
		     if (depth % 7 == 6) {
		          PI->strongReduce();
		     } else {
		          PI->reduce();
		     };

		     // test if graph is divided and start recursion
		     CostsGraph::vertex_matrix_type vertex_matrix = (PI->getGraph()).getConnectedVertices();
		     if (_split && splitWorthIt(vertex_matrix)) {
		          splitSearch(*PI, vertex_matrix, depth);
		     } else {
		          recSearch(*PI, depth+1);
		     }
		} catch ( ProblemInstanceException e ) {
		}
		delete PI;
		
	}
}
