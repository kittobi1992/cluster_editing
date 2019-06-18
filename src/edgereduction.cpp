/*+++++++++++++++++ class EdgeReduction +++++++++++++++++++ */
#include <edgereduction.h>

/* ##################### Constructors ###################### */

// constructor for as many edge lists as new the size of the parameter
// in egde list i are all elements inserted with costs i <= c < (i+1)
EdgeReduction::EdgeReduction(CostsGraph *graph, double parameter) : _graph(graph), _parameter(parameter), _min_insert_costs(_graph->getSize(),0), _min_delete_costs(_graph->getSize(),0), forbidden(_graph->forbidden), permanent(_graph->permanent)
{
	createCostsMatrices();	
}


// copy constructor
EdgeReduction::EdgeReduction(const EdgeReduction& ER, CostsGraph* graph) : _graph(graph), _parameter(ER._parameter), _min_insert_costs(ER._min_insert_costs), _min_delete_costs(ER._min_delete_costs), forbidden(ER.forbidden), permanent(ER.permanent)
{
}


// copy constructor which gets elements as parameters
inline EdgeReduction::EdgeReduction(CostsGraph *graph, double parameter, triangle_matrix_type m_i_c, triangle_matrix_type m_d_c) : _graph(graph), _parameter(parameter), _min_insert_costs(m_i_c), _min_delete_costs(m_d_c), forbidden(_graph->forbidden), permanent(_graph->permanent)
{
}


/* ###################### Destructor ###################### */
EdgeReduction::~EdgeReduction()
{
}

/* ###################### private/help functions ########## */

/* ---------------------- constructor help functions ------ */


// returns min insert and delete costs of (i,k) in combination
// with aditional vertex n by calling getInDelCosts with the 
// two edge values
inline void EdgeReduction::getInDelCosts(int i, int k, int n, double &sum_insert, double &sum_delete)
{
	double e1 = _graph->getEdge(i, n);
	double e2 = _graph->getEdge(k, n);

 	getInDelCosts(e1, e2, sum_insert, sum_delete);
}


// returns min insert and delete costs of an edge incident to the
// two parameter edges e1 and e2, under the condition, that both edges
// e1 and e2 are incident to the same node at there other end
inline void EdgeReduction::getInDelCosts(double e1, double e2, double &sum_insert, double &sum_delete)
{
	// both are positiv --> just delete costs will be added
	if (e1 > 0 && e2 > 0) {
	          sum_delete += std::min(e1, e2);
	// one is negative --> min of add and del operation will be added
	} else if (e1 > 0 || e2 > 0) {
		if (e1 > 0) {
		     if (e2 != forbidden) {
		          sum_insert += std::min(e1,e2 * -1);
		     } else {
		          sum_insert += e1;
		     }
		} else {
		     if (e1 != forbidden) {
		          sum_insert += std::min(e2, e1 * -1);
		     } else {
		          sum_insert += e2;
		     }
		}
	}

}


// returns min insert and delete costs for edge (i,k) with i>k
// therefore it calls getInDelCosts which calculates insert
// and delete costs for each triple combination of (i,k) with
// another vertex n
inline void EdgeReduction::calculateCostsForEdge(int i, int k, double &sum_insert, double &sum_delete)
{
	int size = _graph->getSize();

	for (int n = 0; n < k; n++) {
		getInDelCosts(i, k, n, sum_insert, sum_delete);
	}
	for (int n = k + 1; n < i; n++) {
		getInDelCosts(i, k, n, sum_insert, sum_delete);
	}
	for (int n = i + 1; n < size; n++) {
		getInDelCosts(i, k, n, sum_insert, sum_delete);
	}
}


// fills triangular matrices with insert and delete costs
// calls calculateCostsForEdge to fill each matrix element
inline void EdgeReduction::createCostsMatrices()
{
	int size = _graph->getSize();

	// fill non-common & common-neighbor matrices
	for (int i = 0; i < size; i++) {
		for (int k = 0; k < i; k++) {

		     double sum_insert = 0;
		     double sum_delete = 0;

			     calculateCostsForEdge(i, k, sum_insert, sum_delete);

		     _min_insert_costs.dirpos(i,k) = sum_insert + ((_graph->getEdge(i,k) <= 0) ? (_graph->getEdge(i,k) * -1) : 0);
		     _min_delete_costs.dirpos(i,k) = sum_delete + ((_graph->getEdge(i,k) > 0) ? _graph->getEdge(i,k) : 0);
		
		     if (_min_insert_costs.dirpos(i,k) > _parameter+0.01 && _min_delete_costs.dirpos(i,k) > _parameter+0.01) {
		        //std::cout << _parameter+0.01 << " " << _min_insert_costs.dirpos(i,k) << " " << _min_delete_costs.dirpos(i,k) << std::endl;
			throw ProblemInstanceException( " problem kernel reduction within given parameter range not possible " );
		     }
		}
	}
}



/* ---------------------- general help functions ------------------ */


// since the matrices are triangular this function helps to get the min deletion costs
double EdgeReduction::getMinDeleteCosts(unsigned short i, unsigned short j) const
{
	return _min_delete_costs.pos(i, j);
}


// since the matrices are triangular this function helps to get the min insertion costs
double EdgeReduction::getMinInsertCosts(unsigned short i, unsigned short j) const
{
	return _min_insert_costs.pos(i, j);
}


// since the matrices are triangular this function helps to set the min deletion costs
inline double EdgeReduction::setMinDeleteCosts(unsigned short i, unsigned short j, double value)
{
	_min_delete_costs.pos(i,j) = value;
}


// since the matrices are triangular this function helps to set the min insertion costs
inline double EdgeReduction::setMinInsertCosts(unsigned short i, unsigned short j, double value)
{
	_min_insert_costs.pos(i,j) = value;
}


// updates dependend on variable 'type' min insert or delete costs
// to the given value and will copy edge to the appropriate
// insertion or deletion list
inline void EdgeReduction::updateEdgeReference(int i, int j, char type, double value)
{
	if (i < j) {
		std::swap(i, j);
	}

	if (_graph->getEdge(i, j) != forbidden) {

		// dependend on 'type' its updating min insert or delete costs
		if ( type == 'i') {
		     _min_delete_costs.dirpos(i,j) = value;
		} else {
		     _min_insert_costs.dirpos(i,j) = value;
		}
	}
}


// checks whether its cheaper to delete a or b and return costs for this operation
inline double EdgeReduction::minDelEdge(double a, double b)
{
	if (a == forbidden) {
		return 0.0;
	} else {
		if (b == forbidden) {
		     return 0.0;
		} else {
		     return std::min(std::abs(a), std::abs(b));
		}
	}
}


// check whether its cheaper to insert a, or delete b and returns costs
inline double EdgeReduction::minInDelEdge(double a, double b)
{
	if (a == forbidden) {
		if (b == forbidden) {
		     // if b is already forbidden, the costs are zero
		     return 0.0;
		} else {
		     return std::abs(b);
		}
	} else {
		if (b == forbidden) {
		     // if b is forbidden its cheap = 0 to delete it
		     return 0.0;
		} else {
		     return std::min(std::abs(a), std::abs(b));
		}
	}
}

// set edge to new value (permanent or forbidden)
// in any case the min insert and delete costs will be updated and
// consequently the egdes will be copied updated to other edge lists
// also the parameter will be decreased and list with corresponding
// costs bigger then the new parameter will be concatenated to the
// last deletion and insertion list
void EdgeReduction::setEdgeToValue(int i, int j, double old_value, double parameter, double parameter_before)
{
	// get value and calculate list of a value equal to the current parameter
	double value = _graph->getEdge(i, j);

	_parameter = parameter;

	// case 1: edge was set and will be not_set
	if (old_value > 0 && value <= 0) {
		
		// iterate over all vertices incident to (i,j)
		for (int k = 0; k < _graph->getSize(); k++) {

		      // update for every other node the new neighbor situation for i and j
		      if (k != i && k != j) {

		          // get edge values to neighbor
		          double edge1 = _graph->getEdge(i, k);
		          double edge2 = _graph->getEdge(j, k);

		           // case 2.1: k is a common neighbor of i and j
		          if( edge1 > 0 && edge2 > 0 ){
		               updateEdgeReference(j,k,'i', getMinDeleteCosts(j,k) - minDelEdge( old_value, _graph->getEdge(i,k) ) );
		               updateEdgeReference(j,k,'d', getMinInsertCosts(j,k) + _graph->getEdge(i,k));
		               updateEdgeReference(i,k,'i', getMinDeleteCosts(i,k) - minDelEdge( old_value, _graph->getEdge(j,k) ) );
		               updateEdgeReference(i,k,'d', getMinInsertCosts(i,k) + _graph->getEdge(j,k));
		          } else if( edge1 > 0 && edge2 <= 0 ){
		          // case 2.2: k is neighbor of j but not of i
		               updateEdgeReference(j,k,'i',  getMinDeleteCosts(j,k) - minDelEdge( old_value, _graph->getEdge(i,k) ) );
		               updateEdgeReference(j,k,'d', getMinInsertCosts(j,k) +  _graph->getEdge(i,k));
		               updateEdgeReference(i,k,'d', getMinInsertCosts(i,k) - minInDelEdge(_graph->getEdge(j,k), old_value) );	
		          } else if( edge1 <= 0 && edge2 > 0 ){
		          // case 2.3: k is neighbor of i but not of j
		               updateEdgeReference(i,k,'i', getMinDeleteCosts(i,k) - minDelEdge( old_value, _graph->getEdge(j,k) ) );
		               updateEdgeReference(i,k,'d', getMinInsertCosts(i,k) + _graph->getEdge(j,k));
		               updateEdgeReference(j,k,'d', getMinInsertCosts(j,k) - minInDelEdge(_graph->getEdge(i,k), old_value) );	
		          } else {
		          // case 2.4: k is neither neighbor of i and j
		               updateEdgeReference(j,k,'d', getMinInsertCosts(j,k) - minInDelEdge(_graph->getEdge(i,k), old_value) );
		               updateEdgeReference(i,k,'d', getMinInsertCosts(i,k) - minInDelEdge(_graph->getEdge(j,k), old_value) );
		          }
		      }
		}

	}

	// additional case, which occurs not in unweighted case
	// if edge was not set and just the value will be changed
	if (old_value <= 0 && value == forbidden) {

		// iterate over all k incident to (i,j)
		for (int k = 0; k < _graph->getSize(); k++) {

		     if (k != i && k != j) {

		          double edge1 = _graph->getEdge(i,k);
		          double edge2 = _graph->getEdge(j,k);

		          if (edge1 > 0 && edge2 > 0) {
		               if (old_value * -1 < edge2 && edge2 != permanent && edge1 != forbidden && edge1 != permanent) {
		                    updateEdgeReference(i, k, 'd', getMinInsertCosts(i, k) + old_value + edge2);
		               }
		               if (old_value * -1 < edge1 && edge1 != permanent && edge2 != forbidden && edge2 != permanent) {
		                    updateEdgeReference(j, k, 'd', getMinInsertCosts(j, k) + old_value + edge1);
		               }
		          } else if (edge1 <= 0 && edge2 > 0) {
		               if (old_value * -1 < edge2 && edge2 != permanent && edge1 != forbidden) {
		                    updateEdgeReference(i, k, 'd' , getMinInsertCosts(i, k) + old_value + edge2);
		               }
		          } else if (edge1 > 0 && edge2 <= 0) {
		               if (old_value * -1 < edge1 && edge1 != permanent && edge2 != forbidden ) {
		                    updateEdgeReference(j, k, 'd', getMinInsertCosts(j, k) + old_value + edge1);
		               }
		          }
		     }
		}
	}

}


// will consider the update of the edge (i,j) if a and b were merged
// therefore it calculates the old min insert and delete costs considering a and b
// and calculates the new costs of edge (i,j) by getting the new costs of the
// edge i,b and j,b since a is the vertex which was deleted
// it returns the updated costs by adding the  relative change of 
// the min insert and delete costs for the given edge
inline void EdgeReduction::updateEdgeCase1(int i, int j, int a, int b, double new_costs1, double new_costs2, double &up_insert, double &up_delete)
{
	double old_insert = 0;
	double old_delete = 0;
	getInDelCosts(i, j, a, old_insert, old_delete);
	getInDelCosts(i, j, b, old_insert, old_delete);
	
	double new_insert = 0;
	double new_delete = 0;
	getInDelCosts(new_costs1, new_costs2, new_insert, new_delete);
	
	up_insert += new_insert - old_insert;
	up_delete += new_delete - old_delete;
}

// considers the second case, that edges (i,b) is calculated from scratch
// therefore it calculates the min insert and delete costs of this edge
// considering the vertex j
// therefore it calculates the min costs for del/insert the edges (i,j)
// and (b,j)
inline void EdgeReduction::updateEdgeCase2(int i, int j, int a, int b, double new_costs1, double new_costs2, double &up_insert, double &up_delete)
{
	double new_insert = 0;
	double new_delete = 0;
	getInDelCosts(_graph->getEdge(i,j), new_costs2, new_insert, new_delete);

	up_insert += new_insert;
	up_delete += new_delete;
}


// function updates the min(del/insert)costs matrices after a and b were merged in O(n^2)
// For all i: the edge(i,b) will be calculated completly new = O(n*n)
// For all i,k (not a,b): edge(i,k) will be updated in 2 steps with the new values to a and b = O(n^2 * 2)
// note that size is still the original size
inline void EdgeReduction::updateCostsMatrices(int a, int b, double_array_type new_costs)
{

	// get graph size (should be the old size, before (a,b) is really merged
	int size = _graph->getSize();

	_min_insert_costs.deleteIndex(a);
	_min_delete_costs.deleteIndex(a);

	// set (i,b) to zero or to edge, since it will be calculated again from scratch
	for (int i = 0; i < b; i++) {
		_min_insert_costs.dirpos(b, i) = ((new_costs[i] <= 0) ? (new_costs[i] * -1) : 0);
		_min_delete_costs.dirpos(b, i) = ((new_costs[i] > 0) ? new_costs[i] : 0);
	}
	for (int i = b + 1; i < a; i++) {
		_min_insert_costs.dirpos(i, b) = ((new_costs[i - 1] <= 0) ? (new_costs[i - 1] * -1) : 0);
		_min_delete_costs.dirpos(i, b) = ((new_costs[i - 1] > 0) ? new_costs[i - 1] : 0);
	}
	for (int i = a + 1; i < size; i++) {
		_min_insert_costs.dirpos(i - 1, b) = ((new_costs[i - 2] <= 0) ? (new_costs[i - 2] * -1) : 0);
		_min_delete_costs.dirpos(i - 1, b) = ((new_costs[i - 2] > 0) ? new_costs[i - 2] : 0);
	}

	// note to the indices of new_costs:
	// if x < b then its simple at pos x
	// if b < x < a then its x-1 because position b was skipped in the array
	// if x > a then its x-2 cause 1 position was skipped and a was deleted

	// iterate over all (i,k) and update edges
	// updated (i,k) edges will be copied directly after update to the appropriate list
	for (int i = 0; i < b; i++) {
		for (int k = 0; k < i; k++) {
		     // update (i,k)
		     updateEdgeCase1(i, k, a, b, new_costs[i], new_costs[k], _min_insert_costs.dirpos(i, k), _min_delete_costs.dirpos(i, k));
		     // update i,ab
		     updateEdgeCase2(i ,k, a, b, new_costs[i], new_costs[k], _min_insert_costs.dirpos(b, i), _min_delete_costs.dirpos(b, i));
		     updateEdgeCase2(k ,i, a, b, new_costs[k], new_costs[i], _min_insert_costs.dirpos(b, k), _min_delete_costs.dirpos(b, k));
		}
	}
	for (int i = b+1; i < a; i++) {
		for (int k = 0; k < b; k++) {
		     // update i,k
		     updateEdgeCase1(i, k, a, b, new_costs[i - 1], new_costs[k], _min_insert_costs.dirpos(i, k), _min_delete_costs.dirpos(i, k));
		     // update i,ab
		     updateEdgeCase2(i ,k, a, b, new_costs[i - 1], new_costs[k], _min_insert_costs.dirpos(i, b), _min_delete_costs.dirpos(i, b));
		     updateEdgeCase2(k ,i, a, b, new_costs[k], new_costs[i - 1], _min_insert_costs.dirpos(b, k), _min_delete_costs.dirpos(b, k));
		}
		for (int k = b+1; k < i; k++) {
		     // update i,k
		     updateEdgeCase1(i, k, a, b, new_costs[i - 1], new_costs[k - 1], _min_insert_costs.dirpos(i, k), _min_delete_costs.dirpos(i,k));
		     // update i,ab
		     updateEdgeCase2(i ,k, a, b, new_costs[i - 1], new_costs[k - 1], _min_insert_costs.dirpos(i, b), _min_delete_costs.dirpos(i, b));
		     updateEdgeCase2(k ,i, a, b, new_costs[k - 1], new_costs[i - 1], _min_insert_costs.dirpos(k, b), _min_delete_costs.dirpos(k, b));
		}
	}
	for (int i = a+1; i < size; i++) {
		for (int k = 0; k < b; k++) {
		     // update i,k
		     updateEdgeCase1(i, k, a, b, new_costs[i - 2], new_costs[k], _min_insert_costs.dirpos(i - 1, k), _min_delete_costs.dirpos(i - 1, k));
		     // update i,ab
		     updateEdgeCase2(i ,k, a, b, new_costs[i - 2], new_costs[k], _min_insert_costs.dirpos(i - 1, b), _min_delete_costs.dirpos(i - 1, b));
		     updateEdgeCase2(k ,i, a, b, new_costs[k], new_costs[i - 2], _min_insert_costs.dirpos(b, k), _min_delete_costs.dirpos(b, k));
		}
		for (int k = b + 1; k < a; k++) {
		     // update i,k
		     updateEdgeCase1(i, k, a, b, new_costs[i - 2], new_costs[k - 1], _min_insert_costs.dirpos(i - 1, k), _min_delete_costs.dirpos(i - 1, k));
		     // update i,ab
		     updateEdgeCase2(i ,k, a, b, new_costs[i - 2], new_costs[k - 1], _min_insert_costs.dirpos(i - 1, b), _min_delete_costs.dirpos(i - 1, b));
		     updateEdgeCase2(k ,i, a, b, new_costs[k - 1], new_costs[i - 2], _min_insert_costs.dirpos(k, b), _min_delete_costs.dirpos(k, b));
		}
		for (int k = a + 1; k < i; k++) {
		     // update i,k
		     updateEdgeCase1(i, k, a, b, new_costs[i - 2], new_costs[k - 2], _min_insert_costs.dirpos(i - 1, k - 1), _min_delete_costs.dirpos(i - 1, k - 1));
		     // update i,ab
		     updateEdgeCase2(i ,k, a, b, new_costs[i - 2], new_costs[k - 2], _min_insert_costs.dirpos(i - 1, b), _min_delete_costs.dirpos(i - 1, b));
		     updateEdgeCase2(k ,i, a, b, new_costs[k - 2], new_costs[i - 2], _min_insert_costs.dirpos(k - 1, b), _min_delete_costs.dirpos(k - 1, b));
		}		
	}
}


// will fullfill the complete merge operation of the vertices i and j
// it will update the underlying graph with the new edge values,
// the min insert and delete costs matrices, the conflict tripel and
// might decrease the parameter, since the merge operation can trigger
// new immediate changes
void EdgeReduction::mergeVertices(int i, int j, double parameter, double_array_type new_costs)
{
	double parameter_before = _parameter;
	_parameter = parameter;

	// update min insert and delete costs
	updateCostsMatrices(i, j, new_costs);
	
}


/* ################### main/reduce function ################### */


// if instance was created by empty copy constructor, initiate matrices
void EdgeReduction::init(CostsGraph* graph, double parameter, bool merge)
{
	_graph = graph;
	_parameter = parameter;
	_min_insert_costs = triangle_matrix_type(_graph->getSize(),0);
	_min_delete_costs = triangle_matrix_type(_graph->getSize(),0);
	forbidden = _graph->forbidden;
	permanent = _graph->permanent;

	createCostsMatrices();
}


// check for edges to delete and return them
inline CostsGraph::edge_list_type EdgeReduction::getEdgesToDelete(bool strong, triangle_matrix_type &lower_bound_edge)
{
	CostsGraph::edge_list_type edge_list = CostsGraph::edge_list_type(0);
	int size = _graph->getSize();
	
	// iterate over all edges
	for (int h = 0; h < size; h++) {
		for (int k = 0; k < h; k++) {

		     double edge = _graph->getEdge(h,k);

		     // check simple if icp is greater than parameter
		     if (edge != forbidden && _min_insert_costs.pos(h,k) > _parameter + 0.001) {
		          //std::cout << "(" << h << "," << k << ") simple mi=" << _min_insert_costs.pos(h,k) << std::endl;
		          CostsGraph::pair_type pair = {h, k};
		          edge_list.insert(edge_list.end(), pair);

		     // if not and you use strong reduce, include lower bound in the calculation
		     } else if (strong && edge != forbidden && (_min_insert_costs.pos(h,k) + lower_bound_edge.dirpos(h,k)) > _parameter + 0.001) {
			  //std::cout << "(" << h << "," << k << ") compl mi=" << _min_insert_costs.pos(h,k) << " lb=" << lower_bound_edge.dirpos(h,k) << std::endl;
		          CostsGraph::pair_type pair = {h, k};
		          edge_list.insert(edge_list.end(), pair);
	             };
		}
	}
	return edge_list;
}


// check for edges to insert and return them
inline CostsGraph::edge_list_type EdgeReduction::getEdgesToInsert(double min_cut_threshold, triangle_matrix_type &lower_bound_edge)
{
	CostsGraph::edge_list_type edge_list = CostsGraph::edge_list_type(0);

	int size = _graph->getSize();
	std::vector< Array<int> > fragile_edges = std::vector< Array<int> >(0);

	// iterate over all edges
	for (int h = 0; h < size; h++) {
		for (int k = 0; k < h; k++) {

		     double edge = _graph->getEdge(h,k);

		     if (edge != forbidden) {
		          // simply check if icf is greater than parameter
		          if (_min_delete_costs.pos(h,k) > _parameter + 0.001) {
		               CostsGraph::pair_type pair = {h, k};
		               edge_list.insert(edge_list.end(), pair);

		          // include lower bound in calculation if strong reduce is used
	                  } else if (min_cut_threshold > 0.0 && edge != forbidden && _min_delete_costs.pos(h,k) + lower_bound_edge.dirpos(h,k) > _parameter+0.001) {
		               CostsGraph::pair_type pair = {h, k};
		               edge_list.insert(edge_list.end(), pair);
		          // check for fragile edges, meaning close beyond parameter
		          }/* else if ((_min_delete_costs.pos(h,k)/_parameter) > (1.0 - min_cut_threshold)) {
		               Array<int> edge = Array<int>(2);
		               edge[0] = h;
		               edge[1] = k;
		               fragile_edges.insert(fragile_edges.end(), edge);
		          }*/ 
		     }
		}
	}
	/*
	// finally check fragile edges again using min_cut in icf
	if (edge_list.size() == 0) {
		for (int h = 0; h < fragile_edges.size(); h++) {
		     if (MinCut::getMinSTCut(fragile_edges[h][0], fragile_edges[h][1], *_graph) > _parameter) {
		         CostsGraph::pair_type pair = {fragile_edges[h][0], fragile_edges[h][1]};
		         edge_list.insert(edge_list.end(), pair);
		     }
		}
	}*/

	return edge_list;
}


// fill lower bound matrix
EdgeReduction::triangle_matrix_type EdgeReduction::createLowerBoundMatrix()
{		
	// get all independent conflict triples
	
	std::vector< Array<double> > triple_list = LowerBound::getIndependentKT(*_graph);

	// create triangle matrix to save the lower bound excluding the edge indicated by the index ij
	triangle_matrix_type lower_bound_edge = triangle_matrix_type(_graph->getSize(), 0.0);
	double sum_bound = 0.0;

	// iterate over all conflict triples
	// sum up all conflict triples to achieve a first lower bound
	// beside this delete all conflict triples with edge ij from ij
	for (int i = 0; i < triple_list.size(); i++) {
	     sum_bound += triple_list[i][0];
	     for (int k = 0; k < static_cast<int>(triple_list[i][1]); k++) {
	          lower_bound_edge.pos(static_cast<int>(triple_list[i][1]), k) -= triple_list[i][0];
	     }
	     for (int k = static_cast<int>(triple_list[i][1])+1; k < _graph->getSize(); k++) {
	          lower_bound_edge.pos(static_cast<int>(triple_list[i][1]), k) -= triple_list[i][0];
	     }
	     for (int k = 0; k < static_cast<int>(triple_list[i][2]); k++) {
	          lower_bound_edge.pos(static_cast<int>(triple_list[i][2]), k) -= triple_list[i][0];
	     }
	     for (int k = static_cast<int>(triple_list[i][2])+1; k < _graph->getSize(); k++) {
	          lower_bound_edge.pos(static_cast<int>(triple_list[i][2]), k) -= triple_list[i][0];
	     }
	     for (int k = 0; k < static_cast<int>(triple_list[i][3]); k++) {
	          lower_bound_edge.pos(static_cast<int>(triple_list[i][3]), k) -= triple_list[i][0];
	     }
	     for (int k = static_cast<int>(triple_list[i][3])+1; k < _graph->getSize(); k++) {
	          lower_bound_edge.pos(static_cast<int>(triple_list[i][3]), k) -= triple_list[i][0];
	     }
	}
	
	// than add first lower bound to evey ij
	// ij now indicates all independent conflict triples excluding edge ij
	for (int i = 0; i < _graph->getSize(); i++) {
	     for (int k = 0; k < i; k++) {
	          lower_bound_edge.dirpos(i, k) += sum_bound;
	     }
	}

	return lower_bound_edge;
}


// this is the main function
// it returns edges which need to be set to forbidden or permanent respectivley
bool EdgeReduction::getReduceableEdges(CostsGraph::edge_list_type &permanent_edges, CostsGraph::edge_list_type &forbidden_edges, double fuzzy_threshold)
{
	triangle_matrix_type lower_bound_edge;

	// create lower bound matrix
	if (fuzzy_threshold > 0.0) {
		lower_bound_edge = createLowerBoundMatrix();
	}

	// call above functions to get edges to delete and insert
	permanent_edges = getEdgesToInsert(fuzzy_threshold, lower_bound_edge);
	forbidden_edges = getEdgesToDelete(((fuzzy_threshold>0) ? true : false), lower_bound_edge);
	
	return (permanent_edges.size() + forbidden_edges.size() > 0) ? true : false;
}


// create deep copy of object
inline EdgeReduction* EdgeReduction::clone() const
{
	EdgeReduction* newObjPtr = new EdgeReduction(*this);
	return newObjPtr;
}

// divide this object under the constraint given in the vertex matrix
// normally each row corresponds to one connected component
// always all vertices in one row will packed in one new object
EdgeReduction::edge_red_list_type EdgeReduction::divideInstance(CostsGraph::vertex_matrix_type vertex_matrix, CostsGraph::graph_list_type graph_list)
{
	int cc_nr = vertex_matrix.size();

	// build costs matrices
	Array< triangle_matrix_type > m_i_l = Array< triangle_matrix_type >(cc_nr);
	Array< triangle_matrix_type > m_d_l = Array< triangle_matrix_type >(cc_nr);

	for (int i = 0; i < cc_nr; i++) {
		m_i_l[i] = triangle_matrix_type(vertex_matrix[i].size(), 0);
		m_d_l[i] = triangle_matrix_type(vertex_matrix[i].size(), 0);
		for (int k = 0; k < vertex_matrix[i].size(); k++) {
		     for (int h = 0; h < k; h++) {
		          m_i_l[i].pos(k,h) = getMinInsertCosts(vertex_matrix[i][k], vertex_matrix[i][h]);
		          m_d_l[i].pos(k,h) = getMinDeleteCosts(vertex_matrix[i][k], vertex_matrix[i][h]);
		     }
		}
	}	
	
	// create list of problem instances
	edge_red_list_type er_list = edge_red_list_type(cc_nr);

	// create copies
	for (int i = 0; i < cc_nr; i++) {
		er_list[i] = new EdgeReduction(graph_list[i], _parameter,  m_i_l[i], m_d_l[i]);
	}

	return er_list;
	
}


// set parameter manually to value
// will cause merging of list between old and new parameter
// and a reduction step will be called afterwards
void EdgeReduction::setParameter(double new_parameter)
{
	_parameter = new_parameter;
}


// return _min_insert_costs matrix (necessary for branching numbers)
EdgeReduction::triangle_matrix_type EdgeReduction::getMinInsertMatrix() const
{
	return _min_insert_costs;
}


// return _min_delete_costs matrix (necessary for branching numbers)
EdgeReduction::triangle_matrix_type EdgeReduction::getMinDeleteMatrix() const
{
	return _min_delete_costs;
}


// assignment operator
EdgeReduction& EdgeReduction::operator=(const EdgeReduction& er)
{
	_graph = er._graph;
	_min_insert_costs = er._min_insert_costs;
	_min_delete_costs = er._min_delete_costs;

	forbidden = er.forbidden;
	permanent = er.permanent;

	_parameter = er._parameter;
}

/* ################### Heuristic functions ############################### */

// return nr_edges with greates icp
CostsGraph::edge_list_type EdgeReduction::getTopEdgesToDelete(int nr_edges)
{
	triangle_matrix_type lower_bound_edge = createLowerBoundMatrix();

	int size = _graph->getSize();
	std::vector< edge_type > icp_list = std::vector< edge_type>(0);

	// iterate over all edges
	for (int h = 0; h < size; h++) {
		for (int k = 0; k < h; k++) {
		     double edge = _graph->getEdge(h,k);
		     // include lower bound in the calculation and add icp
		     if (edge != forbidden) {
		          edge_type new_edge = {h, k, _min_insert_costs.pos(h,k) + lower_bound_edge.dirpos(h,k)};
		          icp_list.insert(icp_list.end(), new_edge);
	             };
		}
	}

	sort(icp_list.begin(), icp_list.end(), &compare);

	// create final return list	
	CostsGraph::edge_list_type edge_list = CostsGraph::edge_list_type(0);
	
	// get nr_edges top edges
	if (nr_edges > icp_list.size()) {
		nr_edges = icp_list.size();
	}
	for (int h = 0; h < nr_edges; h++) {	
		CostsGraph::pair_type pair = {icp_list[h].i, icp_list[h].j};
		edge_list.insert(edge_list.end(), pair);
	}

	return edge_list;
}


// return nr_edges with greates icf
CostsGraph::edge_list_type EdgeReduction::getTopEdgesToInsert(int nr_edges)
{
	triangle_matrix_type lower_bound_edge = createLowerBoundMatrix();

	int size = _graph->getSize();
	std::vector< edge_type > icf_list = std::vector< edge_type>(0);

	// iterate over all edges
	for (int h = 0; h < size; h++) {
		for (int k = 0; k < h; k++) {
		     double edge = _graph->getEdge(h,k);
		     if (edge != forbidden) {
 		          edge_type new_edge = {h, k, _min_delete_costs.pos(h,k) + lower_bound_edge.dirpos(h,k)};
		          icf_list.insert(icf_list.end(), new_edge);
		     }
		}
	}

	sort(icf_list.begin(), icf_list.end(), &compare);

	CostsGraph::edge_list_type edge_list = CostsGraph::edge_list_type(0);

	// get nr_edges top edges
	if (nr_edges > icf_list.size()) {
		nr_edges = icf_list.size();
	}
	for (int h = 0; h < nr_edges; h++) {	
		CostsGraph::pair_type pair = {icf_list[h].i, icf_list[h].j};
		edge_list.insert(edge_list.end(), pair);
	}

	return edge_list;
}


/* ################### Programmer's help functions ####################### */


// print min insert and delete lists
void EdgeReduction::printLists()
{
	std::cout << "min_insert: " << std::endl;
	for (int i = 0; i < _graph->getSize(); i++) {
		for (int k = 0; k < i; k++) {
			std::cout << _min_insert_costs.pos(i,k) << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	std::cout << "min_delete: " << std::endl;
	for (int i = 0; i < _graph->getSize(); i++) {
		for (int k = 0; k < i; k++) {
			std::cout << _min_delete_costs.pos(i, k) << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}
