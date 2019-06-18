/*+++++++++++++++++ class MergeReduction +++++++++++++++++++ */
#ifndef __mergereduction_h__
#define __mergereduction_h__

// include for IO for testing !!!
#include <sstream>

// includes from stlib
#include <vector>
#include <math.h>

#include <array.h>
#include <costsgraph.h>

/* This class implements the merging vertices rule.
   It receives an edge which needs to be merged and returns
   the new adjacencies of the new vertex and the modification
   costs of this merging operation.
   The actual graph update is done within the
   weightedprobleminstance class.

   Note that this class is static and is used modular by function call.
*/

class MergeReduction {
public:
	
	MergeReduction();

	~MergeReduction();


	// merges virtually two vertices and
	// returns the new costs to all other vertices and the modifcation costs
	static inline double mergeVertices(int i, int j, const CostsGraph &graph, double forbidden, Array<double> &new_costs)
	{
		// save old graph size
		int old_size = graph.getSize();
	
		// create a vector just for the new edge value from the now merged ij and all other nodes
		new_costs = Array<double>(old_size - 2, 0);
	
		// save parameter decreasment during the update
		double to_pay = 0;
		
		// iterate over all possible other nodes k, to calculate all values (ij,k)
		for (int k = 0; k < j; k++) {
			double edge1 = graph.getEdge(i, k);
			double edge2 = graph.getEdge(j, k);
	
			if (edge1 <= 0 && edge2 <= 0) {
			if (edge1 == forbidden) {
				new_costs[k] = forbidden;
			} else if (edge2 == forbidden) {
				new_costs[k] = forbidden;
			} else {
				new_costs[k] = edge1 + edge2;
			}
			} else if (edge1 > 0 && edge2 > 0) {
			new_costs[k] = edge1 + edge2;
			} else if (edge1 > 0 && edge2 <= 0) {
			if (edge1 > (edge2 * -1)) {
				new_costs[k] = edge1 + edge2;
				to_pay += -1 * edge2;
			} else {
				if (edge2 == forbidden) {
				new_costs[k] = forbidden;
				} else {
				new_costs[k] = edge2 + edge1;
				}
				to_pay += edge1;
			}
			} else {
			if (edge2 > (edge1 * -1)) {
				new_costs[k] = edge2 + edge1;
				to_pay += -1 * edge1;
			} else {
				if (edge1 == forbidden) {
				new_costs[k] = forbidden;
				} else {
				new_costs[k] = edge1 + edge2;
				}
				to_pay += edge2;
			}
			}
		}
	
		// remark here that k is bigger then j and we will not save the incidences to j itself
		// so that all other indices are decreased by one
		for (int k = j + 1; k < i; k++) {
			double edge1 = graph.getEdge(i, k);
			double edge2 = graph.getEdge(j, k);
	
			if (edge1 <= 0 && edge2 <= 0) {
			if (edge1 == forbidden) {
				new_costs[k - 1] = forbidden;
			} else if (edge2 == forbidden) {
				new_costs[k - 1] = forbidden;
			} else {
				new_costs[k - 1] = edge1 + edge2;
			}
			} else if (edge1 > 0 && edge2 > 0) {
			new_costs[k - 1] = edge1 + edge2;
			} else if (edge1 > 0 && edge2 <= 0) {
			if (edge1 > (edge2 * -1)) {
				new_costs[k - 1] = edge1 + edge2;
				to_pay += -1 * edge2;
			} else {
				if (edge2 == forbidden) {
				new_costs[k - 1] = forbidden;
				} else {
				new_costs[k - 1] = edge2 + edge1;
				}
				to_pay += edge1;
			}
			} else {
			if (edge2 > (edge1 * -1)) {
				new_costs[k - 1] = edge2 + edge1;
				to_pay += -1 * edge1;
			} else {
				if (edge1 == forbidden) {
				new_costs[k - 1] = forbidden;
				} else {
				new_costs[k - 1] = edge1 + edge2;
				}
				to_pay += edge2;
			}
			}
		}
	
		// remark here that k is bigger then j and j will be removed, so that indices of nodes after j
		// are decreased by two now
		for (int k = i + 1; k < old_size; k++) {
			double edge1 = graph.getEdge(i, k);
			double edge2 = graph.getEdge(j, k);
	
			if (edge1 <= 0 && edge2 <= 0) {
			if (edge1 == forbidden) {
				new_costs[k - 2] = forbidden;
			} else if (edge2 == forbidden) {
				new_costs[k - 2] = forbidden;
			} else {
				new_costs[k - 2] = edge1 + edge2;
			}
			} else if (edge1 > 0 && edge2 > 0) {
			new_costs[k - 2] = edge1 + edge2;
			} else if (edge1 > 0 && edge2 <= 0) {
			if (edge1 > (edge2 * -1)) {
				new_costs[k - 2] = edge1 + edge2;
				to_pay += -1 * edge2;
			} else {
				if (edge2 == forbidden) {
				new_costs[k - 2] = forbidden;
				} else {
				new_costs[k - 2] = edge2 + edge1;
				}
				to_pay += edge1;
			}
			} else {
			if (edge2 > (edge1 * -1)) {
				new_costs[k - 2] = edge2 + edge1;
				to_pay += -1 * edge1;
			} else {
				if (edge1 == forbidden) {
				new_costs[k - 2] = forbidden;
				} else {
				new_costs[k - 2] = edge1 + edge2;
				}
				to_pay += edge2;
			}
			}
		}
		
		return to_pay;	
	}
};


#endif
