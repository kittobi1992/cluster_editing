/*+++++++++++++++++ class BNManager +++++++++++++++++++ */
#ifndef __bnmanager_h__
#define __bnmanager_h__

// include for IO for testing !!!
#include <sstream>

// includes from stlib
#include <vector>
#include <math.h>

#include <costsgraph.h>
#include <edgereduction.h>

/* This static class implements the search of the edge
   with the smallest branching number.
   It needs an edgereduction object to get the minimum merging costs = icp
   and a costsgraph to obtain the minimum deleting costs = edge itself.

   It finally returns the edge with minimal branching number.
*/

class BNManager {
public:
	/* ####  main functions  #### */

	// return the ege with minimal branching number and the branching number itself
	// -1 is returned if no branching number < inf exists.
	static inline double getMinBNEdge(int &i, int &j, const CostsGraph* graph, const EdgeReduction* edge_reduction) 
	{
		if (graph->getSize() <= 2) {
		     return -1;
		}

		double min = graph->permanent;

		// iterate over all edges
		for (int h = 0; h < graph->getSize(); h++) {
		     for (int k = 0; k < h; k++) {
	
		          // get branching number 
		          double bn = getLogBN(graph->getEdge(h,k)+0.000001, edge_reduction->getMinInsertCosts(h,k) - ((graph->getEdge(h,k)>0) ? 0 : -1 * graph->getEdge(h,k)), graph->permanent);

		          if ( bn > 0 && bn < min) {
		               min = bn;
		               i = h;
		               j = k;
		          }
		     }
		}

		return (min == graph->permanent) ? -1 : min;
	};


private:

	// since the matrices are triangular the swap functions swaps two integer indices
	static inline void swap(double &i, double &j) 
	{
		double help = i;
		i = j;
		j = help;
	}

	// calculate log branching number of (del, in costs)
	static inline double getLogBN(double a, double b, double permanent)
	{
		if (a <= 0.0 || b <= 0.0) {
		     return permanent;
		}

		if (a > b) {
		     swap(a, b);
		};

		double x = b / a;
		double y = 0.0;

		// these functions are fitted cubic splines to already calculated branching numbers, see thesis for details
		if (x <= 10.0) {
		     y = (9.12524 + x*(15.9154 + x*(2.69316 + x * 0.0157118))) / (4.13278 + x*(21.9904 + x*(12.911 + x)));
		} else {
		     y = (5.25471 + x*(74.3963 + x*(3.67704 + x* 0.00435418))) / (-87.1522 + x*(120.572 + x*(41.0281 + x)));
		}

		return (1/a)*y;
	};

};


#endif
