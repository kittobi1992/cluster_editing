/*+++++++++++++++++ class AlmostClique +++++++++++++++++++ */
#ifndef __almostclique_h__
#define __almostclique_h__

// include for IO for testing !!!
#include <sstream>

// includes from stlib
#include <vector>
#include <algorithm>

#include <costsgraph.h>
#include <array.h>
#include <trianglematrix.h>
#include <lowerbound.h>
#include <mincut.h>

/* This static class implements the almost clique rule.
   It therefore uses the min-cut class to obtain the min-cuts
   of almost cliques. It seeks always for highly connected
   components in a greedy way and checks for all of them whether
   we can apply our rule or not.

   It can be used parameter indpendent or also dependent.
   It returns a set of cliques which can be merged without any doubt.
*/

class AlmostClique {
public:

	AlmostClique();

	~AlmostClique();

	// compare function to compare two double arrays by its second entry
	static inline bool compare(Array<double> t1, Array<double> t2) {
		return (t1[1]>t2[1]);
	}


	// compare function to compare two integer arrays by its second entry
	static inline bool compare2(Array<int> t1, Array<int> t2) {
		return (t1[1]>t2[1]);
	}


	// returns a set of cliques and can be used parameter independent or dependent
	// this is the main function an implements the greedy search to highly connected components
	static inline bool get(const CostsGraph &graph, std::vector< CostsGraph::byte_vector_type> &cliques, double parameter, double parameter_dependent)
	{
		// if we use it as a parameter dependend rule then first get a set of independent conflict triples
		std::vector< Array<double> > triple_list;
		if (parameter_dependent) {
		     triple_list = LowerBound::getIndependentKT(graph);
		}

		int size = graph.getSize();

		// init vertex degree counting array/vector
		Array<int> empty = Array<int>(2);
		Array<double> empty2 = Array<double>(2);
		std::vector< Array<int> > vertex_degree = std::vector< Array<int> >(size);
		std::vector< Array<double> > edge_weights = std::vector< Array<double> >(size);
		for (int i = 0; i < size; i++) {
		     empty[0] = i;
		     empty[1] = 0;
		     empty2[0] = i;
		     empty2[1] = 0;
		     vertex_degree[i] = empty;
		     edge_weights[i] = empty2;
		}

		// calculate vertex degrees and maximum vertex
		for (int i = 0; i < size; i++) {
		     for (int k = 0; k < i; k++) {
	                  if (graph.getEdge(i,k) > 0) {
		               vertex_degree[i][1]++;
		               vertex_degree[k][1]++;
		               edge_weights[i][1] += graph.getEdge(i,k);
		               edge_weights[k][1] += graph.getEdge(i,k);
		          }
		     }
		}

		// calculate average edge weight to a vertex
		for (int i = 0; i < size; i++) {
			edge_weights[i][1] = (vertex_degree[i][1] != 0) ? (edge_weights[i][1] / vertex_degree[i][1]) : 0;
		}

		// fill array of average edge weights
		std::vector< Array<double> > ordered_vertex_degree = edge_weights;
		// sort it by their edge weights
		sort(ordered_vertex_degree.begin(), ordered_vertex_degree.end(), &compare);
		
		// save already used vertices and initiate it with false
		std::vector< bool > used_vertices = std::vector< bool >(size, false);

		// get vertex with greates degree=average edge weight and set it as used
		int max_degree_vertex = static_cast<int>(ordered_vertex_degree[0][0]);
		used_vertices[max_degree_vertex] = true;

		// as long we haven't used all vertices (we find a vertex with maximal average degree)
	        while (max_degree_vertex != -1) {

		     // iterate over diferent acceptable not connectebilities
		     // z=x means we accept that a node is not connected to at most x vertics in the already obtained vertex set
		     for (int z = 1; z >= 0; z--) {

		          // set start vertex as initial clique set
		          CostsGraph::byte_vector_type clique = CostsGraph::byte_vector_type(1, max_degree_vertex);

		          // seek for vertices to extend the clique
		          bool extended = true;
		          while (extended) {

		              extended = false;

		              // vector for clique sufficient clique neighbors 0=index 1=degree
		              std::vector< Array<double> > clique_neighbors = std::vector< Array<double> >(0);

		              // fill neighbor set
		              int c_member = 0;
		              // iterate over all vertices which are not yet in clique
		              for (int i = 0; i < size; i++) {
		                   if (c_member < clique.size() && clique[c_member] == i) {
		                       c_member++;
		                   } else {
		                       // check connectivity for vertex i
		                       int not_connected = 0;
		                       for (int k = 0; k < clique.size(); k++) {
		                            if (graph.getEdge(i,clique[k]) <= 0) {
	                                        if (graph.getEdge(i, clique[k]) == graph.forbidden) {
		                                     not_connected += size;
		                                } else {
	                                             not_connected++;
		                                }
		                            }
		                       }

				       // if vertex is connected strong enough to clique (less then z missings)
		                       if (not_connected <= z && not_connected < clique.size()) {
		                            Array<double> neighbor = Array<double>(2);
		                            neighbor[0] = i;

		                            // calculate weight to other clique members
		                            neighbor[1] = 0;
		                            for (int k = 0; k < clique.size(); k++) {
		                                 neighbor[1] += graph.getEdge(i, clique[k]);
		                            }
	                                    clique_neighbors.insert(clique_neighbors.end(), neighbor);
		                       }
		                   }
		              }

		              // if potential new set members are found....
		              if (clique_neighbors.size() > 0) {

		                  // order neighbor vector by weight of connectivity
		                  sort( clique_neighbors.begin(), clique_neighbors.end(), &compare);

		                  // add vertex to set clique and sort clique
		                  clique.insert(clique.end(), static_cast<int>(clique_neighbors[0][0]));
		                  sort( clique.begin(), clique.end() );

		                  // save the extension
		                  extended = true;

		                  // save that we the used vertex
		                  used_vertices[static_cast<int>(clique_neighbors[0][0])] = true;

		                  // if the next vertex in the ordered list was less than 0.9 connected to the set, then test almostCliqueRule
			          if ((clique_neighbors.size() == 1 || clique_neighbors[0][1] * 0.5 > clique_neighbors[1][1]) && almostCliqueRule(clique, graph, parameter, parameter_dependent, triple_list)) {
		                       cliques.insert(cliques.end(), clique);
		                       z = -1;
		                       extended = false;
		                  } 
		              }

		          }
		     }

		     // get new starting vertex
		     max_degree_vertex = -1;
		     for (int i = 0; i < size; i++) {

		         // check if we haven't used this vertex
		         if (!used_vertices[static_cast<int>(ordered_vertex_degree[i][0])]) {
		              max_degree_vertex = static_cast<int>(ordered_vertex_degree[i][0]);
		              used_vertices[static_cast<int>(ordered_vertex_degree[i][0])] = true;
	                      i = size+4;
		         }
		     }
		}

		// if we found cliques then return true, otherwise false
		if (cliques.size() > 0) {
		     return true;
		} else {
		     return false;
		}
	}
	

	// function which actually calculates the almostCliqueRule for a given vertex set
	// can be used parameter-independent and dependent
	inline static bool almostCliqueRule(CostsGraph::byte_vector_type clique, const CostsGraph &graph, double parameter, bool parameter_dependent, std::vector< Array<double> > &triple_list) 
	{
	        int size = graph.getSize();
		if (clique.size() > 1) {

		     // copy the clique into a new graph
		     CostsGraph clique_graph = CostsGraph(clique.size());
		     for (int i = 0; i < clique.size(); i++) {
		          for (int k = 0; k < i; k++) {
		               clique_graph.setEdge(i, k, graph.getEdge(clique[i], clique[k]));
		          }
		     }

		     // obtain minimal separation costs
		     double min_separation = MinCut::getGraphMinCut(clique_graph);

		     // check wether clique is heavy enough or not
		     // get edges to be inserted and cheapest edges in clique
	             //std::vector <double> edges = std::vector< double >(0);
		     //std::vector <int> vertex_degree = std::vector <int>(clique.size(), 0);
		     //std::vector <double> weight_to_neighbors = std::vector <double>(clique.size(), 0);

		     double edges_to_insert = 0.0;
		     double environment_weight = 0.0;

		     for (int i = 0; i < clique.size(); i++) {
		          //for (int k = 0; k < i; k++) {
		             /* if (graph.getEdge(clique[i], clique[k]) <= 0) {
	                          // save edges to be inserted
		                  edges_to_insert += graph.getEdge(clique[i], clique[k]) * -1;
		              } else {
		                  // save present edges to find cheapest edge later on
		                  edges.insert(edges.end(), graph.getEdge(clique[i], clique[k]));
		                  // increase degree
		                  vertex_degree[i]++;
		                  vertex_degree[k]++;
		                  // save weight of neighbors
		                  weight_to_neighbors[i] += graph.getEdge(clique[i],clique[k]);
		                  weight_to_neighbors[k] += graph.getEdge(clique[i],clique[k]);
		              }
		          }*/
		               // check environment
		               int c_member = 0;

	                       //double help = 0.0;
		               for (int k = 0; k < size; k++) {
		                   if (c_member < clique.size() && clique[c_member] == k) {
				       if (graph.getEdge(clique[i], k) <= 0) {
	                                    // save edges to be inserted
		                            edges_to_insert += graph.getEdge(clique[i], k) * -1;
				       }
		                       c_member++;
		                   } else {
		                       if (graph.getEdge(clique[i], k) > 0) {
		                          environment_weight += graph.getEdge(clique[i], k);
		                          //help += graph.getEdge(clique[i], k);
		                          if (!parameter_dependent && min_separation > environment_weight) {
		                              return false;
		                          }
		                       }
		                   }
		               }
		     }


		     // if rule is used as parameter dependent rule
		     if (parameter_dependent && clique.size() > graph.getSize()/10) {

		          // sum of conflict triples outside the clique
		          double bound_after_sep = 0.0;
		          for (int i = 0; i < triple_list.size(); i++) {
		               bool match = false;
		               for (int k = 0; k < clique.size(); k++) {
		                    if (triple_list[i][1] == clique[k] || triple_list[i][2] == clique[k] || triple_list[i][3] == clique[k]) {
		                        match = true;
		                    }
		               }
		               if (!match) {
		                    bound_after_sep += triple_list[i][0];
		               }
		          }

		          // check if lower bound + separations costs is larger than parameter
		          if (min_separation + bound_after_sep > parameter) {
		               return true;
		          }
		     }

		     // now make final almost clique check (greater than to avoid errors due to wrong roundings)
		     if (min_separation > environment_weight + edges_to_insert) {
		          return true;
		     }

		}
		return false;
	}
};


#endif
