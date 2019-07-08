/*+++++++++++++++++ class criticalclique +++++++++++++++++++ */
#ifndef __criticalclique_h__
#define __criticalclique_h__

#include <iostream>
#include <sstream>

#include <algorithm>
#include <vector>
#include <math.h>

#include <sys/types.h>
#include <trianglematrix.h>
#include <Exceptions/graphexception.h>
#include <graphset.h>
#include <napsack.h>

using std::vector;
using std::string;
using std::cout;
using std::endl;

/* This static class implements the critical clique
   adaption to the weighted case, which works by dynamic programming.
   Given a graph the reduce function checks for vertex tuples with
   similar neighborhood and returns them in a "permanent" list.
  
   First the graph will be transformed into a integer weighted graph
   by first multiplying all edges and then rpund them to an int value.
   A Napsack like algorithm in class Napsack does the final
   dynamic programming. 
*/


class CriticalClique {
public:

    CriticalClique() = default;

    ~CriticalClique() = default;

    // function which give back a list of pairs of vertices which can be merged
    // it therefore computes creates a integer trianglematrix from the graph
    // computes a list of tuples and uses then dynamic programming
    inline static CostsGraph::edge_list_type reduce(CostsGraph CG, double roundValue) {
        // transform graph to a triangle matrix with greater values
        TriangleMatrix<long double> G = createTriangleMatrix(CG, roundValue);

        // create emtpy edge list
        CostsGraph::edge_list_type list(0);

        // save forbidden value as long long value to save tuples with inf values
        // divide it by the number of nodes to avoid overflow if forbidden is added
        long double forbidden = ((long long) (roundValue * CG.forbidden)) / G.size();

        // iterate over all edges in G
        for (int u = 0; u < G.size(); u++) {
            for (int v = u + 1; v < G.size(); v++) {

                if (CG.getEdge(u, v) >= 0.0000001) {

                    long double delta_u = 0;
                    long double delta_v = 0;
                    long double rel_differenz = 0;
                    long long X = 0;
                    long long Y = 0;

                    // calculate delta_u and delta_v and the tuple list
                    vector<Napsack::Tupel> B = computeValues(u, v, G, X, Y, delta_u, delta_v, rel_differenz, forbidden,
                                                             G.dirpos(u, v));


                    //std::cout << "check (" << u << "," << v << ") with du=" << delta_u << " dv=" << delta_v << " diff=" << rel_differenz << " tuples: " << std::endl;
                    //for (int i = 0; i < B.size(); i++) {
                    //	std::cout << B[i].x_up << " " << B[i].x_down << " " << B[i].y_up << " " << B[i].y_down << std::endl;
                    //  }

                    // check whether dynamic programming needs to be done
                    if (G.dirpos(v, u) > std::min(delta_u, delta_v)) {

                        // now do dynamic programming...
                        double maxmin = CG.forbidden;
                        long long X_index = 2 * X + 1;
                        long long Y_index = 2 * Y + 1;

                        //std::cout << "test trivial" << std::endl;
                        // if trivial approximation does not already hold
                        if (G.dirpos(v, u) <= 0.5 * (rel_differenz + delta_u + delta_v)) {
                            if (delta_u < CG.permanent && delta_v < CG.permanent) {
                                // calculate a bound for inf to avoid direct equality test
                                long long infinity = (long long) (-1 * forbidden) / 10;
                                // do the adapted napsack dynamic programming algorithm
                                //	    std::cout << "start nap" << std::endl;
                                maxmin = Napsack::getNapsackSolution(B, X_index, Y_index, delta_u, delta_v, G.pos(u, v),
                                                                     infinity);
                            } else {
                                maxmin = CG.permanent;
                            }
                        }

                        // check if uv needs to be set to permanent...if yes add it to the edge list
                        // if equal we get always a problem with two exclusive permanent edges
                        if (G.pos(u, v) > maxmin) {
                            //std::cout << "set to perm...edge is " << G.pos(u,v) << "  maxmin value = " << maxmin << std::endl;
                            CostsGraph::pair_type help = {v, u};
                            list.push_back(help);
                            if (list.size() > G.size() / 10 && list.size() > 10) {
                                return list;
                            }
                        }
                    }
                }

            }
        }

        return list;
    }

private:

    /* #### private help functions #### */


    // function to create a multiply a graph with a certain factor to obtain
    // triangle matrix
    inline static TriangleMatrix<long double> createTriangleMatrix(CostsGraph &CG, double roundValue) {

        // create TriangleMatrix
        long size = CG.getSize();
        TriangleMatrix<long double> G = TriangleMatrix<long double>(size, 0);


        // fill TriangleMatrix with data from CostsGraph
        for (int i = 0; i < size; i++) {
            for (int k = 0; k < i; k++) {
                G.pos(i, k) = roundValue * CG.getEdge(i, k);
            }
        }
        return G;
    }


    // help-function to compute tuple (x,y) and to add it to our tuple list
    // more precise we save not only x and y, but x rounded up and down and analog y
    inline static void
    addTupel(double x, double y, vector<Napsack::Tupel> &B, long long &X, long long &Y, long double forbidden) {
        long long x_up = (long long) ceil(x);
        long long x_down = (long long) floor(x);
        long long y_up = (long long) ceil(y);
        long long y_down = (long long) floor(y);
        Napsack::Tupel value = {x_up, x_down, y_up, y_down};
        B.push_back(value);
        // increase matrix borders if value is not forbidden
        X += (x < forbidden / 10) ? 0 : (long long) ceil(fabs(x));
        Y += (y < forbidden / 10) ? 0 : (long long) ceil(fabs(y));
    }


    // function to compute a list of (x,y) tuples for edge uv
    // beyond this delta_u and deleta_v as written in the thesis is calculated
    // in addition the approximation of the critical clique rule is also calculated
    inline static vector<Napsack::Tupel>
    computeValues(int u, int v, TriangleMatrix<long double> &G, long long &X, long long &Y, long double &delta_u,
                  long double &delta_v, long double &rel_differenz, long double forbidden, long double edge_value) {

        // help variables
        long double x;
        long double y;
        delta_u = 0;
        delta_v = 0;
        X = 0LL;
        Y = 0LL;

        // create empty vector of tuples
        vector<Napsack::Tupel> B = vector<Napsack::Tupel>(0);
        rel_differenz = 0;

        // swap to guarantee that v > u
        if (v < u) {
            std::swap(u, v);
        }

        // neccessary for unweighted instances since real zero weight edges may occur
        // its < 0 since this increases the running times a little (really only little)
        long double zero = -0.00000001;

        // iterate over all vertices which are not u,v in the next 3 loops
        // first always get edges to neighbor i
        // then increase/decrease delta_u and delta_v or add tuple of neighbors of i to tuple list
        for (int i = 0; i < u; i++) {
            x = G.dirpos(u, i);
            y = G.dirpos(v, i);

            if (x > zero) {
                if (y > zero) {
                    addTupel(x, y, B, X, Y, forbidden);
                    if (rel_differenz != -forbidden) {
                        rel_differenz += fabs(x - y);
                    }
                } else {
                    delta_u += x;
                    delta_v -= y;
                }
            } else {
                if (y > zero) {
                    delta_u -= x;
                    delta_v += y;
                } else {
                    if (x <= forbidden) {
                        if (y > forbidden) {
                            addTupel(forbidden, y, B, X, Y, forbidden);
                            rel_differenz = -forbidden;
                        }
                    } else if (y <= forbidden) {
                        addTupel(x, forbidden, B, X, Y, forbidden);
                        rel_differenz += -forbidden;
                    } else {
                        addTupel(x, y, B, X, Y, forbidden);
                        if (rel_differenz != -forbidden) {
                            rel_differenz += fabs(x - y);
                        }
                    }
                }
            }
        }
        for (int i = u + 1; i < v; i++) {
            x = G.dirpos(i, u);
            y = G.dirpos(v, i);

            if (x > zero) {
                if (y > zero) {
                    addTupel(x, y, B, X, Y, forbidden);
                    if (rel_differenz != -forbidden) {
                        rel_differenz += fabs(x - y);
                    }
                } else {
                    delta_u += x;
                    delta_v -= y;
                }
            } else {
                if (y > zero) {
                    delta_u -= x;
                    delta_v += y;
                } else {
                    if (x <= forbidden) {
                        if (y <= forbidden) {
                            if (rel_differenz != -forbidden) {
                                rel_differenz += fabs(x - y);
                            }
                        } else {
                            addTupel(forbidden, y, B, X, Y, forbidden);
                            rel_differenz = -forbidden;
                        }
                    } else if (y <= forbidden) {
                        addTupel(x, forbidden, B, X, Y, forbidden);
                        rel_differenz += -forbidden;
                    } else {
                        addTupel(x, y, B, X, Y, forbidden);
                    }
                }
            }
        }
        for (int i = v + 1; i < G.size(); i++) {
            x = G.dirpos(i, u);
            y = G.dirpos(i, v);

            if (x > zero) {
                if (y > zero) {
                    addTupel(x, y, B, X, Y, forbidden);
                    if (rel_differenz != -forbidden) {
                        rel_differenz += fabs(x - y);
                    }
                } else {
                    delta_u += x;
                    delta_v -= y;
                }
            } else {
                if (y > zero) {
                    delta_u -= x;
                    delta_v += y;
                } else {
                    if (x <= forbidden) {
                        if (y <= forbidden) {
                            if (rel_differenz != -forbidden) {
                                rel_differenz += fabs(x - y);
                            }
                        } else {
                            addTupel(forbidden, y, B, X, Y, forbidden);
                            rel_differenz = -forbidden;
                        }
                    } else if (y <= forbidden) {
                        addTupel(x, forbidden, B, X, Y, forbidden);
                        rel_differenz += -forbidden;
                    } else {
                        addTupel(x, y, B, X, Y, forbidden);
                    }
                }
            }
        }

        return B;
    }

};

#endif
