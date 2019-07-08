/*+++++++++++++++++ class WeightedProblemInstance +++++++++++++++++++ */
#include <weightedprobleminstance.h>

/* ##################### Constructors ###################### */

// constructor for as many edge lists as new the size of the parameter
// in egde list i are all elements inserted with costs i <= c < (i+1)
WeightedProblemInstance::WeightedProblemInstance(CostsGraph graph, double parameter, bool para_independent) : permanent(
        graph.permanent), forbidden(graph.forbidden), _parameter(parameter), _start_parameter(parameter),
                                                                                                              _changed_edges(
                                                                                                                      0),
                                                                                                              _graph(graph),
                                                                                                              _edge_reduction(),
                                                                                                              _pid(para_independent) {
    if (!_pid) {
        _edge_reduction.init(&_graph, parameter, true);
    }
}


// copy constructor
WeightedProblemInstance::WeightedProblemInstance(const WeightedProblemInstance &PI) : _graph(PI._graph),
                                                                                      permanent(PI.permanent),
                                                                                      forbidden(PI.forbidden),
                                                                                      _parameter(PI._parameter),
                                                                                      _start_parameter(
                                                                                              PI._start_parameter),
                                                                                      _changed_edges(PI._changed_edges),
                                                                                      _edge_reduction(), _pid(PI._pid) {
    if (!_pid) {
        _edge_reduction = EdgeReduction(PI._edge_reduction, &_graph);
    }
}


// copy constructor which gets elements as parameters
WeightedProblemInstance::WeightedProblemInstance(CostsGraph &graph, double parameter, double start_parameter,
                                                 const EdgeReduction &edge_reduction, bool para_independent) : _graph(
        graph),
                                                                                                               permanent(
                                                                                                                       graph.permanent),
                                                                                                               forbidden(
                                                                                                                       graph.forbidden),
                                                                                                               _parameter(
                                                                                                                       parameter),
                                                                                                               _start_parameter(
                                                                                                                       start_parameter),
                                                                                                               _changed_edges(
                                                                                                                       0),
                                                                                                               _pid(para_independent) {
    if (!_pid) {
        _edge_reduction = EdgeReduction(edge_reduction, &_graph);
    }
}



/* ###################### private/help functions ########## */


/* ---------------------- general help functions ------------------ */


// set edge to new value (permanent or forbidden)
// in any case the min insert and delete costs will be updated and
// consequently the egdes will be copied updated to other edge lists
// also the parameter will be decreased and list with corresponding
// costs bigger then the new parameter will be concatenated to the
// last deletion and insertion list
inline void WeightedProblemInstance::setEdgeToValue(int i, int j, double value) {
    //std::cout << "(" << i << "," << j << ") to " << value <<  std::endl;

    // get old value and calculate list of a value equal to the current parameter
    double old_value = _graph.getEdge(i, j);
    double parameter_before = _parameter;

    // set egde in graph
    _graph.setEdge(i, j, value);

    // indpendent from merging option
    if (old_value <= 0 && value > 0) {
        // save change
        saveChangedEdge(i, j, old_value * -1, 'a');
        // decrease parameter
        _parameter -= old_value * -1;
    }

    // case 2: edge was set and will be not_set
    if (old_value > 0 && value <= 0) {
        // save change
        saveChangedEdge(i, j, old_value, 'd');
        // decrease parameter
        _parameter -= old_value;
    }

    if (_parameter + 0.01 < 0.0) {
//std::cout << "1 " << _parameter << std::endl;
        throw ProblemInstanceException(" problem kernel reduction within given parameter range not possible ");
    }

    // update the underlying icf/icp costs matrices
    if (!_pid) {
        _edge_reduction.setEdgeToValue(i, j, old_value, _parameter, parameter_before);
    }

}


// merge all vertices given in the input clique vector
inline void WeightedProblemInstance::mergeVertices(CostsGraph::byte_vector_type &clique) {
    // now merge the whole clique
    for (int i = clique.size() - 1; i > 0; i--) {
        setEdgeToValue(clique[i], clique[i - 1], permanent);
        mergeVertices(clique[i], clique[i - 1]);
    }
}


// deletes a created clique by setting all surounding edges to forbidden and merging the clique
inline void WeightedProblemInstance::deleteClique(CostsGraph::byte_vector_type &clique) {
    // test whether i is in clique or not
    if (!clique.empty()) {

        // iterate over all clique members
        for (int k = 0; k < clique.size(); k++) {

            // now iterate over the whole graph
            int clNr = 0;

            for (int h = 0; h < clique[k]; h++) {
                if (clNr < clique.size() && h == clique[clNr]) {
                    // ignore permanent for now
                    clNr++;
                } else {
                    // any other edge to member vertex will be set to forbidden
                    setEdgeToValue(clique[k], h, forbidden);
                }
            }

            clNr++;

            for (int h = clique[k] + 1; h < _graph.getSize(); h++) {
                if (clNr < clique.size() && h == clique[clNr]) {
                    // ignore permanent for now
                    clNr++;
                } else {
                    // any other edge to member vertex will be set to forbidden
                    setEdgeToValue(clique[k], h, forbidden);
                }
            }
        }

        // now finally merge vertices to one vertex
        mergeVertices(clique);
    }
}


// check for vertex i whether it is in a completed clique
// if yes delete it from the rest of the graph and merge it
/* inline void WeightedProblemInstance::deleteClique(int i) {
    // get clique as vector
    CostsGraph::byte_vector_type clique = _graph.getClique(i);
    deleteClique(clique);
} */


// set all edges from the list for forbidden and update all underlying objects
inline void WeightedProblemInstance::setEdgesToForbidden(const CostsGraph::edge_list_type &forbidden_list) {
    // set edges to forbidden
    for (const auto &edge : forbidden_list) {
        setEdgeToValue(edge.i, edge.j, forbidden);
        //std::cout << "set(" << forbidden_list[i].i << "," << forbidden_list[i].j << ") forb -> new para=" << _parameter << std::endl;
    }
}


// sets list of edges to permanent, always the first index of an edge should be 
// greater then the second one !
inline void WeightedProblemInstance::setEdgesToPermanent(CostsGraph::edge_list_type permanent_list) {
    // set edges to permanent and merge, therefore sort them before!
    for (int i = 0; i < permanent_list.size(); i++) {
        if (permanent_list[i].i != permanent_list[i].j &&
            _graph.getEdge(permanent_list[i].i, permanent_list[i].j) != forbidden) {

            setEdgeToValue(permanent_list[i].i, permanent_list[i].j, permanent);
            mergeVertices(permanent_list[i].i, permanent_list[i].j);
            //std::cout << "set(" << permanent_list[i].i << "," << permanent_list[i].j << ") perm -> new para=" << _parameter << std::endl;

            // update vertex indices of other edges
            for (int k = i + 1; k < permanent_list.size(); k++) {
                if (permanent_list[k].i > permanent_list[i].i) {
                    permanent_list[k].i--;
                } else if (permanent_list[k].i == permanent_list[i].i) {
                    permanent_list[k].i = permanent_list[i].j;
                }
                if (permanent_list[k].j > permanent_list[i].i) {
                    permanent_list[k].j--;
                } else if (permanent_list[k].j == permanent_list[i].i) {
                    permanent_list[k].j = permanent_list[i].j;
                }
                if (permanent_list[k].i < permanent_list[k].j) {
                    std::swap<unsigned short>(permanent_list[k].i, permanent_list[k].j);
                }
            }
        }
    }
}


// try to apply simple reduction rules and do the following setting to permanent/forbidden also
inline int WeightedProblemInstance::simpleReduction() {
    int size = _graph.getSize();

    // create containers for edges to change
    std::vector<pair_type> permanent_list;
    std::vector<pair_type> forbidden_list;

    // check for edge to be set to permanent or forbidden
    for (int i = 0; i < size; i++) {
        for (int k = 0; k < i; k++) {
            double edge = _graph.getEdge(i, k);
            //std::cout << "(" << i << "," << k << ")" << std::endl;
            if (edge < 0 && edge != forbidden) {
                // light edge rule
                long double sum1 = 0.0;
                long double sum2 = 0.0;

                for (int h = 0; h < k; h++) {
                    sum1 += std::max(_graph.getEdge(i, h), 0.0);
                    sum2 += std::max(_graph.getEdge(k, h), 0.0);
                }
                for (int h = k + 1; h < i; h++) {
                    sum1 += std::max(_graph.getEdge(i, h), 0.0);
                    sum2 += std::max(_graph.getEdge(k, h), 0.0);
                }
                for (int h = i + 1; h < size; h++) {
                    sum1 += std::max(_graph.getEdge(i, h), 0.0);
                    sum2 += std::max(_graph.getEdge(k, h), 0.0);
                }
                //std::cout << "sum1=" << sum1 << " sum2=" << sum2 << std::endl;
                if (sum1 < -1 * edge || sum2 < -1 * edge) {
                    pair_type e = {i, k};
                    forbidden_list.push_back(e);
                }
            } else if (edge > 0 && edge != forbidden) {
                // heavy edge rule
                double sum1 = 0.0;
                double sum2 = 0.0;
                // heavy clique rule
                double sum3 = 0.0;
                for (int h = 0; h < k; h++) {
                    sum1 += std::abs(_graph.getEdge(i, h));
                    sum2 += std::abs(_graph.getEdge(k, h));
                    sum3 += std::max(_graph.getEdge(i, h), 0.0);
                    sum3 += std::max(_graph.getEdge(k, h), 0.0);
                }
                for (int h = k + 1; h < i; h++) {
                    sum1 += std::abs(_graph.getEdge(i, h));
                    sum2 += std::abs(_graph.getEdge(k, h));
                    sum3 += std::max(_graph.getEdge(i, h), 0.0);
                    sum3 += std::max(_graph.getEdge(k, h), 0.0);
                }
                for (int h = i + 1; h < size; h++) {
                    sum1 += std::abs(_graph.getEdge(i, h));
                    sum2 += std::abs(_graph.getEdge(k, h));
                    sum3 += std::max(_graph.getEdge(i, h), 0.0);
                    sum3 += std::max(_graph.getEdge(k, h), 0.0);
                }

                //std::cout << "sum3=" << sum3 << " sum2=" << sum2 << " sum1=" << sum1 << std::endl;
                if (sum3 < edge) {
                    pair_type e = {i, k};
                    permanent_list.push_back(e);
                }
                if (sum1 < edge || sum2 < edge) {
                    pair_type e = {i, k};
                    permanent_list.push_back(e);
                }
            }

        }
    }

    setEdgesToForbidden(forbidden_list);
    setEdgesToPermanent(permanent_list);

    return permanent_list.size() + forbidden_list.size();
}


// use weighted critical clique reduction as in thesis
// therefore find first a multiplier and then call CriticalClique object
// max_edge indicates the size of the greates edge should be this large after multiplication
inline int WeightedProblemInstance::criticalCliqueReduction(int max_edge, bool check_for_int) {
    // define multiplier
    double multiplier;

    // check whether it is an integer instance or not
    bool integer_instance;

    if (check_for_int) {
        // cut out first digit after the dot
        int zero = static_cast<int>(
                (fabs(_graph.getEdge(1, 0)) * 1 - static_cast<int>(fabs(_graph.getEdge(1, 0)) * 1)) * 10);
        integer_instance = zero == 0;

        // if this zero check all other first digits after the decimal point
        if (zero == 0) {
            // find edge with maximum weight
            for (int i = 0; i < _graph.getSize(); i++) {
                for (int k = 0; k < i; k++) {
                    if (_graph.getEdge(i, k) != _graph.forbidden) {
                        if (zero != static_cast<int>(
                                (fabs(_graph.getEdge(i, k)) * 1 - static_cast<int>(fabs(_graph.getEdge(i, k)) * 1)) *
                                10)) {
                            integer_instance = false;
                        }
                    }
                }
            }
        }
    } else {
        integer_instance = false;
    }

    // if it is an integer instance, calculate maximum weighted instance, and therewith the multiplier
    if (!integer_instance) {
        // find edge with maximum weight
        double max = 0;
        for (int i = 0; i < _graph.getSize(); i++) {
            for (int k = 0; k < i; k++) {
                double x = _graph.getEdge(i, k);
                max = (x != _graph.forbidden && max < std::abs(x)) ? std::abs(x) : max;
            }
        }

        // calculate multiplier
        multiplier = max_edge / max;
    } else {
        // otherwise its one
        multiplier = 1;
    }

    // start reduction
    CostsGraph::edge_list_type permanent_list = CriticalClique::reduce(_graph, multiplier);

    // set edges to permanent
    setEdgesToPermanent(permanent_list);

    // return number of changes
    int changes = permanent_list.size();
    return changes;
}


// use weighted critical clique approximation as suggested in the thesis since
// it is realitvley easy
inline int WeightedProblemInstance::criticalCliqueApproximation() {
    int size = _graph.getSize();

    // create vector for permanent edges
    std::vector<pair_type> permanent_list = std::vector<pair_type>(0);

    // check for every edges...
    for (int i = 0; i < size; i++) {
        for (int k = 0; k < i; k++) {

            if (_graph.getEdge(i, k) > 0) {
                // check for neighbors and increase delta_u, delta_v and so on...see thesis
                double delta_i = 0.0;
                double delta_k = 0.0;
                double neighbor_sum = 0.0;

                double edge1;
                double edge2;
                for (int h = 0; h < k; h++) {
                    edge1 = _graph.getEdge(i, h);
                    edge2 = _graph.getEdge(k, h);
                    if (edge1 > 0) {
                        if (edge2 <= 0) {
                            delta_i += edge1;
                            delta_k -= edge2;
                        } else {
                            neighbor_sum += std::abs(edge1 - edge2);
                        }
                    } else if (edge2 > 0) {
                        delta_k += edge2;
                        delta_i -= edge1;
                    } else {
                        neighbor_sum += std::abs(edge1 - edge2);
                    }
                }
                for (int h = k + 1; h < i; h++) {
                    edge1 = _graph.getEdge(i, h);
                    edge2 = _graph.getEdge(k, h);
                    if (edge1 > 0) {
                        if (edge2 <= 0) {
                            delta_i += edge1;
                            delta_k -= edge2;
                        } else {
                            neighbor_sum += std::abs(edge1 - edge2);
                        }
                    } else if (edge2 > 0) {
                        delta_k += edge2;
                        delta_i -= edge1;
                    } else {
                        neighbor_sum += std::abs(edge1 - edge2);
                    }
                }
                for (int h = i + 1; h < size; h++) {
                    edge1 = _graph.getEdge(i, h);
                    edge2 = _graph.getEdge(k, h);
                    if (edge1 > 0) {
                        if (edge2 <= 0) {
                            delta_i += edge1;
                            delta_k -= edge2;
                        } else {
                            neighbor_sum += std::abs(edge1 - edge2);
                        }
                    } else if (edge2 > 0) {
                        delta_k += edge2;
                        delta_i -= edge1;
                    } else {
                        neighbor_sum += std::abs(edge1 - edge2);
                    }
                }

                // check if rule applies, and edge can be set to permanent
                if (_graph.getEdge(i, k) > 0.5 * (neighbor_sum + delta_i + delta_k)) {

                    /*std::cout << "CC edge (" << i << "," << k << ") =" << _graph.getEdge(i,k) << std::endl;
                    std::cout << "neighbor=" << neighbor_sum << std::endl;
                    std::cout << "delta i=" << delta_i << std::endl;
                    std::cout << "delta j=" << delta_k << std::endl;*/
                    pair_type e = {i, k};
                    permanent_list.push_back(e);
                }
            }
        }
    }

    // set all edges to permanent
    setEdgesToPermanent(permanent_list);

    return permanent_list.size();
}


// use almost clique rule to detect cliques
inline int WeightedProblemInstance::detectCliques() {
    // empty vector for possible cliques
    std::vector<CostsGraph::byte_vector_type> cliques;

    int reduced = 0;

    // apply almost clique rule and get set of cliques
    AlmostClique::get(_graph, cliques, _parameter, false);

    // for every returned clique do...
    for (int i = 0; i < cliques.size(); i++) {

        // increase counter by number of edges in clique
        reduced += cliques[i].size() * (cliques[i].size() - 1) / 2;

        // merge clique
        mergeVertices(cliques[i]);

        // update vertices in other cliques
        for (int k = i + 1; k < cliques.size(); k++) {
            // update clique k
            for (int h = 0; h < cliques[k].size(); h++) {
                // update vertex h of clique k
                int red = 0;
                for (int l = 1; l < cliques[i].size(); l++) {
                    if (cliques[k][h] > cliques[i][l]) {
                        red++;
                    } else if (cliques[k][h] == cliques[i][l]) {
                        cliques[k][h] = cliques[i][0];
                    }
                }
                if (cliques[k][h] != cliques[i][0]) {
                    cliques[k][h] -= red;
                }
            }
            // iterate through clique to delete equal vertices
            for (int h = 0; h < cliques[k].size(); h++) {
                for (auto l = cliques[k].end() - 1;
                     l != cliques[k].begin() + h; l--) {
                    if (*l == cliques[k][h]) {
                        cliques[k].erase(l);
                    }
                }
            }
            // sort clique again so that first vertex is the greates
            // to obtain no conflict with the setEgdesToPermanent function
            sort(cliques[k].begin(), cliques[k].end());
        }
    }


    return reduced;
}


// will fullfill the complete merge operation of the vertices i and j
// it will update the underlying graph with the new edge values,
// the min insert and delete costs matrices and
// might decrease the parameter, since the merge operation can trigger
// new immediate changes
inline void WeightedProblemInstance::mergeVertices(int i, int j) {
    // create a vector just for the new edge value from the now merged ij and all other nodes
    double_array_type new_costs;

    // since the matrices are triangular make i the bigger index
    if (i < j) {
        std::swap(i, j);
    }

    // get new costs adjecencies for new vertex and save modification costs
    double to_pay = MergeReduction::mergeVertices(i, j, _graph, forbidden, new_costs);

    // save changend edges
    saveChangedEdge(i, j, to_pay, 'm');

    // decrease parameter
    _parameter -= to_pay;

    if (_parameter + 0.01 < 0.0) {
//std::cout << "2 " << _parameter << std::endl;
        throw ProblemInstanceException(" problem kernel reduction within given parameter range not possible ");
    }

    // update min insert and delete costs
    if (!_pid) {
        _edge_reduction.mergeVertices(i, j, _parameter, new_costs);
    }

    // merge vertices and add new costs for the new vertex as calculated above
    _graph.mergeVertices(i, j, new_costs);

}


/* ################### main/reduce function ################### */

// set edge in graph to certain value and considers reduction of parameter and
// rules of the problem kernel
void WeightedProblemInstance::setEdge(int i, int j, double value) {

    if (value != permanent && value != forbidden) {
//std::cout << "3 " << value << std::endl;
        throw ProblemInstanceException(" value is not forbidden or permanent ");
    }

    if (value != _graph.getEdge(i, j)) {
        if (_graph.getEdge(i, j) == forbidden || _graph.getEdge(i, j) == permanent) {
//std::cout << "4 " << _graph.getEdge(i,j) << std::endl;
            throw ProblemInstanceException(" edge is already forbidden or permanent ");
        }

        // set edge to its value and update underlying objects
        setEdgeToValue(i, j, value);

        if (value == permanent) {
            // merge vertices and update all underlying objects
            mergeVertices(i, j);
        }

        if (_parameter + 0.01 < 0.0) {
//std::cout << "5 " << _parameter << std::endl;
            throw ProblemInstanceException(" problem kernel reduction within given parameter range not possible ");
        }
    }
}


// save changed edge in list
inline void WeightedProblemInstance::saveChangedEdge(int i, int j, double costs, char operation) {
    edge_type new_edge = {_graph.getVertexIndices(i), _graph.getVertexIndices(j), costs, operation};
    _changed_edges.push_back(new_edge);
}


// reduce parameter if possible and return the number of edges set to permanent or forbidden
int WeightedProblemInstance::reduce() {
    // counter for number of edges set to permanent or forbidden
    int count = 0;


    if (!_pid) {
        CostsGraph::edge_list_type permanent_list;
        CostsGraph::edge_list_type forbidden_list;

        while (_edge_reduction.getReduceableEdges(permanent_list, forbidden_list, 0.0)) {
            count += permanent_list.size() + forbidden_list.size();
            setEdgesToForbidden(forbidden_list);
            setEdgesToPermanent(permanent_list);
        }
    }

    return count;
}


// reduce parameter if possible and return the number of edges set to permanent or forbidden
// this function uses all reduction rules, but not in their full power
int WeightedProblemInstance::strongReduce() {
    // counter for number of edges set to permanent or forbidden
    int count = 0;
    int old_count = -30000;

    while (count - old_count > _graph.getSize()) {

        while (count - old_count > 0) {
            old_count = count;
            //std::cout << "ccaprox" << std::endl;
            int cr = criticalCliqueApproximation();
            while (cr > 0) {
                count += cr;
                cr = criticalCliqueApproximation();
            }
            //std::cout << "sr" << std::endl;
            int sr = simpleReduction();
            count += sr;
        }

        old_count = count;

        /*std::cout << _graph << std::endl;
        std::cout << "pid" << std::endl;*/
        if (!_pid) {
            CostsGraph::edge_list_type permanent_list;
            CostsGraph::edge_list_type forbidden_list;
            while (_edge_reduction.getReduceableEdges(permanent_list, forbidden_list, 0.05)) {
                count += permanent_list.size() + forbidden_list.size();
                setEdgesToForbidden(forbidden_list);
                setEdgesToPermanent(permanent_list);
            }
        }

        //std::cout << "almost" << std::endl;
        int dc = detectCliques();
        count += dc;

        //std::cout << "ccred" << std::endl;
        count += criticalCliqueReduction(10, true);

    }

    return count;
}


// reduce parameter if possible and return the number of edges set to permanent or forbidden
// this function uses all reduction rules in their full power and
// reduces as long as there is nothing to reduce anymore
int WeightedProblemInstance::maxReduce() {
    // counter for number of edges set to permanent or forbidden
    int count = 0;

    //std::cout << _graph << std::endl;
    // if the instance is ment to be unweighted use guos 4k kernalization first
    //std::cout << "cc merging.." << std::endl;
    ccReduce();
    int old_count = -30000;

    while (count - old_count > 0) {

        while (count - old_count > 0) {
            old_count = count;

            /*std::cout << "ccaprox" << std::endl;
            std::cout << _graph << std::endl;
            std::cout << _parameter << std::endl;*/
            int cr = criticalCliqueApproximation();
            while (cr > 0) {
                count += cr;
                /*std::cout << "ccaprox" << std::endl;
                std::cout << _graph << std::endl;
                std::cout << _parameter << std::endl;*/
                cr = criticalCliqueApproximation();
            }
            /*std::cout << "sr" << std::endl;
            std::cout << _graph << std::endl;
            std::cout << _parameter << std::endl;*/
            int sr = simpleReduction();
            count += sr;
        }

        old_count = count;
        /*std::cout << "pid" << std::endl;
        std::cout << _graph << std::endl;
        std::cout << _parameter << std::endl;*/
        if (!_pid) {
            CostsGraph::edge_list_type permanent_list;
            CostsGraph::edge_list_type forbidden_list;

            while (_edge_reduction.getReduceableEdges(permanent_list, forbidden_list, 1)) {
                count += permanent_list.size() + forbidden_list.size();
                setEdgesToForbidden(forbidden_list);
                setEdgesToPermanent(permanent_list);
            }
        }
        //std::cout << "almost" << std::endl;
        int dc = detectCliques();
        count += dc;

        //std::cout << "ccred" << std::endl;
        //int cc = criticalCliqueReduction(1000, false);
        int cc = criticalCliqueReduction(1000, true);
        count += cc;

    }

    return count;
}

// start the merging of critical cliques
int WeightedProblemInstance::ccReduce() {
    // counter for number of edges set to permanent or forbidden
    int count = 0;

    bool valid = true;

    // test whether graph is really unweighted
    for (int i = 0; i < _graph.getSize(); i++) {
        for (int k = 0; k < i; k++) {
            valid = fabs(_graph.getEdge(i, k)) == 1.0;
        }
    }

    // if graph is really unweighted
    if (valid) {
        // create empty lists for permanent & forbidden edges
        CostsGraph::edge_list_type permanent_list;

        // apply critical clique kernalization
        CCKernel::mergeCriticalCliques(_graph, permanent_list);


        // count edge changes
        count += permanent_list.size();

        // correct edges such that first index is always greater than second
        for (auto &edge : permanent_list) {
            if (edge.i < edge.j) {
                std::swap(edge.i, edge.j);
            }
        }

        // set edges to permanent
        setEdgesToPermanent(permanent_list);
    }

    return count;
}

// do the actual 4k kernelization, not just critical cliques
int WeightedProblemInstance::ccKernelization() {
    // counter for number of edges set to permanent or forbidden
    int count = 0;

    bool valid = true;

    // test whether graph is really unweighted
    for (int i = 0; i < _graph.getSize(); i++) {
        for (int k = 0; k < i; k++) {
            valid = fabs(_graph.getEdge(i, k)) == 1.0;
        }
    }

    // if graph is really unweighted
    if (valid) {
        // create empty lists for permanent & forbidden edges
        CostsGraph::edge_list_type permanent_list;
        CostsGraph::edge_list_type forbidden_list;

        // apply critical clique kernelization
        CCKernel::makeCCKernel(_graph, permanent_list, forbidden_list);


        // count edge changes
        //std::cout << "edges to permanent=" << permanent_list.size() << std::endl;
        //std::cout << "edges to forbidden=" << forbidden_list.size() << std::endl;
        count += permanent_list.size() + forbidden_list.size();

        // set edges to forbidden
        setEdgesToForbidden(forbidden_list);

        // correct edges such that first index is always greater than second
        for (auto &edge : permanent_list) {
            if (edge.i < edge.j) {
                std::swap(edge.i, edge.j);
            }
        }

        // set edges to permanent
        setEdgesToPermanent(permanent_list);
    }

    return count;
}


// heuristic approach for a weighted problem instance
// this method simply gets edge with high icp and icf and inserts or deletes them
int WeightedProblemInstance::heuristicSolve() {
    int help_i;
    int help_j;
    int count = 0;

    //std::cout << _graph << std::endl;
    while (getEdgeForBranching(help_i, help_j) != -1) {
        CostsGraph::edge_list_type permanent_list;
        CostsGraph::edge_list_type forbidden_list;

        forbidden_list = _edge_reduction.getTopEdgesToDelete((int) (_graph.getSize() / 10 + 1));
        count += forbidden_list.size();
        for (int i = 0; i < forbidden_list.size(); i++) {
            //	std::cout << "f=(" << forbidden_list[i].i << "," << forbidden_list[i].j << ")" << std::endl;
        }

        setEdgesToForbidden(forbidden_list);

        permanent_list = _edge_reduction.getTopEdgesToInsert((int) (_graph.getSize() / 10 + 1));
        count += permanent_list.size();
        for (int i = 0; i < permanent_list.size(); i++) {
            //	std::cout << "p=(" << permanent_list[i].i << "," << permanent_list[i].j << ")" << std::endl;
        }

        setEdgesToPermanent(permanent_list);
        //std::cout << _graph.getSize() << std::endl;
    }
    //std::cout << "end:" << _graph << std::endl;
    return count; // Added to avoid return-void warning
}

double WeightedProblemInstance::getUpperBound() const {
    WeightedProblemInstance copy_instance = *(this->clone());
    copy_instance.setParameter(1E10);
    copy_instance.heuristicSolve();
    return (1E10 - copy_instance.getParameter());
}


// create deep copy of object
WeightedProblemInstance *WeightedProblemInstance::clone() const {
    auto newObjPtr = new WeightedProblemInstance(*this);
    return newObjPtr;
}


// return edge with minimal branching number
// return -1 if there exits none
double WeightedProblemInstance::getEdgeForBranching(int &i, int &j) const {
    if (!_pid) {
        //std::cout << "check bn" << std::endl;
        return BNManager::getMinBNEdge(i, j, &_graph, &_edge_reduction);
    } else {
        return -1;
    }
}


// return lower bound for this instance
double WeightedProblemInstance::getLowerBound() const {
    // depending on what kind of instance use of the lower bounds
    if (!_pid) {
        return LowerBound::get3(_graph, _edge_reduction.getMinInsertMatrix(), _edge_reduction.getMinDeleteMatrix());
    } else {
        return LowerBound::get2(_graph);
    }
}


// divide this object under the constraint given in the vertex matrix
// normally each row corresponds to one connected component
// always all vertices in one row will packed in one new object
WeightedProblemInstance::pi_list_type
WeightedProblemInstance::divideInstance(const CostsGraph::vertex_matrix_type &vertex_matrix) {
    // build new graphs
    CostsGraph::graph_list_type graph_list = _graph.getConnectedComponents(vertex_matrix);
    int cc_nr = vertex_matrix.size();

    // build reduction object lists
    EdgeReduction::edge_red_list_type edge_red_list;
    if (!_pid) {
        edge_red_list = _edge_reduction.divideInstance(vertex_matrix, graph_list);
    } else {
        edge_red_list = EdgeReduction::edge_red_list_type(cc_nr, NULL);
    }

    // create list of problem instances
    pi_list_type pi_list = pi_list_type(cc_nr);

    // create copies
    for (int i = 0; i < cc_nr; i++) {
        pi_list[i] = new WeightedProblemInstance(*graph_list[i], _parameter, _start_parameter, *edge_red_list[i], _pid);
        delete (graph_list[i]);
        delete (edge_red_list[i]);
    }

    return pi_list;

}


// set parameter manually to value
// will cause merging of list between old and new parameter
// and a reduction step will be called afterwards
void WeightedProblemInstance::setParameter(double new_parameter) {
    if (!_pid) {
        _edge_reduction.setParameter(new_parameter);
    }

    _parameter = new_parameter;

    reduce();
}

