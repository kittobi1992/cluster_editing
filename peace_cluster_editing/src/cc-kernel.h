/*+++++++++++++++++ class CCKernel +++++++++++++++++++ */

#include <fstream>
#include <iostream>
#include <sstream>

#include <vector>

#include <graphset.h>
#include <costsgraph.h>
#include <Exceptions/graphexception.h>

using std::string;
using std::vector;


class CCKernel {
public:
    CCKernel() = default;

    ~CCKernel() = default;

    struct info { // information about the vertices in the Critical Clique Graph
        int number; // number of the vertices in the Critical Clique
        vector<int> names; // "names" of the vertices in the Critical Clique
        info() : number(0), names(0) {}
    };


    inline static void mergeCriticalCliques(CostsGraph &CG, CostsGraph::edge_list_type &permanent) {
        vector<vector<bool>> G = createMatrix(CG);
        vector<info> CC = findCriticalClique(G, permanent);
    }

    inline static void
    makeCCKernel(CostsGraph &CG, CostsGraph::edge_list_type &permanent, CostsGraph::edge_list_type &forbidden) {
        vector<vector<bool>> G = createMatrix(CG);
        vector<info> CC = findCriticalClique(G, permanent);
        vector<vector<bool>> CCG = makeCritcalCliqueGraph(CC, G);
        solveCC(CCG, CC, G, permanent, forbidden);
    }

    inline static vector<int>
    findNeighbours(int u, const vector<vector<bool>> &G) { // do NOT contains u (if G is not reflexive)
        vector<int> neighbours;
        for (int i = 0; i < G.size(); i++) {
            if ((G[u][i])) {
                neighbours.push_back(i);
            }
        }
        return neighbours;
    }

    inline static vector<int> findNeighbours2(int u, const vector<int> &neighbours,
                                              const vector<vector<bool> > &G) { // do not contain u or the neighbours ifselfs
        vector<int> neighbours2;
        for (int i = 0; i < neighbours.size(); i++) {
            for (int j = 0; j < G.size(); j++) {
                if ((G[neighbours[i]][j]) && (j != u) && testIfInNeighbour(j, neighbours) &&
                    testIfInNeighbour(j, neighbours2)) {
                    neighbours2.push_back(j);
                }
            }
        }
        return neighbours2;
    }

    inline static vector<int> closedNeighborhood(int u, const vector<vector<bool>> &G) { // contains u !!!
        vector<int> neighbours;
        for (int i = 0; i < G.size(); i++) { // make a list of the neighbours of u
            if ((G[u][i]) || u == i) {
                neighbours.push_back(i);
            }
        }
        return neighbours;
    }

    // tests if a vertice j for neighbours2 is already in neigbours 
    inline static bool testIfInNeighbour(int j, const vector<int> &neighbours) {
        bool help = true;
        for (int help2 : neighbours) {
            if (help2 > j) {
                break;
            }
            if (j == help2) {
                help = false;
            }
        }
        return help;
    }

    inline static bool testifnew(int number_old, const CostsGraph::edge_list_type &list) {
        bool help1 = true;
        for (const auto& pair_new : list) {
            for (int l = 0; l < number_old; l++) {
                auto pair_org = list[l];
                if (((pair_org.i == pair_new.i) && (pair_org.j == pair_new.j)) ||
                    ((pair_org.i == pair_new.j) && (pair_org.j == pair_new.i))) {
                    help1 = false;
                    break;
                } else {
                    help1 = true;
                }
            }
            if (help1) {
                return true;
            }
        }
        return false;
    }


    inline static void solveCC(vector<vector<bool> > &CCG, vector<info> &CC, vector<vector<bool> > &G,
                               CostsGraph::edge_list_type &permanent, CostsGraph::edge_list_type &forbidden) {
        vector<int> neighbours;
        vector<int> neighbours2;
        bool keep_going = true;
        bool help2;
        bool help3 = true;
        while (keep_going) {
            keep_going = false;
            int i = 0;
            help2 = true;
            while (help2 && (!CCG.empty())) { // do for each vertice in the Criticcal Clique Graph....
                //while (help2){ // do for each vertice in the Criticcal Clique Graph....
                neighbours = findNeighbours(i, CCG);
                neighbours2 = findNeighbours2(i, neighbours, CCG);
                if (neighbours.empty()) {
                    deletion(i, CCG, CC);
                    i--;
                    keep_going = true;
                } else if (testRule2(neighbours, neighbours2, CC[i].number, CC)) {
                    int old_permanent = permanent.size();
                    int old_forbidden = forbidden.size();
                    makeRule2(i, neighbours, neighbours2, permanent, forbidden, G, CC, CCG);
                    if (testifnew(old_permanent, permanent) || testifnew(old_forbidden, forbidden)) {
                        keep_going = true;
                    }
                } else if (!neighbours2.empty()) {
                    for (int j = 0; j < neighbours.size(); j++) {
                        if (testRule3(i, neighbours, CC[i].number, CC[neighbours[j]].number, neighbours[j], CCG, CC,
                                      G)) {
                            makeRule3(i, neighbours[j], neighbours2, permanent, forbidden, G, CC, CCG);
                            int old_forbidden = forbidden.size();
                            if (testifnew(old_forbidden, forbidden)) {
                                keep_going = true;
                            }
                            break;
                        }
                    }
                }
                i++;
                if (i >= CCG.size()) {
                    help2 = false;
                }
            }
        }
    }

    inline static void deletion(int i, vector<vector<bool> > &CCG, vector<info> &CC) {
        for (auto& c : CCG) {
            c.erase(c.begin() + i);
        }
        CCG.erase(CCG.begin() + i);
        CC.erase(CC.begin() + i);
    }

    inline static bool
    testRule2(const vector<int> &neighbours, const vector<int> &neighbours2, int numberOfVerticesInCC, vector<info> &CC) {
        int numberOfNeighbours = 0;
        int numberOf2Neighbours = 0;

        for (int neighbour : neighbours) {
            numberOfNeighbours += CC[neighbour].number;
        }
        for (int neighbour : neighbours2) {
            numberOf2Neighbours += CC[neighbour].number;
        }
        return (numberOfVerticesInCC >= numberOfNeighbours + numberOf2Neighbours);
    }

    inline static bool
    testRule3(int i, const vector<int> &neighbours, int numberOfVerticesInCC, int numberOfVerticesInCC2, int j,
              const vector<vector<bool> > &CCG, const vector<info> &CC, const vector<vector<bool> > &G) {
        int numberOfNeighbours = 0;
        for (int neighbour : neighbours) {
            numberOfNeighbours += CC[neighbour].number;
        }
        if (numberOfVerticesInCC >= numberOfNeighbours) {
            int sum_1 = sum1(i, neighbours, j, CCG, CC, G);
            int sum_2 = sum2(neighbours, j, CC, G);
            return (numberOfVerticesInCC * numberOfVerticesInCC2 >= sum_1 + sum_2);
        } else {
            return false;
        }
    }

    inline static int sum1(int i, const vector<int> &neighbours, int j, const vector<vector<bool>> &CCG, const vector<info> &CC,
                           const vector<vector<bool>> &G) { // i= k and j= k'
        int sum = 0;
        for (int neighbour2 : findNeighbours2(i, neighbours, CCG)) {
            for (int name_l : CC[j].names) {
                for (int name_m : CC[neighbour2].names) {
                    if ((G[name_l][name_m]) && (name_l != name_m)) {
                        sum++;
                    }
                }
            }
        }
        return sum;
    }

    inline static int sum2(const vector<int> &neighbours, int j, const vector<info> &CC, const vector<vector<bool>> &G) {
        int sum = 0;
        for (const auto& n_k : neighbours) {
            if (j != n_k) {
                for (int name_l : CC[j].names) {
                    for (int name_m : CC[n_k].names) {
                        if ((!G[name_l][name_m]) && (name_l != name_m)) {
                            sum++;
                        }
                    }
                }
            }
        }
        return sum;
    }

    inline static void
    makeRule2(int &k, const vector<int> &neighbours, const vector<int> &neighbours2, CostsGraph::edge_list_type &permanent,
              CostsGraph::edge_list_type &forbidden, const vector<vector<bool> > &G, vector<info> &CC,
              vector<vector<bool> > &CCG) {
        // set all edges between k and its neighbours permanent and delete the neighbours
        for (int n_i : neighbours) {
            CostsGraph::pair_type pairs(CC[k].names[0], CC[n_i].names[0]);
            if (pairs.i != pairs.j) {
                permanent.push_back(pairs);
            }
        }

        for (int n_i : neighbours) {
            for (int n_j : neighbours2) {
                for (int pi : CC[n_i].names) {
                    for (int pj : CC[n_j].names) {
                        forbidden.push_back({pi, pj});
                    }
                }
            }
        }
        // delete of the vertex k itself and its neighbours
        deletion(k, CCG, CC);
        k--;
        for (int i = neighbours.size() - 1; i >= 0; i--) {
            if (k < neighbours[i]) {
                deletion(neighbours[i] - 1, CCG, CC);
            } else {
                deletion(neighbours[i], CCG, CC);
                k--;
            }
        }
    }

    inline static void makeRule3(int &k1, int &j, vector<int> &neighbours2, CostsGraph::edge_list_type &permanent,
                                 CostsGraph::edge_list_type &forbidden, const vector<vector<bool>> &G, vector<info> &CC,
                                 vector<vector<bool> > &CCG) {
        vector<int> neighbours;
        bool h = true;
        int i = 0;
        while (h) {
            for (int pi : CC[j].names) {
                for (int pj : CC[neighbours2[i]].names) {
                    forbidden.push_back({pi, pj});
                }
            }
            mergeInCCG(k1, neighbours2[i], j, CCG, CC, permanent);
            i++;
            neighbours = findNeighbours(k1, CCG);
            neighbours2 = findNeighbours2(k1, neighbours, CCG);
            if (i >= neighbours2.size()) {
                h = false;
            }
        }
    }

    inline static void mergeInCCG(int &k1, int &i, int &j, vector<vector<bool>> &CCG, vector<info> &CC,
                                  CostsGraph::edge_list_type &permanent) {
        vector<int> neighbours_a;
        vector<int> neighbours_b;
        bool h = true;
        // Kante zwischen k' und k'' wird gelöscht
        CCG[i][j] = false;
        CCG[j][i] = false;
        neighbours_a = closedNeighborhood(i, CCG);
        int k = 0;
        while (h) {
            neighbours_b = closedNeighborhood(neighbours_a[k], CCG);
            if ((neighbours_a == neighbours_b) &&
                (i != neighbours_a[k])) { // wenn die Nachbarschaft gleich ist, soll gemergt werden
                permanent.push_back({CC[i].names[0], CC[neighbours_a[k]].names[0]});
                // es wird der Knoten gelöscht der mit k' benachbart ist, nicht k'
                for (auto& ccg : CCG) {
                    ccg.erase(ccg.begin() + neighbours_a[k]);
                }
                CCG.erase(CCG.begin() + neighbours_a[k]);
                CC[i].names.insert(CC[i].names.end(), CC[neighbours_a[k]].names.begin(),
                                   CC[neighbours_a[k]].names.end());
                CC[i].number += CC[neighbours_a[k]].number;
                CC.erase(CC.begin() + neighbours_a[k]);
                if (neighbours_a[k] <= i) {
                    i--;
                }
                if (neighbours_a[k] <= k1) {
                    k1--;
                }
                if (neighbours_a[k] <= j) {
                    j--;
                }
            }
            k++;
            neighbours_a = closedNeighborhood(i, CCG);
            if (k >= neighbours_a.size()) {
                h = false;
            }
        }
        // k = 0;
        neighbours_a = closedNeighborhood(j, CCG);
        // unreachable
        /* while (h) {
            neighbours_b = closedNeighborhood(neighbours_a[k], CCG);
            if ((neighbours_a == neighbours_b) && (j != neighbours_a[k])) {
                permanent.push_back({CC[j].names[0], CC[neighbours_a[k]].names[0]});
                for (int l = 0; l < CCG.size(); l++) {
                    CCG[l].erase(CCG[l].begin() + neighbours_a[k]);
                }
                CCG.erase(CCG.begin() + neighbours_a[k]);
                CC[j].names.insert(CC[j].names.end(), CC[neighbours_a[k]].names.begin(),
                                   CC[neighbours_a[k]].names.end());
                CC[j].number += CC[neighbours_a[k]].number;
                CC.erase(CC.begin() + neighbours_a[k]);
                if (neighbours_a[k] <= k1) {
                    k1--;
                }
                if (neighbours_a[k] <= j) {
                    j--;
                }
            }
            k++;
            neighbours_a = closedNeighborhood(j, CCG);
            if (k >= neighbours_a.size()) {
                h = false;
            }
        } */
    }

    inline static vector<vector<bool> > createMatrix(const CostsGraph &CG) {
        // create Matrix
        long size = CG.getSize();
        vector<vector<bool> > G = vector<vector<bool>>(size, vector<bool>(size, false));
        // fill Matrix with data from CostsGraph
        for (int i = 0; i < size; i++) {
            G[i][i] = true;
            for (int k = 0; k < i; k++) {
                if (CG.getEdge(i, k) > 0) {
                    G[i][k] = true;
                    G[k][i] = true;
                }
            }
        }
        return G;
    }

    inline static vector<info> findCriticalClique(vector<vector<bool> > &G, CostsGraph::edge_list_type &permanent) {
        // to set all edges in a Critical Clique on permanent
        vector<bool> vertices = vector<bool>(G.size(), true); // false if the vertice is already in a Critical Clique
        vector<int> neighbours_a;
        vector<int> neighbours_b;
        bool keep_going = true;
        vector<info> CriticalCliquen;
        int u = 0;

        while (keep_going) { // how long there are Vertices, which are in no Critical Clique, do.....
            if (!vertices[u]) { // jump to the next vertex, if this vertex is already in an Critical Clique
                u++;
                continue;
            }
            vertices[u] = false;
            neighbours_a = closedNeighborhood(u, G);// make a list of the neighbours of u
            info CC;
            CC.names = vector<int>({u}); // u is the first member of the new Critical Clique
            CC.number = 1;
            for (int i = 1; i < G.size(); i++) { // run through all left vertices....
                if (!vertices[i]) { // junmp to the next vertice, if this vertice is already in an Critical Clique
                    continue;
                }
                neighbours_b = closedNeighborhood(i, G); // ...and make a list of its neigbours
                if (neighbours_a ==
                    neighbours_b) { // if the neigbourlists are the same, add this vertice to the Critical Clique
                    CC.names.push_back(i);
                    CC.number++;
                    //add a edge between the first member of the Critical Clique and the new member of it
                    permanent.push_back({u, i});
                }
            }
            for (const auto& name : CC.names) { // "removes" the vertices which are in a Critical Clique now
                vertices[name] = false;
            }
            for (int i = vertices.size() - 1; i >= 0; i--) { // tests if verices is "empty"
                keep_going = false;
                if (vertices[i]) {
                    keep_going = true;
                    break;
                }
            }
            CriticalCliquen.push_back(CC); // add the new Critical Clique to the List of Critical Cliques
            u++;
        }
        return CriticalCliquen;
    }

    inline static vector<vector<bool>> makeCritcalCliqueGraph(vector<info> &CCVertexList, vector<vector<bool>> &G) {
        auto n = CCVertexList.size();
        vector<vector<bool>> CCG(n, vector<bool>(n, false)); // Critical Clique Graph
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                bool b = true;
                for (int cc1 : CCVertexList[i].names) {
                    for (int cc2 : CCVertexList[j].names) {
                        if (!G[cc1][cc2]) {
                            b = false;
                            break;
                        }
                    }
                    if (!b) {
                        break;
                    }
                }
                if (b) {
                    CCG[i][j] = true;
                    CCG[j][i] = true;
                }
            }
        }
        return CCG;

    }

};

