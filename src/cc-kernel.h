/*+++++++++++++++++ class CCKernel +++++++++++++++++++ */

#include <fstream>
#include <iostream>
#include <sstream>

#include <vector>

#include <graphset.h>
#include <costsgraph.h>
#include <graphexception.h>

using std::string;
using std::vector;


class CCKernel {
public:
	CCKernel (){};
	~CCKernel(){};

	struct info { // informations about the vertices in the Critical Clique Graph
		int number; // number of the vertices in the Critical Clique
		vector <int> names; // "names" of the vertices in the Critical Clique
	};


	inline static void mergeCriticalCliques (CostsGraph &CG, CostsGraph::edge_list_type &permanent){
		vector <vector <bool> > G = createMatrix (CG);
		vector <info> CC = findCriticalClique(G,permanent);
	}

	inline static void makeCCKernel (CostsGraph &CG, CostsGraph::edge_list_type &permanent, CostsGraph::edge_list_type &forbidden){
		vector <vector <bool> > G = createMatrix (CG);
		vector <info> CC = findCriticalClique(G,permanent);
		vector <vector <bool> > CCG = makeCritcalCliqueGraph (CC,G);
		solveCC(CCG, CC, G, permanent, forbidden);
	}

	inline static vector <int> findNeighbours (int u, vector <vector <bool> > &G){ // do NOT contains u (if G is not reflexiv)
		vector <int> neighbours;
		for (int i=0; i < G.size(); i++){ 
			if ((G[u][i])){
				neighbours.push_back(i);
			}
		}	
		return neighbours;
	}

	inline static vector <int> findNeighbours2 (int u, vector <int> &neighbours, vector <vector <bool> > &G){ // do not contain u or the neighbours ifselfs
		vector <int> neighbours2;
		for (int i=0; i < neighbours.size(); i++){
			for (int j=0; j < G.size(); j++){
				if ((G[neighbours[i]][j]) && (j!=u) && testIfInNeighbour(j, neighbours) && testIfInNeighbour(j, neighbours2)){
					neighbours2.push_back(j);
				}
			}
		}
		return neighbours2;	
	}
	inline static vector <int> closedNeighborhood (int u, vector <vector <bool> > &G){ // contains u !!!
		vector <int> neighbours;
		for (int i=0; i < G.size(); i++){ // make a list of the neighbours of u
			if ((G[u][i]) || u==i){
				neighbours.push_back(i);
			}
		}
		return neighbours;
	}

    // tests if a vertice j for neighbours2 is already in neigbours 
	inline static bool testIfInNeighbour (int j,vector <int> &neighbours){
		bool help = true;
		int help2;
		for (int i=0; i < neighbours.size(); i++){
			help2 = neighbours[i];
			if (help2 > j){
				break;
			}
			if (j==help2){
				help = false;
			}
		}
		return help;
	}

	inline static bool testifnew (int number_old, CostsGraph::edge_list_type list){
		CostsGraph::pair_type pair_org;
		CostsGraph::pair_type pair_new;
		bool help1=true;
		bool help2=false;
		for (int k=number_old; k < list.size(); k++){
			pair_new=list[k];
			for (int l=0; l < number_old; l++){
				pair_org=list[l];
				if (((pair_org.i==pair_new.i) && (pair_org.j==pair_new.j)) || ((pair_org.i==pair_new.j) && (pair_org.j==pair_new.i))){
					help1=false;
					break;	
				}
				else{
					help1=true;
				}
			}
			help2= help2 || help1;
			if (help2){
				return true;
			}
		}
		return false;
	} 


	inline static void solveCC(vector <vector <bool> > &CCG, vector <info> &CC,vector <vector <bool> > &G, CostsGraph::edge_list_type &permanent, CostsGraph::edge_list_type &forbidden){
		vector <int> neighbours;
		vector <int> neighbours2;
		bool help=true;
		bool help2=true;
		bool help3=true;
		while (help){
			help=false;
			int i=0;
			help2=true;
			while (help2 && (CCG.size()>0)){ // do for each vertice in the Criticcal Clique Graph....
			//while (help2){ // do for each vertice in the Criticcal Clique Graph....
				neighbours = findNeighbours (i, CCG);
				neighbours2 = findNeighbours2 (i, neighbours, CCG);
				if (neighbours.size()==0){ 
					deletion (i, CCG, CC);
					i--;
					help=true;
				}
				else if (testRule2(neighbours,neighbours2, CC[i].number, CC)){
					int old_permanent = permanent.size();
					int old_forbidden = forbidden.size();
					makeRule2(i, neighbours, neighbours2, permanent, forbidden, G, CC, CCG);
					if (testifnew(old_permanent, permanent) || testifnew(old_forbidden, forbidden) ){
						help=true;
					}
					
				}
				else if (neighbours2.size() > 0){
					for (int j=0; j < neighbours.size(); j++){
						if (testRule3(i, neighbours, CC[i].number, CC[neighbours[j]].number, neighbours[j], CCG, CC, G)){
							makeRule3 (i, neighbours[j], neighbours2, permanent, forbidden, G, CC, CCG);
							int old_forbidden = forbidden.size();
							if (testifnew(old_forbidden, forbidden)){
								help=true;
							}
							break;
						}
					}
				}
				i++;
				if (i >= CCG.size()){
					help2=false;
				}
			}
		}
	}

	inline static void deletion (int i, vector <vector <bool> > &CCG, vector <info> &CC){
		for (int j=0; j < CCG.size(); j++){
			CCG[j].erase (CCG[j].begin()+i);
		}
		CCG.erase(CCG.begin()+i);
		CC.erase(CC.begin()+i);
	}

	inline static bool testRule2 (vector <int> neighbours, vector <int> neighbours2, int numberOfVerticesInCC, vector <info> &CC){
		int numberOfNeighbours=0;
		int numberOf2Neighbours=0;
		
		for (int i=0; i < neighbours.size(); i++){
			numberOfNeighbours+=CC[neighbours[i]].number;
		}
		for (int i=0; i < neighbours2.size(); i++){
			numberOf2Neighbours+=CC[neighbours2[i]].number;
		}
		if (numberOfVerticesInCC >= numberOfNeighbours + numberOf2Neighbours){
			return true;
		}
		else{
			return false;
		}
	}

	inline static bool testRule3 (int i, vector <int> &neighbours, int numberOfVerticesInCC, int numberOfVerticesInCC2, int j, vector <vector <bool> > &CCG,vector <info> &CC, vector <vector <bool> > &G){
		int numberOfNeighbours=0;
		for (int k=0; k < neighbours.size(); k++){
			numberOfNeighbours+=CC[neighbours[k]].number;
		}
		if (numberOfVerticesInCC >= numberOfNeighbours){
			int sum_1 = sum1(i, neighbours, j, CCG, CC, G);
			int sum_2 = sum2(neighbours, j, CC, G);
			if (numberOfVerticesInCC * numberOfVerticesInCC2 >= sum_1+sum_2){
				return true;
			}
			else{
				return false;
			}
		}
		else{
			return false;
		}
	}

	inline static int sum1(int i, vector <int> &neighbours, int j, vector <vector <bool> > &CCG, vector <info> &CC,vector <vector <bool> > &G){ // i= k and j= k'
		int help=0;
		vector <int> help1;
		vector <int> help2;
		vector <int> neighbours2 = findNeighbours2 (i, neighbours, CCG);
		for (int k=0; k < neighbours2.size(); k++){
			help1=CC[j].names;
			help2=CC[neighbours2[k]].names;
			for (int l=0; l < help1.size(); l++){
					for (int m=0; m < help2.size(); m++){
						if((G[help1[l]][help2[m]])&& (help1[l]!=help2[m])){
							help++;
						}
					}
			}
		}
		return help;
	}

	inline static int sum2(vector <int> &neighbours, int j, vector <info> &CC,vector <vector <bool> > &G){
		int help=0;
		vector <int> help1;
		vector <int> help2;
		for (int k=0; k < neighbours.size(); k++){
			if(j!=neighbours[k]){
				help1=CC[j].names;
				help2=CC[neighbours[k]].names;
				for (int l=0; l < help1.size(); l++){
					for (int m=0; m < help2.size(); m++){
						if((!G[help1[l]][help2[m]]) && (help1[l]!=help2[m])){
							help++;
						}
					}
				}
			}
		}
		return help;
	}

	inline static void makeRule2(int &k, vector <int> &neighbours, vector <int> &neighbours2, CostsGraph::edge_list_type &permanent, CostsGraph::edge_list_type &forbidden, vector <vector <bool> > &G, vector <info> &CC,  vector <vector <bool> > &CCG){
		CostsGraph::pair_type pairs;
		int help;
		vector <int> help1;
		vector <int> help2;
		help1 = CC[k].names;
		help=help1[0];
		// set all edges between k and its neighbours permanent and delete the neighbours
		for (int i=0; i < neighbours.size(); i++){ 
			help2=CC[neighbours[i]].names;
			pairs.i = help;
			pairs.j = help2[0];
			if(pairs.i!=pairs.j) {
				permanent.push_back(pairs);	
			}	
		}

		for (int i=0; i < neighbours.size(); i++){ 
			for (int j=0; j < neighbours2.size(); j++){
				help1=CC[neighbours[i]].names;
				help2=CC[neighbours2[j]].names;
				for (int k=0; k < help1.size(); k++){
					for (int l=0; l < help2.size(); l++){
						pairs.i = help1[k];
						pairs.j = help2[l];
						forbidden.push_back(pairs);
					}
				}
			}
		}
		// delete of the vertice k itself and its neighbours
		deletion(k, CCG, CC);
		k--;
		for (int i=neighbours.size()-1; i >= 0; i--){ 
			if (k < neighbours[i]){
				deletion(neighbours[i]-1, CCG, CC);   
			}
			else{
				deletion(neighbours[i], CCG, CC);
				k--;
			}
		}
	}

	inline static void makeRule3(int &k1, int &j, vector <int> &neighbours2, CostsGraph::edge_list_type &permanent, CostsGraph::edge_list_type &forbidden, vector <vector <bool> > &G, vector <info> &CC, vector <vector <bool> > &CCG){
		CostsGraph::pair_type pairs;
		vector <int> help1;
		vector <int> help2;
		vector <int> neighbours;
        	bool h=true;
        	int i=0;
		while (h){
			help1=CC[j].names;
			help2=CC[neighbours2[i]].names;
			for (int k=0; k < help1.size(); k++){
				for (int l=0; l < help2.size(); l++){
					pairs.i = help1[k];
					pairs.j = help2[l];
					forbidden.push_back(pairs);
				}
			}
			mergeInCCG (k1, neighbours2[i], j, CCG, CC, permanent);
			i++;
			neighbours = findNeighbours (k1, CCG);
			neighbours2 = findNeighbours2 (k1, neighbours, CCG);
			if (i >= neighbours2.size()){
				h=false;
			}
		}
	}

	inline static void mergeInCCG (int &k1,  int &i, int &j, vector <vector <bool> > &CCG, vector <info> &CC, CostsGraph::edge_list_type &permanent){
		CostsGraph::pair_type pairs;
		vector <int> neighbours_a;
		vector <int> neighbours_b;
		bool h=true;
		// Kante zwischen k' und k'' wird gelöscht
		CCG[i][j]=false;
		CCG[j][i]=false;
		neighbours_a = closedNeighborhood (i, CCG);
		int k=0;
		while(h){
			neighbours_b = closedNeighborhood (neighbours_a[k], CCG);
			if ((neighbours_a==neighbours_b) && (i!=neighbours_a[k])){ // wenn die Nachbarschaft gleich ist, soll gemergt werden
				pairs.i = CC[i].names[0];
				pairs.j = CC[neighbours_a[k]].names[0];
				permanent.push_back(pairs);
				// es wird der Knoten gelöscht der mit k' benachbar ist, nicht k'
				for (int l=0; l < CCG.size(); l++){
					 CCG[l].erase (CCG[l].begin()+neighbours_a[k]);
				}
				CCG.erase(CCG.begin()+neighbours_a[k]);
				CC[i].names.insert(CC[i].names.end(), CC[neighbours_a[k]].names.begin(), CC[neighbours_a[k]].names.end());
				CC[i].number+=CC[neighbours_a[k]].number;
				CC.erase(CC.begin()+neighbours_a[k]);
				if (neighbours_a[k] <= i){
					i--;
				}
				if (neighbours_a[k] <= k1){
					k1--;
				}
				if (neighbours_a[k] <= j){
					j--;
				}
			}
			k++;
			neighbours_a = closedNeighborhood (i, CCG);
			if (k >= neighbours_a.size()){
				h=false;
			}
		}
		k=0;
		neighbours_a = closedNeighborhood (j, CCG);
		while(h){
			neighbours_b = closedNeighborhood (neighbours_a[k], CCG);
			if ((neighbours_a==neighbours_b) && (j!=neighbours_a[k])){
				pairs.i = CC[j].names[0];
				pairs.j = CC[neighbours_a[k]].names[0];
				permanent.push_back(pairs);
				for (int l=0; l < CCG.size(); l++){
					CCG[l].erase (CCG[l].begin()+neighbours_a[k]);
				}
				CCG.erase(CCG.begin()+neighbours_a[k]);
				CC[j].names.insert(CC[j].names.end(), CC[neighbours_a[k]].names.begin(), CC[neighbours_a[k]].names.end());
				CC[j].number+=CC[neighbours_a[k]].number;
				CC.erase(CC.begin()+neighbours_a[k]);
				if (neighbours_a[k] <= k1){
					k1--;
				}
			        if (neighbours_a[k] <= j){
					j--;
		                }
			}
			k++;
			neighbours_a = closedNeighborhood (j, CCG);
			if (k >= neighbours_a.size()){
				h=false;
			}
		}
	}

	inline static vector <vector <bool> > createMatrix (CostsGraph &CG){
		// create Matrix		
		long size = CG.getSize();
		vector <vector <bool> > G = vector <vector <bool> > (size,vector <bool> (size,false));
		double help;
		// fill Matrix with data from CostsGraph
		G[0][0] = true;
		for (int i = 0; i < size; i++) {
			for (int k = 0; k < i; k++) {
				help= CG.getEdge(i,k);
				if (help > 0){
					G[i][k] = true; 
					G[k][i] = true; 
				}
				G[i][i] = true;
			}
		}
		return G; 
	} 

	inline static vector <info> findCriticalClique (vector <vector <bool> >  &G, CostsGraph::edge_list_type &permanent){
        // to set all edges in a Critical Clique on permanent
		CostsGraph::pair_type pairs; 
		vector <bool> vertices = vector <bool> (G.size(),true); // false if the vertice is already in a Critical Clique
		vector <int> neighbours_a;
		vector <int> neighbours_b;
		vector <int> help;
		int help2;
		bool help3=true;
		vector <info> CriticalCliquen;
		info CC;
		int u=0;

		while (help3){ // how long there are Vertices, which are in no Critical Clique, do.....
			if (!vertices[u]){ // junmp to the next vertice, if this vertice is already in an Critical Clique
				u++;
				continue;
			}
			vertices[u] = false;	
			neighbours_a = closedNeighborhood(u,G);// make a list of the neighbours of u
			CC.names.clear();
			help = CC.names;
			help.push_back(u);
			help2 = 1;
			CC.names = help; // u is the first member of the new Critical Clique
			CC.number = help2;
			for (int i = 1; i < G.size(); i++){ // run through all left vertices....
				if (!vertices[i]){ // junmp to the next vertice, if this vertice is already in an Critical Clique
					continue;
				} 
				neighbours_b = closedNeighborhood(i,G); // ...and make a list of its neigbours
				if (neighbours_a==neighbours_b){ // if the neigbourlists are the same, add this vertice to the Critical Clique
					help=CC.names;
					help.push_back(i);
					help2=CC.number+1;
					CC.names=help; 
					CC.number=help2;
					//add a edge between the first member of the Critical Clique and the new member of it 
					pairs.i = u;
					pairs.j = i;
					permanent.push_back(pairs);
				}
			}
			help=CC.names;
			for (int i=0; i < help.size(); i++){ // "removes" the vertices which are in a Critical Clique now
				vertices[help[i]]=false;
			}
			for (int i=vertices.size()-1; i >= 0; i--){ // tests if verices is "empty"
				help3=false;
				if (vertices[i]){
					help3=true;
					break;
				}
			}
			CriticalCliquen.push_back(CC); // add the new Critical Clique to the List of Critical Cliques
			u++;
		}
		return CriticalCliquen;
	}

	inline static vector <vector <bool> > makeCritcalCliqueGraph (vector <info> &CCVerticeList ,vector <vector <bool> > &G){
		vector <vector <bool> > CCG = vector <vector <bool> > (CCVerticeList.size(),vector <bool> (CCVerticeList.size(), false)); // Critical Clique Graph
		vector <int>  CC1;
		vector <int>  CC2;
		bool help;
		for (int i=0; i < CCVerticeList.size(); i++){
			CC1 = CCVerticeList[i].names;
			for (int j=i+1; j < CCVerticeList.size(); j++){
				CC2 = CCVerticeList[j].names;
				help=true;
				for (int k=0; k < CC1.size(); k++){
					for (int l=0; l < CC2.size(); l++){
						if(!G[CC1[k]][CC2[l]]){
							help=false;
							break;
						}
					}
					if (!help){
						break;
					}
				}
				if (help){
					CCG[i][j]=true;
					CCG[j][i]=true;
				}
			}
		}
		return CCG;
		
	}

};

