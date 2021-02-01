
#include <algorithm>
#include <iostream>
#include <numeric>
#include <sstream>
#include <cassert>

using namespace std;

// read graph into matrix
vector<vector<int>> readGraph() {

  std::istringstream sstream;
  auto getline = [&]() {
    std::string line;
    do {
      std::getline(cin, line);
    } while ( line[0] == 'c' );
    sstream = std::istringstream(line);
  };


  getline();
  std::string skip;
  int n, m;
  sstream >> skip >> skip >> n >> m;

  vector<vector<int>> res(n, vector<int>(n, -1));

  for ( int i = 0; i < m; ++i ) {
    getline();

    int u, v;
    sstream >> u >> v;
    --u; --v;
    res[u][v] = 1;
    res[v][u] = 1;
  }

  return res;
}


int main(int argc, char* argv[]) {
  string file = "../../../instances/exact/exact005.gr";

  if(argc>1)
    file = argv[1];
  freopen(file.c_str(), "r", stdin);
  cout << "reading file " << file << endl;
  auto edges = readGraph();
  cout << "n = " << size(edges) << endl;

  int n = size(edges);

  vector<vector<int>> clust{
    {1,6,7,8,9,10,13,14,15,19},
    {2,16,17,18,20},
    {3,4,5},
    {11,12},
  };

  vector solution(n,0);
  for(int c=0; c<size(clust); ++c)
    for(auto v : clust[c])
      solution[v-1] = c;


  int detail = 19;
  int plus = 0, neg = 0;

  int cost = 0;
  for(int i=0; i<n; ++i) 
    for(int j=i+1; j<n; ++j) {
      if(solution[i]!=solution[j]) cost += max(0, edges[i][j]);
      else cost += max(0,-edges[i][j]);
      if(i==detail-1||j==detail-1) {
        if(solution[i]!=solution[j] && edges[i][j]>0) cout << "+ " << i+1 << ' ' << j+1 << endl, plus++;
        if(solution[i]==solution[j] && edges[i][j]<0) cout << "- " << i+1 << ' ' << j+1 << endl, neg++;

      }
    }

  cout << "pos: " << plus<< endl;
  cout << "neg: " << neg << endl;

  cout << cost << endl;

  return 0;
}
