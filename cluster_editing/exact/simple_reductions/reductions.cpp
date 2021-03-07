#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int n = 0;

bool heavy_non_edge_rule() {


  return false;
}

int** makeAdjacencyMatrix() {
  string line;
  ifstream inputfile("graphs/exact001.gr");
  int** Adj = 0;

  if (inputfile.is_open()) { 
    getline(inputfile, line);

    string edges = line.substr(line.find_last_of(" ") + 1);
    string rest = line.substr(0,line.find_last_of(" "));
    string vertices = rest.substr(rest.find_last_of(" ") + 1);

    n = stoi(vertices);
    int e = stoi(edges);

    cout << n << e << endl;

    // int Adj[n+1][n+1];
    Adj = new int*[n];

    for (int i = 0; i < n + 1; i++) {
      Adj[i] = new int[n];
      for (int j = 0; j < n + 1; j++) {
        Adj[i][j] = 0;
      }
    }

    while (getline(inputfile, line)) {
      int v1 = stoi(line.substr(0, line.find(" "))) - 1;
      int v2 = stoi(line.substr(line.find(" ") + 1)) - 1;
      cout << v1 << " " << v2 << endl;
      Adj[v1][v2] = 1;
      Adj[v2][v1] = 1;

    }
    inputfile.close();
  }
  return Adj;
}


int main(int argc,  char **argv) {
  std::cout << "Implementing data reduction rules..." << std::endl;
  int** Adj = makeAdjacencyMatrix();
  
  // Traverse the Adj[][] 
  for (int i = 1; i < n + 1; i++) { 
      for (int j = 1; j < n + 1; j++) { 

          // Print the value at Adj[i][j] 
          printf("%d ", Adj[i][j]); 
      } 
      printf("\n"); 
  } 
  return 0;
}
