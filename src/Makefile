ClusterEditing: array.o lowerbound.o graphexception.o costsgraph.o vertexlistsexception.o vertexlists.o  graphset.o probleminstanceexception.o edgereduction.o mergereduction.h bnmanager.o weightedprobleminstance.o searchtreeweighted.o  mincut.o almostclique.o criticalclique.o napsack.o trianglematrix.o cc-kernel.o

searchtreeweighted.o: searchtreeweighted.h searchtreeweighted.cpp graphset.o weightedprobleminstance.o
	g++ -c -O2 -I. searchtreeweighted.cpp

weightedprobleminstance.o: weightedprobleminstance.h weightedprobleminstance.cpp edgereduction.o  probleminstanceexception.o mergereduction.o bnmanager.o costsgraph.o array.o trianglematrix.o lowerbound.o mincut.o almostclique.o criticalclique.o cc-kernel.o
	g++ -c -O2 -I. weightedprobleminstance.cpp

cc-kernel.o: cc-kernel.h graphset.o costsgraph.o graphexception.o
	g++ -c -O2 -I. cc-kernel.h

criticalclique.o: criticalclique.h costsgraph.o trianglematrix.h Exceptions/graphexception.h graphset.h napsack.o
	g++ -c -O2 -I. criticalclique.h

napsack.o: napsack.h
	g++ -c -O2 -I. napsack.h

almostclique.o: almostclique.h array.h costsgraph.o
	g++ -c -O2 -I. almostclique.h

mincut.o: mincut.h array.h costsgraph.o
	g++ -c -O2 -I. mincut.h

bnmanager.o: bnmanager.h costsgraph.o edgereduction.o
	g++ -c -O2 -I. bnmanager.h

mergereduction.o: mergereduction.h costsgraph.o edgereduction.o
	g++ -c -O2 -I. mergereduction.h

edgereduction.o: edgereduction.cpp edgereduction.h costsgraph.o probleminstanceexception.o array.h trianglematrix.h mincut.h lowerbound.h
	g++ -c -O2 -I. edgereduction.cpp

lowerbound.o: costsgraph.o array.o lowerbound.h edgereduction.o
	g++ -c -O2 -I. lowerbound.h

probleminstanceexception.o: Exceptions/probleminstanceexception.cpp Exceptions/probleminstanceexception.h
	g++ -c -O2 -I. probleminstanceexception.cpp

graphset.o: graphset.cpp graphset.h vertexlists.o
	g++ -c -O2 -I. graphset.cpp

#edgelists.o: edgelists.cpp edgelists.h edgelistsexception.o array.h trianglematrix.h
#	g++ -c -O2 -I. edgelists.cpp
	
#edgelistsexception.o: edgelistsexception.cpp edgelistsexception.h
#	g++ -c -O2 -I. edgelistsexception.cpp
		
costsgraph.o: costsgraph.cpp costsgraph.h vertexlists.o graphexception.o array.h trianglematrix.h
	g++ -c -O2 -I. costsgraph.cpp

graphexception.o: Exceptions/graphexception.cpp Exceptions/graphexception.h
	g++ -c -O2 -I. graphexception.cpp

vertexlists.o: vertexlists.cpp vertexlists.h vertexlistsexception.o
	g++ -c -O2 -I. vertexlists.cpp
	
vertexlistsexception.o: Exceptions/vertexlistsexception.cpp Exceptions/vertexlistsexception.h
	g++ -c -O2 -I. vertexlistsexception.cpp

trianglematrix.o: trianglematrix.h array.h
	g++ -c -O2 -I. trianglematrix.h

array.o: array.h
	g++ -c -O2 -I. array.h