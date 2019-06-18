all: WeightedClusterEditing WeightedClusterReduction CriticalCliqueReduction GraphVizCreator GraphVizCreator_YED GraphProducer GraphProducer_iwpec


WeightedClusterEditing: weightedclusterediting.cpp stopwatch.h
	g++ ../src/graphexception.o ../src/costsgraph.o ../src/vertexlists.o ../src/vertexlistsexception.o ../src/graphset.o ../src/GraphParser/edgefileparser.o ../src/GraphParser/matrixparser.o ../src/CostsParser/blastparser.o ../src/CostsParser/doubleparser.o ../src/probleminstanceexception.o ../src/edgereduction.o  ../src/weightedprobleminstance.o ../src/searchtreeweighted.o -o weightedclusterediting -O2 -I. -I../src/ -I../src/GraphParser/ -I../src/CostsParser/ -I . weightedclusterediting.cpp

WeightedClusterReduction:  weightedclusterediting.cpp stopwatch.h
	g++ ../src/graphexception.o ../src/costsgraph.o ../src/vertexlists.o ../src/vertexlistsexception.o ../src/graphset.o ../src/GraphParser/edgefileparser.o ../src/GraphParser/matrixparser.o ../src/CostsParser/blastparser.o ../src/probleminstanceexception.o ../src/edgereduction.o  ../src/weightedprobleminstance.o ../src/searchtreeweighted.o -o weightedclusterreduction -O2 -I. -I../src/ -I../src/GraphParser/ -I../src/CostsParser/ -I . weightedclusterreduction.cpp

CriticalCliqueReduction: criticalcliquereduction.cpp stopwatch.h
	g++ ../src/graphexception.o ../src/costsgraph.o ../src/vertexlists.o ../src/vertexlistsexception.o ../src/graphset.o ../src/GraphParser/edgefileparser.o ../src/GraphParser/matrixparser.o ../src/CostsParser/blastparser.o ../src/probleminstanceexception.o ../src/edgereduction.o  ../src/weightedprobleminstance.o ../src/searchtreeweighted.o -o criticalcliquereduction -O2 -I. -I../src/ -I../src/GraphParser/ -I../src/CostsParser/ -I . criticalcliquereduction.cpp

GraphVizCreator: graphvizcreator.cpp
	g++ ../src/graphexception.o ../src/costsgraph.o ../src/vertexlists.o ../src/vertexlistsexception.o ../src/GraphParser/matrixparser.o -o graphvizcreator_neato -O2 -I../src/ -I../src/GraphParser/ graphvizcreator.cpp
GraphVizCreator_YED: graphvizcreator_yed.cpp
	g++ ../src/graphexception.o ../src/costsgraph.o ../src/vertexlists.o ../src/vertexlistsexception.o ../src/GraphParser/matrixparser.o -o graphvizcreator_yep -O2 -I../src/ -I../src/GraphParser/ graphvizcreator_yed.cpp

GraphProducer: ReducedGraphProducer.cpp stopwatch.h
	g++ ../src/graphexception.o ../src/costsgraph.o ../src/vertexlists.o ../src/vertexlistsexception.o ../src/graphset.o ../src/GraphParser/edgefileparser.o ../src/GraphParser/matrixparser.o ../src/CostsParser/blastparser.o ../src/probleminstanceexception.o ../src/edgereduction.o  ../src/weightedprobleminstance.o ../src/searchtreeweighted.o -o reducedgraphroducer -O2 -I. -I../src/ -I../src/GraphParser/ -I../src/CostsParser/ -I .  ReducedGraphProducer.cpp

GraphProducer_iwpec: ReducedGraphProducer_iwpec.cpp stopwatch.h
	g++ ../src/graphexception.o ../src/costsgraph.o ../src/vertexlists.o ../src/vertexlistsexception.o ../src/graphset.o ../src/GraphParser/edgefileparser.o ../src/GraphParser/matrixparser.o ../src/CostsParser/blastparser.o ../src/probleminstanceexception.o ../src/edgereduction.o  ../src/weightedprobleminstance.o ../src/searchtreeweighted.o -o reducedgraphroducer_ipec -O2 -I. -I../src/ -I../src/GraphParser/ -I../src/CostsParser/ -I .  ReducedGraphProducer_iwpec.cpp