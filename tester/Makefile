all: ClusterEditingTester CCTester ClusterEditingTesterIWPEC ReductionTester

CCTester: cckerneltester.cpp ../tools/stopwatch.h
	g++ ../src/graphexception.o ../src/costsgraph.o ../src/vertexlists.o ../src/vertexlistsexception.o ../src/graphset.o ../src/GraphParser/edgefileparser.o ../src/GraphParser/matrixparser.o ../src/CostsParser/blastparser.o ../src/probleminstanceexception.o ../src/edgereduction.o  ../src/weightedprobleminstance.o ../src/searchtreeweighted.o -o cckerneltester -O2 -I. -I../src/ -I../src/GraphParser/ -I../src/CostsParser/ -I . cckerneltester.cpp

ClusterEditingTester: weightedclustereditingtester.cpp ../tools/stopwatch.h
	g++ ../src/graphexception.o ../src/costsgraph.o ../src/vertexlists.o ../src/vertexlistsexception.o ../src/graphset.o ../src/GraphParser/edgefileparser.o ../src/GraphParser/matrixparser.o ../src/CostsParser/blastparser.o ../src/probleminstanceexception.o ../src/edgereduction.o  ../src/weightedprobleminstance.o ../src/searchtreeweighted.o -o weightedclustereditingtester -O2 -I. -I../src/ -I../src/GraphParser/ -I../src/CostsParser/ -I . weightedclustereditingtester.cpp

ClusterEditingTesterIWPEC: weightedclustereditingtesteriwpec.cpp ../tools/stopwatch.h
	g++ ../src/graphexception.o ../src/costsgraph.o ../src/vertexlists.o ../src/vertexlistsexception.o ../src/graphset.o ../src/GraphParser/edgefileparser.o ../src/GraphParser/matrixparser.o ../src/CostsParser/blastparser.o ../src/probleminstanceexception.o ../src/edgereduction.o  ../src/weightedprobleminstance.o ../src/searchtreeweighted.o -o weightedclustereditingtesteriwpec -O2 -I. -I../src/ -I../src/GraphParser/ -I../src/CostsParser/ -I . weightedclustereditingtesteriwpec.cpp

#ReductionTester: weightedclusterreductiontester.cpp ../tools/stopwatch.h
#	g++ ../src/graphexception.o ../src/costsgraph.o ../src/vertexlists.o ../src/vertexlistsexception.o ../src/graphset.o ../src/GraphParser/edgefileparser.o ../src/GraphParser/matrixparser.o ../src/CostsParser/blastparser.o ../src/probleminstanceexception.o ../src/edgereduction.o  ../src/weightedprobleminstance.o ../src/searchtreeweighted.o -o weightedclusterreductiontester -O2 -I. -I../src/ -I../src/GraphParser/ -I../src/CostsParser/ -I . weightedclusterreductiontester.cpp
