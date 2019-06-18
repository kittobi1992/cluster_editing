GraphParser: matrixparser.o edgefileparser.o

matrixparser.o: matrixparser.h matrixparser.cpp ../costsgraph.o
	g++ -c -O2 -I. -I../ matrixparser.cpp

edgefileparser.o: edgefileparser.h edgefileparser.cpp ../costsgraph.o ../graphset.o
	g++ -c -O2 -I. -I../ edgefileparser.cpp
