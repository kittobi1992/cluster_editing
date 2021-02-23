cmake_minimum_required(VERSION 3.14)
project(peace_cluster_editing)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wold-style-cast")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wcast-align -Wunused -Woverloaded-virtual -Wpedantic")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsign-conversion -Wmisleading-indentation -Wduplicated-cond")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wduplicated-branches -Wlogical-op -Wnull-dereference")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} -Wuseless-cast -Wdouble-promotion -Wformat=2")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -DNDEBUG")


message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS_RELEASE: ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "CMAKE_CXX_FLAGS_DEBUG: ${CMAKE_CXX_FLAGS_DEBUG}")


include_directories(src)
include_directories(src/CostsParser)
include_directories(src/GraphParser)
include_directories(tester)
include_directories(tools)

add_executable(weighted_cluster_editing
        src/CostsParser/blastparser.cpp
        src/CostsParser/blastparser.h
        src/CostsParser/doubleparser.cpp
        src/CostsParser/doubleparser.h
        src/GraphParser/edgefileparser.cpp
        src/GraphParser/edgefileparser.h
        src/GraphParser/matrixparser.cpp
        src/GraphParser/matrixparser.h
        src/almostclique.h
        src/array.h
        src/bnmanager.h
        src/cc-kernel.h
        src/costsgraph.cpp
        src/costsgraph.h
        src/criticalclique.h
        src/edgereduction.cpp
        src/edgereduction.h
        src/Exceptions/graphexception.cpp
        src/Exceptions/graphexception.h
        src/graphset.cpp
        src/graphset.h
        src/lowerbound.h
        src/mergereduction.h
        src/mincut.h
        src/napsack.h
        src/Exceptions/probleminstanceexception.cpp
        src/Exceptions/probleminstanceexception.h
        src/searchtreeweighted.cpp
        src/searchtreeweighted.h
        src/trianglematrix.h
        src/vertexlists.cpp
        src/vertexlists.h
        src/Exceptions/vertexlistsexception.cpp
        src/Exceptions/vertexlistsexception.h
        src/weightedprobleminstance.cpp
        src/weightedprobleminstance.h
        tools/weightedclusterediting.cpp
        )