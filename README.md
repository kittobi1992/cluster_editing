# KaPoCE - An Exact and Heuristic Solver for the Cluster Editing Problem

The **Ka**rlsuhe and **Po**tsdam **C**luster **E**diting framework is a software for
solving the cluster editing problem. The cluster editing problem is to transform an input graph into a cluster graph (a disjoint union of complete graphs) by performing a minimum number of edge editing operations.  An edit operation can be either adding a new edge or removing an existing edge.
This repository provides an exact and heuristic solver for the cluster editing problem.
For a detailed overview of techniques used in our solvers, we refer to the following resources:

- Exact Solver ([PDF][exact_description])
- Heuristic Solver ([PDF][heuristic_description])

![cluster-editing](https://user-images.githubusercontent.com/9654047/119774492-88069e00-bec2-11eb-8800-c4abfcacb82f.png)

Requirements
-----------

  - A 64-bit Linux operating system.
  - A modern, ![C++17](https://img.shields.io/badge/C++-17-blue.svg?style=flat)-ready compiler
 - The [cmake][cmake] build system (>= 3.16).
 - The [Boost - Program Options][Boost.Program_options] library and the boost header files (>= 1.48).

Build Application
-----------

1. Clone the repository including submodules:

   ```git clone --recursive git@github.com:kittobi1992/cluster_editing.git```
2. Create a build directory: `mkdir build && cd build`
3. Run cmake: `cmake .. -DCMAKE_BUILD_TYPE=RELEASE` (or `cmake .. -DCMAKE_BUILD_TYPE=DEBUG` for debug build)
4. Run make: `make -j 8`

Run Application
-----------

There are two applications: `heuristic` and `exact`.
Both read a cluster editing instance from stdin and print the edits to stdout.
For the input and output format, please refer to the [PACE challenge web page](https://pacechallenge.org/2021/).

A usage example from the `build` folder of the exact solver would be:

    cat ../instances/exact/exact015.gr | ./exact > edits.txt

and of the heuristic solver (also from the `build` folder):

    cat ../instances/heur/heur011.gr | ./heuristic > edits.txt

Per default, the heuristic solver runs for 10 minutes and the exact solver until an optimal solution is found.
The time limit of the solvers can be adjusted via the `--time-limit` flag (in seconds):

    cat ../instances/heur/heur011.gr | ./heuristic --time-limit=100 > edits.txt

Detailed output can be enabled via the `--enable-logging` flag:

    cat ../instances/heur/heur011.gr | ./heuristic --enable-logging=true > edits.txt


Run Tests
-----------

To run the tests, use the following command (from `build` folder):

    ./ClusterEditingTests

or

    ./ClusterEditingTests --gtest_filter=<test-regex>

to run a specific test.

Submit
-----------
Prepare CMake submission for optil.io with

    tar -czvf exact.tgz cluster_editing/ cmake/ CMakeLists.txt config/ googletest/ tests/
    tar -czvf heuristic.tgz cluster_editing/ cmake/ CMakeLists.txt config/ googletest/ tests/


[cmake]: http://www.cmake.org/ "CMake tool"
[Boost.Program_options]: http://www.boost.org/doc/libs/1_58_0/doc/html/program_options.html
[heuristic_description]: http://algo2.iti.kit.edu/heuer/kapoce/kapoce_heuristic.pdf "Heuristic Solver Description"
[exact_description]: http://algo2.iti.kit.edu/heuer/kapoce/kapoce_heuristic.pdf "Exact Solver Description"