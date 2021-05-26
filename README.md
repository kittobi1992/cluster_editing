# PACE Challenge 2021 - Cluster Editing

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
