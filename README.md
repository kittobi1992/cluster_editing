# PACE Challenge 2021 - Cluster Editing Playground

Requirements
-----------

  - A 64-bit Linux operating system.
  - A modern, ![C++14](https://img.shields.io/badge/C++-17-blue.svg?style=flat)-ready compiler such as `g++` version 7 or higher or `clang` version 11.0.3 or higher.
 - The [cmake][cmake] build system (>= 3.16).
 - The [Boost - Program Options][Boost.Program_options] library and the boost header files (>= 1.48).

Build Application
-----------

1. Clone the repository including submodules:

   ```git clone --depth=1 --recursive git@github.com:kittobi1992/cluster_editing.git```
2. Create a build directory: `mkdir build && cd build`
3. Run cmake: `cmake .. -DCMAKE_BUILD_TYPE=RELEASE` (or `cmake .. -DCMAKE_BUILD_TYPE=DEBUG` for debug build)
4. Run make: `make ClusterEditing -j8`

The binary will be located in `build/cluster_editing/application/`.

Run Application
-----------

To run the application, you can use the following command (from `build` folder):

    ./cluster_editing/application/ClusterEditing -g <path-to-graph> -p <path-to-config> --seed=<some-random-seed>

There are additional command line arguments (see `--help`) which are passed to the application with a config file (see `config/preset.ini`).

Run Tests
-----------

To build the tests, you can use the following command:

    make ClusterEditingTests -j8

To run the tests, use the following command (from `build` folder):

    ./tests/ClusterEditingTests

or

    ./tests/ClusterEditingTests --gtest_filter=<test-regex>

to run a specific test.

Assertions
-----------

Assertions are defined in the file `cluster_editing/macros.h` and enabled in debug mode per default. Usage:

```cpp
ASSERT(<cond>, "Assertion Failed!"); // or
ASSERT([&] {
  if ( <cond> ) {
    LOG << "Some detailed error message";
    return false;
  }
  ...
  return true;
}(), "Assertion failed");
```

Timings
-----------

I have added a special class to measure timings (see `cluster_editing/utils/timer.h`). Usage:

```cpp
// simple measurement
utils::Timer::instance().start_timer("<key>", "<human-readable-description>");
// DO SOME WORK
...
utils::Timer::instance().stop_timer("<key>");

// nested timings
utils::Timer::instance().start_timer("measurement_1", "Measurement 1");
utils::Timer::instance().start_timer("measurement_2", "Measurement 2");
// DO SOME WORK
...
utils::Timer::instance().stop_timer("measurement_2");

// NOTE: Timings are aggregated, if called multiple times within same scope
for ( int i = 0; i < n; ++i ) {
  utils::Timer::instance().start_timer("measurement_3", "Measurement 3");
  // DO SOME WORK
  ...
  utils::Timer::instance().stop_timer("measurement_3");
}
utils::Timer::instance().stop_timer("measurement_1");
```

At the end of each run a runtime summary is printed (currently looks like this):

```bash
Timings:
 + Import Graph                               = 0.0527741 s
 + Preprocessing                              = 8.32e-07 s
 + Multilevel Solver                          = 8.9257e-05 s
    + Coarsening                              = 1.5631e-05 s
    + Unoarsening                             = 6.07e-07 s
 + Undo Preprocessing                         = 1.049e-06 s
```


 [cmake]: http://www.cmake.org/ "CMake tool"
[Boost.Program_options]: http://www.boost.org/doc/libs/1_58_0/doc/html/program_options.html