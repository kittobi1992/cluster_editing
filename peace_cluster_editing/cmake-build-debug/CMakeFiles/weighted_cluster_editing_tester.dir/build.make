# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/clion-2019.3.2/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /opt/clion-2019.3.2/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jonas/data/projects/hiwi/peace_cluster_editing

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jonas/data/projects/hiwi/peace_cluster_editing/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/weighted_cluster_editing_tester.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/weighted_cluster_editing_tester.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/weighted_cluster_editing_tester.dir/flags.make

CMakeFiles/weighted_cluster_editing_tester.dir/src/CostsParser/blastparser.cpp.o: CMakeFiles/weighted_cluster_editing_tester.dir/flags.make
CMakeFiles/weighted_cluster_editing_tester.dir/src/CostsParser/blastparser.cpp.o: ../src/CostsParser/blastparser.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jonas/data/projects/hiwi/peace_cluster_editing/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/weighted_cluster_editing_tester.dir/src/CostsParser/blastparser.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/weighted_cluster_editing_tester.dir/src/CostsParser/blastparser.cpp.o -c /home/jonas/data/projects/hiwi/peace_cluster_editing/src/CostsParser/blastparser.cpp

CMakeFiles/weighted_cluster_editing_tester.dir/src/CostsParser/blastparser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/weighted_cluster_editing_tester.dir/src/CostsParser/blastparser.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jonas/data/projects/hiwi/peace_cluster_editing/src/CostsParser/blastparser.cpp > CMakeFiles/weighted_cluster_editing_tester.dir/src/CostsParser/blastparser.cpp.i

CMakeFiles/weighted_cluster_editing_tester.dir/src/CostsParser/blastparser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/weighted_cluster_editing_tester.dir/src/CostsParser/blastparser.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jonas/data/projects/hiwi/peace_cluster_editing/src/CostsParser/blastparser.cpp -o CMakeFiles/weighted_cluster_editing_tester.dir/src/CostsParser/blastparser.cpp.s

CMakeFiles/weighted_cluster_editing_tester.dir/src/CostsParser/doubleparser.cpp.o: CMakeFiles/weighted_cluster_editing_tester.dir/flags.make
CMakeFiles/weighted_cluster_editing_tester.dir/src/CostsParser/doubleparser.cpp.o: ../src/CostsParser/doubleparser.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jonas/data/projects/hiwi/peace_cluster_editing/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/weighted_cluster_editing_tester.dir/src/CostsParser/doubleparser.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/weighted_cluster_editing_tester.dir/src/CostsParser/doubleparser.cpp.o -c /home/jonas/data/projects/hiwi/peace_cluster_editing/src/CostsParser/doubleparser.cpp

CMakeFiles/weighted_cluster_editing_tester.dir/src/CostsParser/doubleparser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/weighted_cluster_editing_tester.dir/src/CostsParser/doubleparser.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jonas/data/projects/hiwi/peace_cluster_editing/src/CostsParser/doubleparser.cpp > CMakeFiles/weighted_cluster_editing_tester.dir/src/CostsParser/doubleparser.cpp.i

CMakeFiles/weighted_cluster_editing_tester.dir/src/CostsParser/doubleparser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/weighted_cluster_editing_tester.dir/src/CostsParser/doubleparser.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jonas/data/projects/hiwi/peace_cluster_editing/src/CostsParser/doubleparser.cpp -o CMakeFiles/weighted_cluster_editing_tester.dir/src/CostsParser/doubleparser.cpp.s

CMakeFiles/weighted_cluster_editing_tester.dir/src/GraphParser/edgefileparser.cpp.o: CMakeFiles/weighted_cluster_editing_tester.dir/flags.make
CMakeFiles/weighted_cluster_editing_tester.dir/src/GraphParser/edgefileparser.cpp.o: ../src/GraphParser/edgefileparser.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jonas/data/projects/hiwi/peace_cluster_editing/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/weighted_cluster_editing_tester.dir/src/GraphParser/edgefileparser.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/weighted_cluster_editing_tester.dir/src/GraphParser/edgefileparser.cpp.o -c /home/jonas/data/projects/hiwi/peace_cluster_editing/src/GraphParser/edgefileparser.cpp

CMakeFiles/weighted_cluster_editing_tester.dir/src/GraphParser/edgefileparser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/weighted_cluster_editing_tester.dir/src/GraphParser/edgefileparser.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jonas/data/projects/hiwi/peace_cluster_editing/src/GraphParser/edgefileparser.cpp > CMakeFiles/weighted_cluster_editing_tester.dir/src/GraphParser/edgefileparser.cpp.i

CMakeFiles/weighted_cluster_editing_tester.dir/src/GraphParser/edgefileparser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/weighted_cluster_editing_tester.dir/src/GraphParser/edgefileparser.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jonas/data/projects/hiwi/peace_cluster_editing/src/GraphParser/edgefileparser.cpp -o CMakeFiles/weighted_cluster_editing_tester.dir/src/GraphParser/edgefileparser.cpp.s

CMakeFiles/weighted_cluster_editing_tester.dir/src/GraphParser/matrixparser.cpp.o: CMakeFiles/weighted_cluster_editing_tester.dir/flags.make
CMakeFiles/weighted_cluster_editing_tester.dir/src/GraphParser/matrixparser.cpp.o: ../src/GraphParser/matrixparser.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jonas/data/projects/hiwi/peace_cluster_editing/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/weighted_cluster_editing_tester.dir/src/GraphParser/matrixparser.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/weighted_cluster_editing_tester.dir/src/GraphParser/matrixparser.cpp.o -c /home/jonas/data/projects/hiwi/peace_cluster_editing/src/GraphParser/matrixparser.cpp

CMakeFiles/weighted_cluster_editing_tester.dir/src/GraphParser/matrixparser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/weighted_cluster_editing_tester.dir/src/GraphParser/matrixparser.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jonas/data/projects/hiwi/peace_cluster_editing/src/GraphParser/matrixparser.cpp > CMakeFiles/weighted_cluster_editing_tester.dir/src/GraphParser/matrixparser.cpp.i

CMakeFiles/weighted_cluster_editing_tester.dir/src/GraphParser/matrixparser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/weighted_cluster_editing_tester.dir/src/GraphParser/matrixparser.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jonas/data/projects/hiwi/peace_cluster_editing/src/GraphParser/matrixparser.cpp -o CMakeFiles/weighted_cluster_editing_tester.dir/src/GraphParser/matrixparser.cpp.s

CMakeFiles/weighted_cluster_editing_tester.dir/src/costsgraph.cpp.o: CMakeFiles/weighted_cluster_editing_tester.dir/flags.make
CMakeFiles/weighted_cluster_editing_tester.dir/src/costsgraph.cpp.o: ../src/costsgraph.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jonas/data/projects/hiwi/peace_cluster_editing/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/weighted_cluster_editing_tester.dir/src/costsgraph.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/weighted_cluster_editing_tester.dir/src/costsgraph.cpp.o -c /home/jonas/data/projects/hiwi/peace_cluster_editing/src/costsgraph.cpp

CMakeFiles/weighted_cluster_editing_tester.dir/src/costsgraph.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/weighted_cluster_editing_tester.dir/src/costsgraph.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jonas/data/projects/hiwi/peace_cluster_editing/src/costsgraph.cpp > CMakeFiles/weighted_cluster_editing_tester.dir/src/costsgraph.cpp.i

CMakeFiles/weighted_cluster_editing_tester.dir/src/costsgraph.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/weighted_cluster_editing_tester.dir/src/costsgraph.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jonas/data/projects/hiwi/peace_cluster_editing/src/costsgraph.cpp -o CMakeFiles/weighted_cluster_editing_tester.dir/src/costsgraph.cpp.s

CMakeFiles/weighted_cluster_editing_tester.dir/src/edgereduction.cpp.o: CMakeFiles/weighted_cluster_editing_tester.dir/flags.make
CMakeFiles/weighted_cluster_editing_tester.dir/src/edgereduction.cpp.o: ../src/edgereduction.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jonas/data/projects/hiwi/peace_cluster_editing/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/weighted_cluster_editing_tester.dir/src/edgereduction.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/weighted_cluster_editing_tester.dir/src/edgereduction.cpp.o -c /home/jonas/data/projects/hiwi/peace_cluster_editing/src/edgereduction.cpp

CMakeFiles/weighted_cluster_editing_tester.dir/src/edgereduction.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/weighted_cluster_editing_tester.dir/src/edgereduction.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jonas/data/projects/hiwi/peace_cluster_editing/src/edgereduction.cpp > CMakeFiles/weighted_cluster_editing_tester.dir/src/edgereduction.cpp.i

CMakeFiles/weighted_cluster_editing_tester.dir/src/edgereduction.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/weighted_cluster_editing_tester.dir/src/edgereduction.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jonas/data/projects/hiwi/peace_cluster_editing/src/edgereduction.cpp -o CMakeFiles/weighted_cluster_editing_tester.dir/src/edgereduction.cpp.s

CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/graphexception.cpp.o: CMakeFiles/weighted_cluster_editing_tester.dir/flags.make
CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/graphexception.cpp.o: ../src/Exceptions/graphexception.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jonas/data/projects/hiwi/peace_cluster_editing/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/graphexception.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/graphexception.cpp.o -c /home/jonas/data/projects/hiwi/peace_cluster_editing/src/Exceptions/graphexception.cpp

CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/graphexception.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/graphexception.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jonas/data/projects/hiwi/peace_cluster_editing/src/Exceptions/graphexception.cpp > CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/graphexception.cpp.i

CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/graphexception.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/graphexception.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jonas/data/projects/hiwi/peace_cluster_editing/src/Exceptions/graphexception.cpp -o CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/graphexception.cpp.s

CMakeFiles/weighted_cluster_editing_tester.dir/src/graphset.cpp.o: CMakeFiles/weighted_cluster_editing_tester.dir/flags.make
CMakeFiles/weighted_cluster_editing_tester.dir/src/graphset.cpp.o: ../src/graphset.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jonas/data/projects/hiwi/peace_cluster_editing/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/weighted_cluster_editing_tester.dir/src/graphset.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/weighted_cluster_editing_tester.dir/src/graphset.cpp.o -c /home/jonas/data/projects/hiwi/peace_cluster_editing/src/graphset.cpp

CMakeFiles/weighted_cluster_editing_tester.dir/src/graphset.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/weighted_cluster_editing_tester.dir/src/graphset.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jonas/data/projects/hiwi/peace_cluster_editing/src/graphset.cpp > CMakeFiles/weighted_cluster_editing_tester.dir/src/graphset.cpp.i

CMakeFiles/weighted_cluster_editing_tester.dir/src/graphset.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/weighted_cluster_editing_tester.dir/src/graphset.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jonas/data/projects/hiwi/peace_cluster_editing/src/graphset.cpp -o CMakeFiles/weighted_cluster_editing_tester.dir/src/graphset.cpp.s

CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/probleminstanceexception.cpp.o: CMakeFiles/weighted_cluster_editing_tester.dir/flags.make
CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/probleminstanceexception.cpp.o: ../src/Exceptions/probleminstanceexception.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jonas/data/projects/hiwi/peace_cluster_editing/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/probleminstanceexception.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/probleminstanceexception.cpp.o -c /home/jonas/data/projects/hiwi/peace_cluster_editing/src/Exceptions/probleminstanceexception.cpp

CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/probleminstanceexception.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/probleminstanceexception.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jonas/data/projects/hiwi/peace_cluster_editing/src/Exceptions/probleminstanceexception.cpp > CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/probleminstanceexception.cpp.i

CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/probleminstanceexception.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/probleminstanceexception.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jonas/data/projects/hiwi/peace_cluster_editing/src/Exceptions/probleminstanceexception.cpp -o CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/probleminstanceexception.cpp.s

CMakeFiles/weighted_cluster_editing_tester.dir/src/searchtreeweighted.cpp.o: CMakeFiles/weighted_cluster_editing_tester.dir/flags.make
CMakeFiles/weighted_cluster_editing_tester.dir/src/searchtreeweighted.cpp.o: ../src/searchtreeweighted.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jonas/data/projects/hiwi/peace_cluster_editing/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/weighted_cluster_editing_tester.dir/src/searchtreeweighted.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/weighted_cluster_editing_tester.dir/src/searchtreeweighted.cpp.o -c /home/jonas/data/projects/hiwi/peace_cluster_editing/src/searchtreeweighted.cpp

CMakeFiles/weighted_cluster_editing_tester.dir/src/searchtreeweighted.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/weighted_cluster_editing_tester.dir/src/searchtreeweighted.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jonas/data/projects/hiwi/peace_cluster_editing/src/searchtreeweighted.cpp > CMakeFiles/weighted_cluster_editing_tester.dir/src/searchtreeweighted.cpp.i

CMakeFiles/weighted_cluster_editing_tester.dir/src/searchtreeweighted.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/weighted_cluster_editing_tester.dir/src/searchtreeweighted.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jonas/data/projects/hiwi/peace_cluster_editing/src/searchtreeweighted.cpp -o CMakeFiles/weighted_cluster_editing_tester.dir/src/searchtreeweighted.cpp.s

CMakeFiles/weighted_cluster_editing_tester.dir/src/vertexlists.cpp.o: CMakeFiles/weighted_cluster_editing_tester.dir/flags.make
CMakeFiles/weighted_cluster_editing_tester.dir/src/vertexlists.cpp.o: ../src/vertexlists.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jonas/data/projects/hiwi/peace_cluster_editing/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/weighted_cluster_editing_tester.dir/src/vertexlists.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/weighted_cluster_editing_tester.dir/src/vertexlists.cpp.o -c /home/jonas/data/projects/hiwi/peace_cluster_editing/src/vertexlists.cpp

CMakeFiles/weighted_cluster_editing_tester.dir/src/vertexlists.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/weighted_cluster_editing_tester.dir/src/vertexlists.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jonas/data/projects/hiwi/peace_cluster_editing/src/vertexlists.cpp > CMakeFiles/weighted_cluster_editing_tester.dir/src/vertexlists.cpp.i

CMakeFiles/weighted_cluster_editing_tester.dir/src/vertexlists.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/weighted_cluster_editing_tester.dir/src/vertexlists.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jonas/data/projects/hiwi/peace_cluster_editing/src/vertexlists.cpp -o CMakeFiles/weighted_cluster_editing_tester.dir/src/vertexlists.cpp.s

CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/vertexlistsexception.cpp.o: CMakeFiles/weighted_cluster_editing_tester.dir/flags.make
CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/vertexlistsexception.cpp.o: ../src/Exceptions/vertexlistsexception.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jonas/data/projects/hiwi/peace_cluster_editing/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/vertexlistsexception.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/vertexlistsexception.cpp.o -c /home/jonas/data/projects/hiwi/peace_cluster_editing/src/Exceptions/vertexlistsexception.cpp

CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/vertexlistsexception.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/vertexlistsexception.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jonas/data/projects/hiwi/peace_cluster_editing/src/Exceptions/vertexlistsexception.cpp > CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/vertexlistsexception.cpp.i

CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/vertexlistsexception.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/vertexlistsexception.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jonas/data/projects/hiwi/peace_cluster_editing/src/Exceptions/vertexlistsexception.cpp -o CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/vertexlistsexception.cpp.s

CMakeFiles/weighted_cluster_editing_tester.dir/src/weightedprobleminstance.cpp.o: CMakeFiles/weighted_cluster_editing_tester.dir/flags.make
CMakeFiles/weighted_cluster_editing_tester.dir/src/weightedprobleminstance.cpp.o: ../src/weightedprobleminstance.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jonas/data/projects/hiwi/peace_cluster_editing/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object CMakeFiles/weighted_cluster_editing_tester.dir/src/weightedprobleminstance.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/weighted_cluster_editing_tester.dir/src/weightedprobleminstance.cpp.o -c /home/jonas/data/projects/hiwi/peace_cluster_editing/src/weightedprobleminstance.cpp

CMakeFiles/weighted_cluster_editing_tester.dir/src/weightedprobleminstance.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/weighted_cluster_editing_tester.dir/src/weightedprobleminstance.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jonas/data/projects/hiwi/peace_cluster_editing/src/weightedprobleminstance.cpp > CMakeFiles/weighted_cluster_editing_tester.dir/src/weightedprobleminstance.cpp.i

CMakeFiles/weighted_cluster_editing_tester.dir/src/weightedprobleminstance.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/weighted_cluster_editing_tester.dir/src/weightedprobleminstance.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jonas/data/projects/hiwi/peace_cluster_editing/src/weightedprobleminstance.cpp -o CMakeFiles/weighted_cluster_editing_tester.dir/src/weightedprobleminstance.cpp.s

CMakeFiles/weighted_cluster_editing_tester.dir/tester/weightedclustereditingtester.cpp.o: CMakeFiles/weighted_cluster_editing_tester.dir/flags.make
CMakeFiles/weighted_cluster_editing_tester.dir/tester/weightedclustereditingtester.cpp.o: ../tester/weightedclustereditingtester.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jonas/data/projects/hiwi/peace_cluster_editing/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building CXX object CMakeFiles/weighted_cluster_editing_tester.dir/tester/weightedclustereditingtester.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/weighted_cluster_editing_tester.dir/tester/weightedclustereditingtester.cpp.o -c /home/jonas/data/projects/hiwi/peace_cluster_editing/tester/weightedclustereditingtester.cpp

CMakeFiles/weighted_cluster_editing_tester.dir/tester/weightedclustereditingtester.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/weighted_cluster_editing_tester.dir/tester/weightedclustereditingtester.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jonas/data/projects/hiwi/peace_cluster_editing/tester/weightedclustereditingtester.cpp > CMakeFiles/weighted_cluster_editing_tester.dir/tester/weightedclustereditingtester.cpp.i

CMakeFiles/weighted_cluster_editing_tester.dir/tester/weightedclustereditingtester.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/weighted_cluster_editing_tester.dir/tester/weightedclustereditingtester.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jonas/data/projects/hiwi/peace_cluster_editing/tester/weightedclustereditingtester.cpp -o CMakeFiles/weighted_cluster_editing_tester.dir/tester/weightedclustereditingtester.cpp.s

# Object files for target weighted_cluster_editing_tester
weighted_cluster_editing_tester_OBJECTS = \
"CMakeFiles/weighted_cluster_editing_tester.dir/src/CostsParser/blastparser.cpp.o" \
"CMakeFiles/weighted_cluster_editing_tester.dir/src/CostsParser/doubleparser.cpp.o" \
"CMakeFiles/weighted_cluster_editing_tester.dir/src/GraphParser/edgefileparser.cpp.o" \
"CMakeFiles/weighted_cluster_editing_tester.dir/src/GraphParser/matrixparser.cpp.o" \
"CMakeFiles/weighted_cluster_editing_tester.dir/src/costsgraph.cpp.o" \
"CMakeFiles/weighted_cluster_editing_tester.dir/src/edgereduction.cpp.o" \
"CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/graphexception.cpp.o" \
"CMakeFiles/weighted_cluster_editing_tester.dir/src/graphset.cpp.o" \
"CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/probleminstanceexception.cpp.o" \
"CMakeFiles/weighted_cluster_editing_tester.dir/src/searchtreeweighted.cpp.o" \
"CMakeFiles/weighted_cluster_editing_tester.dir/src/vertexlists.cpp.o" \
"CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/vertexlistsexception.cpp.o" \
"CMakeFiles/weighted_cluster_editing_tester.dir/src/weightedprobleminstance.cpp.o" \
"CMakeFiles/weighted_cluster_editing_tester.dir/tester/weightedclustereditingtester.cpp.o"

# External object files for target weighted_cluster_editing_tester
weighted_cluster_editing_tester_EXTERNAL_OBJECTS =

weighted_cluster_editing_tester: CMakeFiles/weighted_cluster_editing_tester.dir/src/CostsParser/blastparser.cpp.o
weighted_cluster_editing_tester: CMakeFiles/weighted_cluster_editing_tester.dir/src/CostsParser/doubleparser.cpp.o
weighted_cluster_editing_tester: CMakeFiles/weighted_cluster_editing_tester.dir/src/GraphParser/edgefileparser.cpp.o
weighted_cluster_editing_tester: CMakeFiles/weighted_cluster_editing_tester.dir/src/GraphParser/matrixparser.cpp.o
weighted_cluster_editing_tester: CMakeFiles/weighted_cluster_editing_tester.dir/src/costsgraph.cpp.o
weighted_cluster_editing_tester: CMakeFiles/weighted_cluster_editing_tester.dir/src/edgereduction.cpp.o
weighted_cluster_editing_tester: CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/graphexception.cpp.o
weighted_cluster_editing_tester: CMakeFiles/weighted_cluster_editing_tester.dir/src/graphset.cpp.o
weighted_cluster_editing_tester: CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/probleminstanceexception.cpp.o
weighted_cluster_editing_tester: CMakeFiles/weighted_cluster_editing_tester.dir/src/searchtreeweighted.cpp.o
weighted_cluster_editing_tester: CMakeFiles/weighted_cluster_editing_tester.dir/src/vertexlists.cpp.o
weighted_cluster_editing_tester: CMakeFiles/weighted_cluster_editing_tester.dir/src/Exceptions/vertexlistsexception.cpp.o
weighted_cluster_editing_tester: CMakeFiles/weighted_cluster_editing_tester.dir/src/weightedprobleminstance.cpp.o
weighted_cluster_editing_tester: CMakeFiles/weighted_cluster_editing_tester.dir/tester/weightedclustereditingtester.cpp.o
weighted_cluster_editing_tester: CMakeFiles/weighted_cluster_editing_tester.dir/build.make
weighted_cluster_editing_tester: CMakeFiles/weighted_cluster_editing_tester.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jonas/data/projects/hiwi/peace_cluster_editing/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Linking CXX executable weighted_cluster_editing_tester"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/weighted_cluster_editing_tester.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/weighted_cluster_editing_tester.dir/build: weighted_cluster_editing_tester

.PHONY : CMakeFiles/weighted_cluster_editing_tester.dir/build

CMakeFiles/weighted_cluster_editing_tester.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/weighted_cluster_editing_tester.dir/cmake_clean.cmake
.PHONY : CMakeFiles/weighted_cluster_editing_tester.dir/clean

CMakeFiles/weighted_cluster_editing_tester.dir/depend:
	cd /home/jonas/data/projects/hiwi/peace_cluster_editing/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jonas/data/projects/hiwi/peace_cluster_editing /home/jonas/data/projects/hiwi/peace_cluster_editing /home/jonas/data/projects/hiwi/peace_cluster_editing/cmake-build-debug /home/jonas/data/projects/hiwi/peace_cluster_editing/cmake-build-debug /home/jonas/data/projects/hiwi/peace_cluster_editing/cmake-build-debug/CMakeFiles/weighted_cluster_editing_tester.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/weighted_cluster_editing_tester.dir/depend

