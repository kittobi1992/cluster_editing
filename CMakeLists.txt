cmake_minimum_required(VERSION 3.16)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

project(ClusterEditing CXX C)
set(PROJECT_VENDOR "Tobias Heuer")
set(PROJECT_CONTACT "tobias.heuer@kit.edu")
set(PROJECT_URL "https://github.com/kittobi1992/cluster_editing")
set(PROJECT_DESCRIPTION "Cluster Editing")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include_directories(${PROJECT_SOURCE_DIR})

option(DISABLE_ASSERTIONS
  "Disable KaHyPar's internal assertions." OFF)

option(USE_STANDARD_ASSERTIONS
  "Use standard C++ asserts instead of custom assertions." OFF)

if(DISABLE_ASSERTIONS)
  add_compile_definitions(DISABLE_ASSERTIONS)
endif(DISABLE_ASSERTIONS)

if(USE_STANDARD_ASSERTIONS)
  add_compile_definitions(USE_STANDARD_ASSERTIONS)
endif(USE_STANDARD_ASSERTIONS)

#ensure that gmock is built before tests are executed
add_subdirectory(googletest EXCLUDE_FROM_ALL)
include_directories(SYSTEM ${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
include_directories(SYSTEM ${gtest_SOURCE_DIR}/../googlemock/include ${gtest_SOURCE_DIR}/../googlemock/)

set(BOOST_MIN_VERSION "1.48.0")
find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS program_options REQUIRED)
if(Boost_FOUND)
  include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
  message(STATUS "Boost Include: ${Boost_INCLUDE_DIR}")
  message(STATUS "Boost Library Dirs: ${Boost_LIBRARY_DIRS}")
  message(STATUS "Boost Libraries: ${Boost_LIBRARIES}")
endif()

if(NOT CMAKE_BUILD_TYPE)
  set( CMAKE_BUILD_TYPE Debug CACHE STRING
       "Choose the type of build, options are: Debug Release, RelWithDebInfo"
       FORCE )
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wuninitialized")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wfatal-errors")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wcast-qual")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Woverloaded-virtual")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wredundant-decls")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Winit-self")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPARANOID ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-function")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ltbb")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ltbbmalloc_proxy")

# lld is way faster than ld. If you have it, use it!
find_program(LLD_BIN lld)
if (LLD_BIN)
  message(STATUS "Found and will use LLVM linker " ${LLD_BIN})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fuse-ld=lld")
else()
  message(STATUS "Default linker")
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g3")
set(CMAKE_CXX_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=undefined -fsanitize=address -fno-omit-frame-pointer")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -mtune=native -march=native")
set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -mtune=native -march=native -g3 ")

if(ENABLE_PROFILE MATCHES ON)
  message(STATUS "Profiling activated")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DENABLE_PROFILE")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g3 -DENABLE_PROFILE -fno-omit-frame-pointer")
  set(PROFILE_FLAGS "-lprofiler")
endif()

message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS_RELEASE: ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "CMAKE_CXX_FLAGS_DEBUG: ${CMAKE_CXX_FLAGS_DEBUG}")

include(gmock)
enable_testing()

add_subdirectory(cluster_editing/application)
add_subdirectory(tests)

# This adds the source files. It's important that this happens after the compile targets have been added
add_subdirectory(cluster_editing)