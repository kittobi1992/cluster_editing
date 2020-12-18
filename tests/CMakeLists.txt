add_executable(ClusterEditingTests run_tests.cpp)
target_link_libraries(ClusterEditingTests gtest gtest_main ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES})
set_property(TARGET ClusterEditingTests PROPERTY CXX_STANDARD 17)
set_property(TARGET ClusterEditingTests PROPERTY CXX_STANDARD_REQUIRED ON)

add_subdirectory(datastructures)

set(TARGETS_WANTING_ALL_SOURCES ${TARGETS_WANTING_ALL_SOURCES} ClusterEditingTests PARENT_SCOPE)
