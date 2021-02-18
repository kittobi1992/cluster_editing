add_executable(ClusterEditing cluster_editing.cc)
target_link_libraries(ClusterEditing ${Boost_LIBRARIES})

set_property(TARGET ClusterEditing PROPERTY CXX_STANDARD 17)
set_property(TARGET ClusterEditing PROPERTY CXX_STANDARD_REQUIRED ON)

if(ENABLE_PROFILE MATCHES ON)
  target_link_libraries(ClusterEditing ${PROFILE_FLAGS})
endif()

add_executable(DrawGraph draw_graph.cc)
target_link_libraries(DrawGraph ${Boost_LIBRARIES})

set_property(TARGET DrawGraph PROPERTY CXX_STANDARD 17)
set_property(TARGET DrawGraph PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(exact exact.cc)
target_link_libraries(exact ${Boost_LIBRARIES})

set_property(TARGET exact PROPERTY CXX_STANDARD 17)
set_property(TARGET exact PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(calc_cost calc_cost.cc)
target_link_libraries(calc_cost ${Boost_LIBRARIES})
set_property(TARGET calc_cost PROPERTY CXX_STANDARD 17)
set_property(TARGET calc_cost PROPERTY CXX_STANDARD_REQUIRED ON)

set(TARGETS_WANTING_ALL_SOURCES ${TARGETS_WANTING_ALL_SOURCES} ClusterEditing DrawGraph exact calc_cost PARENT_SCOPE)
