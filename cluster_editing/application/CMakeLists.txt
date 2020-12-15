add_executable(ClusterEditing cluster_editing.cc)
target_link_libraries(ClusterEditing ${Boost_LIBRARIES})

set_property(TARGET ClusterEditing PROPERTY CXX_STANDARD 17)
set_property(TARGET ClusterEditing PROPERTY CXX_STANDARD_REQUIRED ON)

if(ENABLE_PROFILE MATCHES ON)
  target_link_libraries(ClusterEditing ${PROFILE_FLAGS})
endif()

set(TARGETS_WANTING_ALL_SOURCES ${TARGETS_WANTING_ALL_SOURCES} ClusterEditing PARENT_SCOPE)