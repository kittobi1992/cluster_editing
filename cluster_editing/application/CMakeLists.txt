
# heuristic solver
add_executable(heuristic heuristic.cc)
target_link_libraries(heuristic ${Boost_LIBRARIES})
set_property(TARGET heuristic PROPERTY CXX_STANDARD 17)
set_property(TARGET heuristic PROPERTY CXX_STANDARD_REQUIRED ON)

if(ENABLE_PROFILE MATCHES ON)
  target_link_libraries(ClusterEditing ${PROFILE_FLAGS})
endif()


add_executable(exact exact.cc)
target_link_libraries(exact ${Boost_LIBRARIES})
set_property(TARGET exact PROPERTY CXX_STANDARD 17)
set_property(TARGET exact PROPERTY CXX_STANDARD_REQUIRED ON)

set(TARGETS_WANTING_ALL_SOURCES ${TARGETS_WANTING_ALL_SOURCES} heuristic exact PARENT_SCOPE)

