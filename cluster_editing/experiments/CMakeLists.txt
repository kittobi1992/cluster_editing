
find_package(girgs)
if(NOT girgs_FOUND)
  return()
endif()

find_package(OpenMP REQUIRED) # needed for girgs

add_executable(experiments main.cc)
target_link_libraries(experiments girgs::girgs ${Boost_LIBRARIES})
set_property(TARGET experiments PROPERTY CXX_STANDARD 17)
set_property(TARGET experiments PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(gexp gexp.cc)
target_link_libraries(gexp girgs::girgs ${Boost_LIBRARIES})
set_property(TARGET gexp PROPERTY CXX_STANDARD 17)
set_property(TARGET gexp PROPERTY CXX_STANDARD_REQUIRED ON)

set(TARGETS_WANTING_ALL_SOURCES ${TARGETS_WANTING_ALL_SOURCES} experiments gexp PARENT_SCOPE)
