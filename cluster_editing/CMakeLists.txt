add_subdirectory(context)
add_subdirectory(datastructures)
add_subdirectory(io)
add_subdirectory(refinement)
add_subdirectory(utils)

add_subdirectory(exact)

configure_file(data_path.h.in data_path.h)

set(ClusterEditingSources
      flat.cpp
      data_path.h
      )

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
      target_include_directories(${modtarget} PRIVATE ${CMAKE_BINARY_DIR}) # to get the configured file 'data_path.h'
      target_sources(${modtarget} PRIVATE ${ClusterEditingSources})
endforeach()