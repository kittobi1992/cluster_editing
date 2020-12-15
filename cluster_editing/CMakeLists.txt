add_subdirectory(context)
add_subdirectory(datastructures)
add_subdirectory(io)

set(ClusterEditingSources
      preprocessing.cpp
      multilevel.cpp
      )

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
      target_sources(${modtarget} PRIVATE ${ClusterEditingSources})
endforeach()