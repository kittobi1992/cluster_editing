set(DatastructureSources
      ${PROJECT_SOURCE_DIR}/cluster_editing/datastructures/graph_factory.cpp
      ${PROJECT_SOURCE_DIR}/cluster_editing/datastructures/graph.cpp
      # ${PROJECT_SOURCE_DIR}/cluster_editing/datastructures/adjacency_row.cpp
      )

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
      target_sources(${modtarget} PRIVATE ${DatastructureSources})
endforeach()