set(IOSources
      ${PROJECT_SOURCE_DIR}/cluster_editing/io/command_line_options.cpp
      ${PROJECT_SOURCE_DIR}/cluster_editing/io/output.cpp
      ${PROJECT_SOURCE_DIR}/cluster_editing/io/graph_io.cpp
      ${PROJECT_SOURCE_DIR}/cluster_editing/io/ipe.cpp
      ${PROJECT_SOURCE_DIR}/cluster_editing/io/colors.cpp
      # ${PROJECT_SOURCE_DIR}/cluster_editing/io/draw_graph.cpp
      )

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
    target_sources(${modtarget} PRIVATE ${IOSources})
endforeach()
