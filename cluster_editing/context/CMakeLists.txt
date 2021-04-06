set(ContextSources
      ${PROJECT_SOURCE_DIR}/cluster_editing/context/context.cpp
      ${PROJECT_SOURCE_DIR}/cluster_editing/context/context_enum_classes.cpp
      )

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
      target_sources(${modtarget} PRIVATE ${ContextSources})
endforeach()