set(DatastructureSources
      graph_factory.cpp
      )

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
      target_sources(${modtarget} PRIVATE ${DatastructureSources})
endforeach()