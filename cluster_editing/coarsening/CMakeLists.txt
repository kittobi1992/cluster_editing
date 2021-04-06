set(CoarseningSources
      ${PROJECT_SOURCE_DIR}/cluster_editing/coarsening/lp_coarsener.cpp
      )

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
      target_sources(${modtarget} PRIVATE ${CoarseningSources})
endforeach()