set(RefinementSources
      ${PROJECT_SOURCE_DIR}/cluster_editing/refinement/clique_remover.cpp
      ${PROJECT_SOURCE_DIR}/cluster_editing/refinement/clique_splitter.cpp
      ${PROJECT_SOURCE_DIR}/cluster_editing/refinement/node_swapper.cpp
      ${PROJECT_SOURCE_DIR}/cluster_editing/refinement/pseudo_fm.cpp
      ${PROJECT_SOURCE_DIR}/cluster_editing/refinement/evolutionary.cpp
      ${PROJECT_SOURCE_DIR}/cluster_editing/refinement/localized_evo.cpp
      ${PROJECT_SOURCE_DIR}/cluster_editing/refinement/lp_refiner.cpp
      ${PROJECT_SOURCE_DIR}/cluster_editing/refinement/fm_refiner.cpp
      )

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
      target_sources(${modtarget} PRIVATE ${RefinementSources})
endforeach()