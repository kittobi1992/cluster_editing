set(RefinementSources
      clique_remover.cpp
      clique_splitter.cpp
      node_swapper.cpp
      localized_evo.cpp
      evolutionary.cpp
      lp_refiner.cpp
      fm_refiner.cpp
      )

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
      target_sources(${modtarget} PRIVATE ${RefinementSources})
endforeach()