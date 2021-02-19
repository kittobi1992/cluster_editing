set(RefinementSources
      lp_refiner.cpp
        fm_refiner.cpp
      )

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
      target_sources(${modtarget} PRIVATE ${RefinementSources})
endforeach()