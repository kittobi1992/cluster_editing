set(RefinementSources
      evolutionary.cpp
      lp_refiner.cpp
      fm_refiner.cpp
      exact_refiner.cpp
      )

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
      target_sources(${modtarget} PRIVATE ${RefinementSources})
endforeach()