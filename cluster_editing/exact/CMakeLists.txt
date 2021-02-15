set(ExactSources
      instance.h
      reductions.h
      reductions.cpp
      thomas.h
      thomas.cpp
      solver.h
      solver.cpp
)


foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
    target_sources(${modtarget} PRIVATE ${ExactSources})
endforeach()
