set(ExactSources
      instance.h
      instance.cpp
      reductions.h
      reductions.cpp
      thomas.h
      thomas.cpp
      solver.h
      solver.cpp
    lower_bounds.h
    lower_bounds.cpp
    star_bound.h
    star_bound.cpp
)


foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
    target_sources(${modtarget} PRIVATE ${ExactSources})
endforeach()
