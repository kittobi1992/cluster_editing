set(ContextSources
      context.cpp
      context_enum_classes.cpp
      )

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
      target_sources(${modtarget} PRIVATE ${ContextSources})
endforeach()