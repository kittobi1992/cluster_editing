# set(ExactSources
#       ${PROJECT_SOURCE_DIR}/cluster_editing/exact/instance.h
#       ${PROJECT_SOURCE_DIR}/cluster_editing/exact/instance.cpp
#       ${PROJECT_SOURCE_DIR}/cluster_editing/exact/reductions.h
#       ${PROJECT_SOURCE_DIR}/cluster_editing/exact/reductions.cpp
#       ${PROJECT_SOURCE_DIR}/cluster_editing/exact/thomas.h
#       ${PROJECT_SOURCE_DIR}/cluster_editing/exact/thomas.cpp
#       ${PROJECT_SOURCE_DIR}/cluster_editing/exact/solver.h
#       ${PROJECT_SOURCE_DIR}/cluster_editing/exact/solver.cpp
#       ${PROJECT_SOURCE_DIR}/cluster_editing/exact/lower_bounds.h
#       ${PROJECT_SOURCE_DIR}/cluster_editing/exact/lower_bounds.cpp
#       ${PROJECT_SOURCE_DIR}/cluster_editing/exact/star_bound.h
#       ${PROJECT_SOURCE_DIR}/cluster_editing/exact/star_bound.cpp
#       ${PROJECT_SOURCE_DIR}/cluster_editing/exact/kab_bounds.h
#       ${PROJECT_SOURCE_DIR}/cluster_editing/exact/kab_bounds.cpp
# )
#
#
# foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
#     target_sources(${modtarget} PRIVATE ${ExactSources})
# endforeach()
