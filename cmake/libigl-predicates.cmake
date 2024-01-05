if (TARGET igl::predicates)
    return()
endif()

message(STATUS "Third-party (external): creating target 'igl::predicates'")

include(CPM)
CPMAddPackage(
    NAME predicates
    GITHUB_REPOSITORY libigl/libigl-predicates
    GIT_TAG 50c2149e7a520d13cd10e9aeff698bd68edd5a4f
)

set_target_properties(predicates PROPERTIES FOLDER third_party)
target_compile_options(predicates PRIVATE "-Wno-deprecated-non-prototype")

add_library(igl::predicates ALIAS predicates)
