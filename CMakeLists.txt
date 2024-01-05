cmake_minimum_required(VERSION 3.14)

project(convex_hull_membership)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(libigl-predicates)

add_library(convex_hull_membership INTERFACE)
target_include_directories(convex_hull_membership INTERFACE include)
target_link_libraries(convex_hull_membership INTERFACE igl::predicates)
target_compile_features(convex_hull_membership INTERFACE cxx_std_20)

option(CONVEX_HULL_MEMBERSHIP_BUILD_TESTS "Build tests" OFF)
if (CONVEX_HULL_MEMBERSHIP_BUILD_TESTS)
    include(CTest)
    enable_testing()
    include(catch2)

    file(GLOB TEST_FILES "${CMAKE_CURRENT_LIST_DIR}/tests/*.cpp")
    add_executable(convex_hull_membership_tests ${TEST_FILES})
    target_link_libraries(convex_hull_membership_tests PRIVATE
        convex_hull_membership
        Catch2::Catch2WithMain)
    target_compile_definitions(convex_hull_membership_tests PRIVATE CATCH_CONFIG_ENABLE_BENCHMARKING)
    target_compile_features(convex_hull_membership_tests PRIVATE cxx_std_20)
endif()
