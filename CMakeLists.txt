cmake_minimum_required(VERSION 3.1)
project(matchmaker)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_CGAL              "Use CGAL"           ON)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"       ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          ON)

find_package(LIBIGL REQUIRED QUIET)

set(THIRD_PARTY_DIR ${CMAKE_SOURCE_DIR}/ext)
include(DownloadProject)
download_project(
    PROJ         tbb
    SOURCE_DIR   ${THIRD_PARTY_DIR}/tbb
    DOWNLOAD_DIR ${THIRD_PARTY_DIR}/.cache/tbb
    QUIET
    GIT_REPOSITORY https://github.com/wjakob/tbb.git
    GIT_TAG        08b4341a1893a72656467e96137f1f99d0112547
)

set(TBB_BUILD_STATIC ON CACHE BOOL " " FORCE)
set(TBB_BUILD_SHARED OFF CACHE BOOL " " FORCE)
set(TBB_BUILD_TBBMALLOC OFF CACHE BOOL " " FORCE)
set(TBB_BUILD_TBBMALLOC_PROXY OFF CACHE BOOL " " FORCE)
set(TBB_BUILD_TESTS OFF CACHE BOOL " " FORCE)

add_subdirectory(${THIRD_PARTY_DIR}/tbb tbb)
target_compile_definitions(tbb_static PUBLIC -DUSE_TBB)

add_library(progressive_embedding 
local_operation.cpp
validity_check.cpp
progressive_embedding.cpp
local_smooth/auto_grad.cpp
local_smooth/coloring_mesh.cpp
local_smooth/local_smooth.cpp)

add_library(shor_van_wyck 
shor.cpp
is_simple_polygon.cpp
plot.cpp)

set(MATCHMAKER_FILES
    main.cpp
    edge_split.cpp
    path_tracing.cpp
    target_polygon.cpp
    matchmaker.cpp
    mst.cpp
    decompose_polygon.cpp
    embed_points.cpp
    plot.cpp
    loader.cpp
)

set(GENUSZERO_FILES
  genus_zero_tutte.cpp
  plot.cpp
  loader.cpp
)

set(UNTANGLE_FILES
  untangling.cpp
  plot.cpp
  loader.cpp
)

set(RANDOMINIT_FILES
  random_init.cpp
  plot.cpp
  loader.cpp
)

set(TEST_LOCAL_SMOOTHING

  test_local_smoothing.cpp
  plot.cpp
  loader.cpp
  
)

add_executable(matchmaker_bin ${MATCHMAKER_FILES})
add_executable(genus_zero_tutte_bin ${GENUSZERO_FILES})
add_executable(untangle_bin ${UNTANGLE_FILES})
add_executable(random_init_bin ${RANDOMINIT_FILES})
add_executable(test_local_smoothing_bin ${TEST_LOCAL_SMOOTHING})

target_link_libraries(progressive_embedding igl::core igl::triangle igl::opengl_glfw igl::cgal igl::opengl_glfw_imgui tbb_static)
target_link_libraries(shor_van_wyck igl::core igl::triangle igl::cgal igl::opengl_glfw igl::opengl_glfw_imgui)
target_include_directories(progressive_embedding SYSTEM PUBLIC ${THIRD_PARTY_DIR}/tbb/include)

target_link_libraries(matchmaker_bin progressive_embedding shor_van_wyck)
target_link_libraries(genus_zero_tutte_bin progressive_embedding shor_van_wyck)
target_link_libraries(untangle_bin progressive_embedding shor_van_wyck)
target_link_libraries(random_init_bin progressive_embedding shor_van_wyck)
target_link_libraries(test_local_smoothing_bin progressive_embedding shor_van_wyck)

