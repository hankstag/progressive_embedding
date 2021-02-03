cmake_minimum_required(VERSION 3.1)
project(matchmaker)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(THIRD_PARTY_DIR ${CMAKE_SOURCE_DIR}/ext)
include(PEDependencies)

add_library(progressive_embedding 
  local_operation.cpp
  validity_check.cpp
  progressive_embedding.cpp
  local_smooth/auto_grad.cpp
  local_smooth/coloring_mesh.cpp
  local_smooth/local_smooth.cpp
)

add_library(shor_van_wyck 
  shor.cpp
  is_simple_polygon.cpp
  plot.cpp
)

add_library(matchmaker
  edge_split.cpp path_tracing.cpp target_polygon.cpp matchmaker.cpp 
  mst.cpp decompose_polygon.cpp embed_points.cpp plot.cpp
)

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

add_executable(matchmaker_bin ${MATCHMAKER_FILES})
add_executable(genus_zero_tutte_bin ${GENUSZERO_FILES})
add_executable(untangle_bin ${UNTANGLE_FILES})
add_executable(random_init_bin ${RANDOMINIT_FILES})

target_link_libraries(progressive_embedding igl::core igl::triangle igl::opengl_glfw igl::cgal igl::predicates igl::opengl_glfw_imgui tbb_static)
target_link_libraries(shor_van_wyck igl::core igl::triangle igl::cgal igl::opengl_glfw igl::opengl_glfw_imgui igl::predicates)
target_link_libraries(matchmaker shor_van_wyck)
target_include_directories(matchmaker PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(progressive_embedding SYSTEM PUBLIC ${THIRD_PARTY_DIR}/tbb/include)

target_link_libraries(matchmaker_bin progressive_embedding shor_van_wyck)
target_link_libraries(genus_zero_tutte_bin progressive_embedding shor_van_wyck)
target_link_libraries(untangle_bin progressive_embedding shor_van_wyck)
target_link_libraries(random_init_bin progressive_embedding shor_van_wyck)

