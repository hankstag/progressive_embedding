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

add_compile_definitions(CGAL_DISABLE_ROUNDING_MATH_CHECK)

add_library(matchmaker
  edge_split.cpp path_tracing.cpp target_polygon.cpp matchmaker.cpp 
  mst.cpp decompose_polygon.cpp embed_points.cpp plot.cpp validity_check.cpp
)

add_library(property
  Property.cpp
)

target_link_libraries(progressive_embedding igl::core igl::triangle igl::opengl_glfw igl::cgal igl::predicates igl::opengl_glfw_imgui tbb_static)
target_link_libraries(shor_van_wyck igl::core igl::triangle igl::cgal igl::opengl_glfw igl::opengl_glfw_imgui igl::predicates)
target_link_libraries(property igl::core)
target_link_libraries(matchmaker shor_van_wyck progressive_embedding property)
target_include_directories(matchmaker PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(property PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
