cmake_minimum_required(VERSION 3.5...4.4)

project(geodesic_depth)

add_executable(
  GeodesicDepthMain
  GeodesicDepthMain.cpp
)

install(
  TARGETS GeodesicDepthMain
  RUNTIME DESTINATION bin
)

target_link_libraries(
  GeodesicDepthMain
  PRIVATE
  MeshAnalyser
)