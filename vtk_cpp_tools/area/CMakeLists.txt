cmake_minimum_required(VERSION 3.5...4.4)

project(point_area)

add_executable(
  PointAreaMain
  PointAreaMain.cpp
)

install(
  TARGETS PointAreaMain
  RUNTIME DESTINATION bin
)

target_link_libraries(
  PointAreaMain
  PRIVATE
  MeshAnalyser
  PointAreaComputer
)