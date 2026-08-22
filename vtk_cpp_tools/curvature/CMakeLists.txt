cmake_minimum_required(VERSION 3.5...4.4)

project(curvature)

add_executable(
  CurvatureMain
  CurvatureMain.cpp
)

install(
  TARGETS CurvatureMain
  RUNTIME DESTINATION bin
)

target_link_libraries(
  CurvatureMain
  PRIVATE
  MeshAnalyser
)