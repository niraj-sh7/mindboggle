cmake_minimum_required(VERSION 3.5...4.4)

project(gradient)

add_executable(
  GradientMain
  GradientMain.cpp
)

install(
  TARGETS GradientMain
  RUNTIME DESTINATION bin
)

target_link_libraries(
  GradientMain
  PRIVATE
  GradientComputer
  MeshAnalyser
  PointAreaComputer
)