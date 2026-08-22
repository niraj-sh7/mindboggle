cmake_minimum_required(VERSION 3.5...4.4)

project(travel_depth)

add_executable(
  TravelDepthMain
  TravelDepthMain.cpp
)

install(
  TARGETS TravelDepthMain
  RUNTIME DESTINATION bin
)

target_link_libraries(
  TravelDepthMain
  PRIVATE
  TravelDepth
  MeshAnalyser
)