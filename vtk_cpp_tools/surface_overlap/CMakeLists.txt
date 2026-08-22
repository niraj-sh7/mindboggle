cmake_minimum_required(VERSION 3.5...4.4)

project(surface_overlap)

add_executable(
  SurfaceOverlapMain
  SurfaceOverlapMain.cpp
)

install(
  TARGETS SurfaceOverlapMain
  RUNTIME DESTINATION bin
)

target_link_libraries(
  SurfaceOverlapMain
  PRIVATE
  Overlap
  VTK::IOLegacy
  VTK::IOGeometry
  VTK::IOMINC
)