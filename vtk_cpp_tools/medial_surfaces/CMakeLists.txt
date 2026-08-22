cmake_minimum_required(VERSION 3.5...4.4)

project(medial_surface)

add_library(
  MedialSurfaceComputer
  MedialSurfaceComputer.cpp
)

add_executable(
  MedialSurfaceMain
  MedialSurfaceMain.cpp
)

install(
  TARGETS MedialSurfaceMain
  RUNTIME DESTINATION bin
)

target_link_libraries(
  MedialSurfaceMain
  PRIVATE
  MedialSurfaceComputer
)