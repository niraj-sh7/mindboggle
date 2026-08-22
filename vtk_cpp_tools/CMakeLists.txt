cmake_minimum_required(VERSION 3.5...4.4)
project(mindboggle_surfaces)

find_package(VTK REQUIRED COMPONENTS
  CommonCore
  CommonDataModel
  CommonExecutionModel
  CommonMisc
  FiltersCore
  FiltersGeneral
  FiltersModeling
  FiltersSources
  IOLegacy
  IOGeometry
  IOMINC
  ImagingCore
  ImagingGeneral
  ImagingStencil
  ImagingMorphological
)

add_subdirectory(
travel_depth
)

add_subdirectory(
curvature
)

add_subdirectory(
surface_overlap
)

add_subdirectory(
area
)

add_subdirectory(
gradient
)

add_subdirectory(
geodesic_depth
)

#add_subdirectory(
#medial_surfaces
#)

add_library(MeshAnalyser MeshAnalyser.cpp)
target_link_libraries(MeshAnalyser PUBLIC
  VTK::CommonCore
  VTK::CommonDataModel
  VTK::CommonExecutionModel
  VTK::CommonMisc
  VTK::FiltersCore
  VTK::FiltersGeneral
  VTK::FiltersModeling
  VTK::FiltersSources
  VTK::IOLegacy
  VTK::IOGeometry
  VTK::IOMINC
  VTK::ImagingCore
  VTK::ImagingGeneral
  VTK::ImagingStencil
  VTK::ImagingMorphological
)
add_library(FsSurfaceReader FsSurfaceReader.cpp)
target_link_libraries(FsSurfaceReader PUBLIC
  VTK::CommonCore
  VTK::CommonDataModel
)

add_library(VtkFileEditor VtkFileEditor.cpp)
target_link_libraries(VtkFileEditor PUBLIC
  VTK::CommonCore
  VTK::CommonDataModel
  VTK::IOLegacy
)

add_library(TravelDepth TravelDepth.cpp)
target_link_libraries(TravelDepth PUBLIC
  VTK::CommonCore
  VTK::CommonDataModel
  MeshAnalyser
)

add_library(GradientComputer GradientComputer.cpp)
target_link_libraries(GradientComputer PUBLIC
  VTK::CommonCore
  VTK::CommonDataModel
)

add_library(PointAreaComputer PointAreaComputer.cpp)
target_link_libraries(PointAreaComputer PUBLIC
  VTK::CommonCore
  VTK::CommonDataModel
  VTK::IOLegacy
)

add_library(Overlap Overlap.cpp)
target_link_libraries(Overlap PUBLIC
  PointAreaComputer
  VTK::CommonCore
  VTK::CommonDataModel
)