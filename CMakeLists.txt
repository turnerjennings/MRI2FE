cmake_minimum_required(VERSION 3.18)

#project name
project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

#c++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

INCLUDE_DIRECTORIES(
    include
    eigen
)

#find required extensions
find_package(python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 REQUIRED)
find_package(CGAL REQUIRED COMPONENTS Core)

#define source c++ files
set(SRC_FILES
src/cpp/meshing.cpp
)

#link libraries
Python_add_library(_MESHUTILS MODULE ${SRC_FILES} WITH_SAOBI)
target_link_libraries(_MESHUTILS CGAL::CGAL)
target_link_libraries(_MESHUTILS PRIVATE pybind11::module)

#compile and install
target_compile_definitions(_MESHUTILS PRIVATE VERSION_INFO=${PROJECT_VERSION})

install(TARGETS _MESHUTILS DESTINATION ${SKBUILD_PROJECT_NAME})