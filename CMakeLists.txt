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

include(T:/codes/vcpkg/scripts/buildsystems/vcpkg.cmake)


#find required extensions
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 REQUIRED)
find_package(CGAL REQUIRED)

#define source c++ files
set(SRC_FILES
src/cpp/meshing.cpp
)

#link libraries
pybind11_add_module(_MESHUTILS MODULE ${SRC_FILES})

target_link_libraries(_MESHUTILS PRIVATE CGAL::CGAL pybind11::module)

#compile and install
target_compile_definitions(_MESHUTILS PRIVATE VERSION_INFO=${PROJECT_VERSION})

install(TARGETS _MESHUTILS DESTINATION ${SKBUILD_PROJECT_NAME})