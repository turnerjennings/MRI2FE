cmake_minimum_required(VERSION 3.18)

#project name
project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

#c++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

INCLUDE_DIRECTORIES(include)


#find required extensions
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)
find_package(CGAL CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)

include(CGAL_Eigen3_support)
if(NOT TARGET CGAL::Eigen3_support)
  message("NOTICE: All examples require the Eigen3 library, and will not be compiled.")
  return()
endif()


#define source c++ files
set(SRC_FILES
src/cpp/meshing.cpp
)

set(PYBIND11_FINDPYTHON ON)

#link libraries
pybind11_add_module(_MESHUTILS MODULE ${SRC_FILES})

target_link_libraries(_MESHUTILS PRIVATE 
    Eigen3::Eigen
    CGAL::CGAL
    ${GMP_LIBRARIES}
    ${MPFR_LIBRARIES}
)

target_include_directories(_MESHUTILS PRIVATE 
    ${GMP_INCLUDE_DIR}
    ${MPFR_INCLUDE_DIR}
)

#compile and install
target_compile_definitions(_MESHUTILS PRIVATE VERSION_INFO=${PROJECT_VERSION})

install(TARGETS _MESHUTILS DESTINATION ${SKBUILD_PROJECT_NAME})