cmake_minimum_required(VERSION 3.16)
project(MyCatch2Tests)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add Catch2 (assuming it's cloned into extern/Catch2)

INCLUDE_DIRECTORIES(include)


#find the python extensions package
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 REQUIRED)
find_package(Catch2 REQUIRED)
find_package(CGAL CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)

#check for Eigen3 support
include(CGAL_Eigen3_support)
if(NOT TARGET CGAL::Eigen3_support)
  message("NOTICE: Requires the Eigen3 Library")
  return()
endif()

#include both the test files and the .cpp files that the tests are pulling from
set(SRC_FILES
../src/cpp/meshing.cpp
Test_meshing.cpp
)

#create executable from source files defined above
add_executable(test_util ${SRC_FILES})

# Use precompiled headers for Catch2
target_precompile_headers(test_util PRIVATE 
    <catch2/catch_all.hpp>
)

target_compile_definitions(test_util PUBLIC TEST_RESOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/test_data/")

#define include directory with headers in the root file
target_include_directories(test_util PRIVATE 
    "${CMAKE_CURRENT_LIST_DIR}/../include"
    "${CMAKE_CURRENT_LIST_DIR}/../eigen")

#link catch2
target_link_libraries(test_util PRIVATE 
    Catch2::Catch2WithMain
    pybind11::module
    Eigen3::Eigen
    CGAL::CGAL
)

#add catch2extras
list(APPEND CMAKE_MODULE_PATH "../Catch2/extras")

#compile and create tests
include(CTest)
include(Catch)
catch_discover_tests(test_util)