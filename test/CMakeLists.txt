cmake_minimum_required(VERSION 3.11)
project(MyCatch2Tests)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add Catch2 (assuming it's cloned into extern/Catch2)
add_subdirectory(../Catch2 ${CMAKE_BINARY_DIR}/Catch2-build)

INCLUDE_DIRECTORIES(include)
INCLUDE_DIRECTORIES(eigen)

#find the python extensions package
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 REQUIRED)

#include both the test files and the .cpp files that the tests are pulling from
set(SRC_FILES
../src/cpp/meshing.cpp
Test_meshing.cpp
)

#create executable from source files defined above
add_executable(test_util ${SRC_FILES})

#define include directory with headers in the root file
target_include_directories(test_util PRIVATE 
    "${CMAKE_CURRENT_LIST_DIR}/../include"
    "${CMAKE_CURRENT_LIST_DIR}/../eigen")
    
#link catch2
target_link_libraries(test_util PRIVATE 
    Catch2::Catch2WithMain
    pybind11::module
)



#compile and create tests
include(CTest)

add_test(NAME "test_debug_add" COMMAND test_util)