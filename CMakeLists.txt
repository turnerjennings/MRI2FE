cmake_minimum_required(VERSION 3.18)

#project name
project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

#c++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add C++ standard library include paths for macOS
if(APPLE)
    set(CMAKE_OSX_SYSROOT "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem /usr/local/include")
    include_directories(SYSTEM "/usr/local/include")
    link_directories("/usr/local/lib")
endif()

INCLUDE_DIRECTORIES(include)


#find required extensions
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# Set up vcpkg paths if available
if(DEFINED ENV{VCPKG_TOOLCHAIN_FILE})
    message(STATUS "Using vcpkg toolchain: $ENV{VCPKG_TOOLCHAIN_FILE}")
    
    # Check for different possible vcpkg directories (both installed and packages)
    set(VCPKG_INSTALLED_DIRS
        "${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/installed/x64-windows"
        "${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/installed/x64-windows-static"
        "${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/installed/x86-windows"
        "${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/installed/x86-windows-static"
        "${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/packages/cgal_x64-windows"
        "${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/packages/eigen3_x64-windows"
        "${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/packages/gmp_x64-windows"
        "${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/packages/mpfr_x64-windows"
    )
    
    message(STATUS "Checking vcpkg directories:")
    foreach(VCPKG_DIR ${VCPKG_INSTALLED_DIRS})
        if(EXISTS "${VCPKG_DIR}")
            list(APPEND CMAKE_PREFIX_PATH "${VCPKG_DIR}")
            message(STATUS "  Found and added: ${VCPKG_DIR}")
            
            # Set package directories if they exist
            if(EXISTS "${VCPKG_DIR}/share/cgal")
                set(CGAL_DIR "${VCPKG_DIR}/share/cgal")
                message(STATUS "  Set CGAL_DIR to: ${CGAL_DIR}")
            endif()
            if(EXISTS "${VCPKG_DIR}/share/eigen3")
                set(Eigen3_DIR "${VCPKG_DIR}/share/eigen3")
                message(STATUS "  Set Eigen3_DIR to: ${Eigen3_DIR}")
            endif()
            
            # Also check for packages directory structure
            if(EXISTS "${VCPKG_DIR}/share/cgal/usage")
                set(CGAL_DIR "${VCPKG_DIR}/share/cgal")
                message(STATUS "  Set CGAL_DIR to: ${CGAL_DIR} (from packages)")
            endif()
        else()
            message(STATUS "  Not found: ${VCPKG_DIR}")
        endif()
    endforeach()
else()
    message(STATUS "VCPKG_TOOLCHAIN_FILE not defined")
endif()

# Try to find packages - let vcpkg handle the finding
find_package(CGAL CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)

# Print debug information
message(STATUS "CGAL_FOUND: ${CGAL_FOUND}")
message(STATUS "Eigen3_FOUND: ${Eigen3_FOUND}")
message(STATUS "GMP_FOUND: ${GMP_FOUND}")
message(STATUS "MPFR_FOUND: ${MPFR_FOUND}")
message(STATUS "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")

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

# Add additional include paths for macOS
if(APPLE)
    target_include_directories(_MESHUTILS PRIVATE SYSTEM
        "/usr/local/include"
    )
endif()

#compile and install
target_compile_definitions(_MESHUTILS PRIVATE VERSION_INFO=${PROJECT_VERSION})

install(TARGETS _MESHUTILS DESTINATION ${SKBUILD_PROJECT_NAME})