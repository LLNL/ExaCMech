#
# see README.md for copyright and license information
# 

cmake_minimum_required(VERSION 3.1)
if("${CMAKE_VERSION}" VERSION_GREATER 3.3)
    cmake_policy(SET CMP0057 NEW)
endif()

message(STATUS "Using cmake version ${CMAKE_VERSION}")

project(ECMECH LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

################################
# BLT
################################
if (DEFINED BLT_SOURCE_DIR)
    # Support having a shared BLT outside of the repository if given a BLT_SOURCE_DIR

    if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
        message(FATAL_ERROR "Given BLT_SOURCE_DIR does not contain SetupBLT.cmake")
    endif()
else()
    # Use internal BLT if no BLT_SOURCE_DIR is given

    set(BLT_SOURCE_DIR "${PROJECT_SOURCE_DIR}/cmake/blt" CACHE PATH "")
    if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
        message(FATAL_ERROR
            "The BLT submodule is not present. "
            "Run the following two commands in your git repository: \n"
            "    git submodule init\n"
            "    git submodule update\n"
            "Or you can point to an outside of source BLT repository by "
            "defining the variable BLT_SOURCE_DIR on your command line " )
    endif()
endif()

set(ENABLE_FRUIT OFF CACHE BOOL "")
if(ENABLE_TESTS)
    set(ENABLE_GTEST ON CACHE BOOL "" FORCE)
endif(ENABLE_TESTS)

include(${BLT_SOURCE_DIR}/SetupBLT.cmake)

set(BUILD_STATIC_LIBS ON CACHE BOOL "Build static libraries")
set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries")

################################
# Include standard build system logic and options
################################
include(cmake/CMakeBasics.cmake)

if(ENABLE_CUDA)
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -restrict -arch ${CUDA_ARCH} --expt-extended-lambda --expt-relaxed-constexpr")
endif()

################################
# SNLS
################################
set(USE_BUILT_SNLS OFF CACHE BOOL "")
if (DEFINED SNLS_DIR)
    # Support a built SNLS outside of the repository
    message(STATUS "Using pre-built SNLS: ${SNLS_DIR}")
    set(USE_BUILT_SNLS ON CACHE BOOL "" FORCE)

else()
    # Use internal SNLS submodule if no SNLS_DIR is given
    set(SNLS_SOURCE_DIR "${PROJECT_SOURCE_DIR}/snls" CACHE PATH "")
    message(STATUS "Using submodule SNLS: ${SNLS_SOURCE_DIR}")
    if (NOT EXISTS ${SNLS_SOURCE_DIR}/CMakeLists.txt)
        message(FATAL_ERROR
            "The SNLS submodule is not present. "
            "Run the following two commands in your git repository: \n"
            "    git submodule init\n"
            "    git submodule update\n"
            "Or you can point to an already built SNLS with SNLS_DIR." )
    endif()

    add_subdirectory(${SNLS_SOURCE_DIR})
endif()
#
set(SNLS_SHARED_TARGET snls CACHE STRING "")
set(SNLS_STATIC_TARGET snls CACHE STRING "")

################################
# Add source subdirectories
################################
add_subdirectory(src/ecmech)

if(ENABLE_TESTS)
    add_subdirectory(test)
endif(ENABLE_TESTS)

if(ENABLE_MINIAPPS)
    add_subdirectory(miniapp)
endif(ENABLE_MINIAPPS)
