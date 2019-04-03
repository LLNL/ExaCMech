#
# see README.md for copyright and license information
# 

cmake_minimum_required(VERSION 3.1)
if("${CMAKE_VERSION}" VERSION_GREATER 3.3)
    cmake_policy(SET CMP0057 NEW)
endif()

message(STATUS "Using cmake version ${CMAKE_VERSION}")

project(ECMECH LANGUAGES CXX)

################################
# BLT
################################
if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
     message(FATAL_ERROR "Given BLT_SOURCE_DIR does not contain SetupBLT.cmake")
endif()

set(ENABLE_GTEST OFF CACHE BOOL "")
set(ENABLE_FRUIT OFF CACHE BOOL "")

include(${BLT_SOURCE_DIR}/SetupBLT.cmake)

set(BUILD_STATIC_LIBS ON CACHE BOOL "Build static libraries")
set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries")

################################
# Include standard build system logic and options
################################
include(cmake/CMakeBasics.cmake)

################################
# Add source subdirectories
################################
add_subdirectory(src/ecmech)

if(ENABLE_TESTS)
    add_subdirectory(test)
endif(ENABLE_TESTS)
