###############################################################################
#
# Setup RAJA
# This file defines:
#  RAJA_FOUND - If RAJA was found
#  RAJA_INCLUDE_DIRS - The RAJA include directories
#  RAJA_LIBRARY - The RAJA library

# first Check for RAJA_DIR

if(NOT RAJA_DIR)
    MESSAGE(FATAL_ERROR "Could not find RAJA. RAJA support needs explicit RAJA_DIR")
endif()

if (NOT RAJA_CONFIG_CMAKE)
   set(RAJA_CONFIG_CMAKE "${RAJA_DIR}/share/raja/cmake/raja-config.cmake")
endif()
if (EXISTS "${RAJA_CONFIG_CMAKE}")
   include("${RAJA_CONFIG_CMAKE}")
endif()
if (NOT RAJA_RELEASE_CMAKE)
   set(RAJA_RELEASE_CMAKE "${RAJA_DIR}/share/raja/cmake/raja-release.cmake")
endif()
if (EXISTS "${RAJA_RELEASE_CMAKE}")
   include("${RAJA_RELEASE_CMAKE}")
endif()

find_package(RAJA REQUIRED)

# message("-- RAJA version information: ${RAJA_VERSION_MAJOR}.${RAJA_VERSION_MINOR}")

if( RAJA_VERSION_MINOR GREATER 10 OR RAJA_VERSION_MAJOR GREATER 0)
   find_package(camp REQUIRED)
   blt_register_library(NAME camp
                        INCLUDES ${camp_INSTALL_PREFIX}/include)
   set(ENABLE_CAMP ON CACHE BOOL "")
endif()

# A more elaborate version of CAMP
#
# ################################
# # CAMP (required)
# ################################
# if (NOT TARGET camp)
#    if (CAMP_DIR)
#       set(camp_DIR ${CAMP_DIR}/lib/cmake/camp)
#    endif ()
# 
#    find_package(camp QUIET)
# 
#    if (camp_FOUND)
#       message(STATUS "CARE: Using external CAMP")
#       set_target_properties(camp PROPERTIES IMPORTED_GLOBAL TRUE)
#    else ()
#       message(STATUS "CARE: Using CAMP submodule")
# 
#       if (NOT EXISTS ${PROJECT_SOURCE_DIR}/tpl/camp/CMakeLists.txt)
#          message(FATAL_ERROR "CARE: CAMP submodule not initialized. Run 'git submodule update --init' in the git repository or set camp_DIR or CAMP_DIR to use an external build of CAMP.")
#       else ()
#          add_subdirectory(${PROJECT_SOURCE_DIR}/tpl/camp)
#       endif ()
#    endif ()
# 
#    if (ENABLE_CUDA)
#       blt_add_target_definitions(TO camp
#                                  SCOPE INTERFACE
#                                  TARGET_DEFINITIONS CAMP_HAVE_CUDA)
#    endif ()
# endif ()


if(RAJA_CONFIG_LOADED)
   if(ENABLE_OPENMP)
      set(BLT_CXX_FLAGS "${BLT_CXX_FLAGS} -fopenmp" CACHE PATH "")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp" CACHE STRING "" FORCE)
   endif()
endif()

get_property(RAJA_INCLUDE_DIRS TARGET RAJA PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
include_directories(${RAJA_INCLUDE_DIRS})

find_library( RAJA_LIBRARY NAMES RAJA libRAJA
              PATHS ${RAJA_LIB_DIR} ${RAJA_DIR}/../../../lib/
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set RAJA_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(RAJA  DEFAULT_MSG
                                  RAJA_INCLUDE_DIRS
                                  RAJA_LIBRARY )
