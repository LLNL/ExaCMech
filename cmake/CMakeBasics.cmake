################################
# Version
################################
set(PACKAGE_BUGREPORT "carson16@llnl.gov")

set(ECMech_VERSION_MAJOR 0)
set(ECMech_VERSION_MINOR 4)
set(ECMech_VERSION_PATCH \"3\")

set(ECMECH_HEADER_INCLUDE_DIR
    ${PROJECT_BINARY_DIR}/include/ecmech
    CACHE PATH
    "Directory where all generated headers will go in the build tree")

################################
# Setup build options and their default values
################################

include(cmake/ECMechOptions.cmake)

##############################
# settings into ECMech_config.h
##############################

set(HAVE_ECMECH "1" CACHE STRING "")

if(USE_DPEFF)
    set(ECMECH_USE_DPEFF "1" CACHE STRING "")
endif()

if(ENABLE_PYTHON)
    set(ECMECH_PY "1" CACHE STRING "")
endif()

if(ENABLE_EXTRA_SOLVERS)
    set(ECMECH_EXTRA_SOLVERS "1" CACHE STRING "1")
endif()

if(CMAKE_BUILD_TYPE MATCHES DEBUG)
    set(ECMECH_DEBUG "1" CACHE STRING "")
endif()

configure_file( src/ecmech/ECMech_config.h.in
                ${ECMECH_HEADER_INCLUDE_DIR}/ECMech_config.h )

################################
# Third party library setup
################################
include(cmake/thirdpartylibraries/SetupThirdPartyLibraries.cmake)

##------------------------------------------------------------------------------
## ecmech_fill_depends_on_list
##
## This macro adds a dependency to the list if ENABLE_<dep name> or <dep name>_FOUND
##------------------------------------------------------------------------------
macro(ecmech_fill_depends_list)

    set(options )
    set(singleValueArgs LIST_NAME)
    set(multiValueArgs DEPENDS_ON)

    # Parse the arguments to the macro
    cmake_parse_arguments(arg
        "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

    if (NOT arg_LIST_NAME)
        message(FATAL_ERROR "ecmech_fill_depends_list requires argument LIST_NAME")
    endif()

    if (NOT arg_DEPENDS_ON)
        message(FATAL_ERROR "ecmech_fill_depends_list requires argument DEPENDS_ON")
    endif()

    foreach( _dep ${arg_DEPENDS_ON})
        string(TOUPPER ${_dep} _ucdep)

        if (ENABLE_${_ucdep} OR ${_ucdep}_FOUND OR ${_dep}_FOUND)
            list(APPEND ${arg_LIST_NAME} ${_dep})
        endif()
    endforeach()

    unset(_dep)
    unset(_ucdep)

endmacro(ecmech_fill_depends_list)
