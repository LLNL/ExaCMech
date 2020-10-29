################################
# Version
################################
set(PACKAGE_BUGREPORT "barton22@llnl.gov")

set(ECMech_VERSION_MAJOR 0)
set(ECMech_VERSION_MINOR 3)
set(ECMech_VERSION_PATCH \"0\")

set(HEADER_INCLUDE_DIR
    ${PROJECT_BINARY_DIR}/include
    CACHE PATH
    "Directory where all generated headers will go in the build tree")

configure_file( src/ecmech/ECMech_Version.h.in
                ${HEADER_INCLUDE_DIR}/ECMech_Version.h )

################################
# Setup build options and their default values
################################
include(cmake/ECMechOptions.cmake)

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

        if (ENABLE_${_ucdep} OR ${_ucdep}_FOUND)
            list(APPEND ${arg_LIST_NAME} ${_dep})
        endif()
    endforeach()

    unset(_dep)
    unset(_ucdep)

endmacro(ecmech_fill_depends_list)
