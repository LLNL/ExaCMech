# Provide backwards compatibility for *_PREFIX options
set(_tpls 
    raja
    snls)

foreach(_tpl ${_tpls})
    string(TOUPPER ${_tpl} _uctpl)
    if (${_uctpl}_PREFIX)
        set(${_uctpl}_DIR ${${_uctpl}_PREFIX} CACHE PATH "")
        mark_as_advanced(${_uctpl}_PREFIX)
    endif()
endforeach()

################################
# CUB (required for CUDA build)
################################
if (ENABLE_CUDA AND CUB_DIR)
   include(cmake/thirdpartylibraries/FindCUB.cmake)
endif ()

################################
# RAJA
################################

if (RAJA_DIR)
    include(cmake/thirdpartylibraries/FindRAJA.cmake)
    if (RAJA_FOUND)
        blt_register_library( NAME       raja
                              TREAT_INCLUDES_AS_SYSTEM ON
                              INCLUDES   ${RAJA_INCLUDE_DIRS}
                              LIBRARIES  ${RAJA_LIBRARY})
    else()
        message(FATAL_ERROR "Unable to find RAJA with given path ${RAJA_DIR}")
    endif()
endif()

################################
# SNLS
################################

if (SNLS_DIR)
    include(cmake/thirdpartylibraries/FindSNLS.cmake)
    if (SNLS_FOUND)
        blt_register_library( NAME       snls
                              TREAT_INCLUDES_AS_SYSTEM ON
                              INCLUDES   ${SNLS_INCLUDE_DIRS}
                              LIBRARIES  ${SNLS_LIBRARIES}
                              DEPENDS_ON ${SNLS_DEPENDS})
    else()
        message(FATAL_ERROR "Unable to find SNLS with given path ${SNLS_DIR}")
    endif()
endif()


