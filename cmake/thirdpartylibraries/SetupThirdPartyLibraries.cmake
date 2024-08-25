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
# RAJA
################################

if (RAJA_DIR)
   find_package(RAJA REQUIRED CONFIG PATHS ${RAJA_DIR})
else()
   message(FATAL_ERROR "RAJA_DIR was not provided. It is needed to find RAJA.")
endif()


################################
# SNLS
################################

if (SNLS_DIR)
    find_package(snls REQUIRED CONFIG PATHS ${SNLS_DIR})
    set_target_properties(snls PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${SNLS_INCLUDE_DIRS}")
endif()


