# Provide backwards compatibility for *_PREFIX options
set(_tpls 
    raja
    chai
    umpire
    camp
    fmt
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
    find_package(SNLS REQUIRED CONFIG PATHS ${SNLS_DIR})
    set_target_properties(snls PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${SNLS_INCLUDE_DIRS}")
endif()

if(SNLS_USE_RAJA_PORT_SUITE)
################################
# camp
################################

    if (CAMP_DIR)
        find_package(camp REQUIRED CONFIG PATHS ${CAMP_DIR})
    else()
        message(FATAL_ERROR "CAMP_DIR was not provided. It is needed to find CAMP.")
    endif()

################################
# chai
################################

    if (CHAI_DIR)
        set(umpire_DIR ${UMPIRE_DIR})
        set(raja_DIR ${RAJA_DIR})
        set(fmt_DIR ${FMT_DIR})
        find_package(chai REQUIRED CONFIG PATHS ${CHAI_DIR})
    else()
        message(FATAL_ERROR "CHAI_DIR was not provided. It is needed to find CHAI.")
    endif()

################################
# fmt
################################

    if (FMT_DIR)
        find_package(fmt CONFIG PATHS ${FMT_DIR})
    else()
        message(WARNING "FMT_DIR was not provided. This is a requirement for camp as of v2024.02.0. Ignore this warning if using older versions of the RAJA Portability Suite")
    endif()
    

################################
# UMPIRE
################################

    if (DEFINED UMPIRE_DIR)
        find_package(umpire REQUIRED CONFIG PATHS ${UMPIRE_DIR})
    else()
        message(FATAL_ERROR "UMPIRE_DIR was not provided. It is needed to find UMPIRE.")
    endif()
endif() # End SNLS_USE_RAJA_PORT_SUITE check

