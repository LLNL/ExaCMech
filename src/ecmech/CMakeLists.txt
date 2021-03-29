#------------------------------------------------------------------------------
# Sources
#------------------------------------------------------------------------------
set(ECMECH_HEADERS
    ${HEADER_INCLUDE_DIR}/ECMech_config.h
    ECMech_cases.h
    ECMech_const.h
    ECMech_core.h
    ECMech_cuda_portability.h
    ECMech_eosSimple.h
    ECMech_evptn.h
    ECMech_kinetics.h
    ECMech_kinetics_KMBalD.h
    ECMech_kinetics_VocePL.h
    ECMech_port.h
    ECMech_slipgeom.h
    ECMech_util.h
    ECMech_matModelBase.h
    ECMech_evptnWrap.h
    mc_vars_set.h
    vad_vars_set.h
    vadl_vars_set.h
    vb_d_vars_set.h
    vw_vars_set.h)

set(ECMECH_SOURCES
    ECMech_kinetics.cxx
    ECMech_cases.cxx)
#    ECMech_slipgeom.cxx)

#------------------------------------------------------------------------------
# Dependencies
#------------------------------------------------------------------------------
set(ECMECH_DEPENDS)

ecmech_fill_depends_list(LIST_NAME  ECMECH_DEPENDS
                        DEPENDS_ON raja snls cub cuda openmp hip rocPRIM)

message("-- ECMECH_DEPENDS: ${ECMECH_DEPENDS}")

#------------------------------------------------------------------------------
# Defines
#------------------------------------------------------------------------------
set(ECMECH_DEFINES HAVE_ECMECH)

SET(CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}   -DECMECH_DEBUG")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DECMECH_DEBUG")

#------------------------------------------------------------------------------
# Includes
#------------------------------------------------------------------------------
set(ECMECH_INCLUDES ${CMAKE_SOURCE_DIR}/src/interface
                    ${HEADER_INCLUDE_DIR}
                    ${CMAKE_SOURCE_DIR}/src/ecmech)

if(NOT USE_BUILT_SNLS)
    list(APPEND ECMECH_INCLUDES ${CMAKE_SOURCE_DIR}/snls)
endif()

#------------------------------------------------------------------------------
# Build Targets
#------------------------------------------------------------------------------
if(BUILD_SHARED_LIBS)

    blt_add_library(NAME        ecmech_shared
                    OUTPUT_NAME ecmech
                    DEFINES     ${ECMECH_DEFINES}
                    INCLUDES    ${ECMECH_INCLUDES}
                    DEPENDS_ON  ${ECMECH_DEPENDS} # ${SNLS_SHARED_TARGET}
                    HEADERS     ${ECMECH_HEADERS}
                    SOURCES     ${ECMECH_SOURCES}
                    SHARED      TRUE)
endif(BUILD_SHARED_LIBS)

if(BUILD_STATIC_LIBS)

    blt_add_library(NAME        ecmech_static
                    OUTPUT_NAME ecmech
                    DEFINES     ${ECMECH_DEFINES}
                    INCLUDES    ${ECMECH_INCLUDES}
                    DEPENDS_ON  ${ECMECH_DEPENDS} # ${SNLS_STATIC_TARGET}
                    HEADERS     ${ECMECH_HEADERS}
                    SOURCES     ${ECMECH_SOURCES}
                    SHARED      FALSE)
endif(BUILD_STATIC_LIBS)

#------------------------------------------------------------------------------
# Install files
#------------------------------------------------------------------------------

set(ECMECH_INSTALLED_HEADERS
    ${ECMECH_HEADERS})

install(FILES ${ECMECH_INSTALLED_HEADERS} DESTINATION include)

if(BUILD_SHARED_LIBS)
    install(TARGETS ecmech_shared DESTINATION lib)
endif(BUILD_SHARED_LIBS)

if(BUILD_STATIC_LIBS)
    install(TARGETS ecmech_static DESTINATION lib)
endif(BUILD_STATIC_LIBS)
