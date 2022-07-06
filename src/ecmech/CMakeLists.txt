#------------------------------------------------------------------------------
# Sources
#------------------------------------------------------------------------------
set(ECMECH_HEADERS
    ${ECMECH_HEADER_INCLUDE_DIR}/ECMech_config.h
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
                        DEPENDS_ON raja snls cub cuda openmp blt::hip blt::hip_runtime)

if(NOT USE_BUILT_SNLS)
    list(APPEND ECMECH_DEPENDS snls)
endif()

message("-- ECMECH_DEPENDS: ${ECMECH_DEPENDS}")

#------------------------------------------------------------------------------
# Defines
#------------------------------------------------------------------------------
set(ECMECH_DEFINES)

#------------------------------------------------------------------------------
# Includes
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Build Targets
#------------------------------------------------------------------------------
if(BUILD_SHARED_LIBS)

    blt_add_library(NAME        ecmech_shared
                    OUTPUT_NAME ecmech
                    DEFINES     ${ECMECH_DEFINES}
                    DEPENDS_ON  ${ECMECH_DEPENDS}
                    HEADERS     ${ECMECH_HEADERS}
                    SOURCES     ${ECMECH_SOURCES}
                    SHARED      TRUE)
endif(BUILD_SHARED_LIBS)

if(BUILD_STATIC_LIBS)

    blt_add_library(NAME        ecmech_static
                    OUTPUT_NAME ecmech
                    DEFINES     ${ECMECH_DEFINES}
                    DEPENDS_ON  ${ECMECH_DEPENDS}
                    HEADERS     ${ECMECH_HEADERS}
                    SOURCES     ${ECMECH_SOURCES}
                    SHARED      FALSE)
endif(BUILD_STATIC_LIBS)

#------------------------------------------------------------------------------
# Install files
#------------------------------------------------------------------------------

install(FILES ${ECMECH_HEADERS} DESTINATION include/ecmech)


if(BUILD_SHARED_LIBS)
   target_include_directories(
      ecmech_shared
      PUBLIC
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/ecmech/>
      $<BUILD_INTERFACE:${ECMECH_HEADER_INCLUDE_DIR}>
      $<INSTALL_INTERFACE:include>)
    install(TARGETS ecmech_shared DESTINATION lib)
endif(BUILD_SHARED_LIBS)

if(BUILD_STATIC_LIBS)
   target_include_directories(
      ecmech_static
      PUBLIC
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/ecmech/>
      $<BUILD_INTERFACE:${ECMECH_HEADER_INCLUDE_DIR}>
      $<INSTALL_INTERFACE:include>)
    install(TARGETS ecmech_static DESTINATION lib)
endif(BUILD_STATIC_LIBS)

install(
  FILES ${ECMECH_HEADER_INCLUDE_DIR}/ECMech_config.h
  DESTINATION include/ecmech/)

install(
  FILES ${ECMECH_HEADERS}
  DESTINATION include/ecmech/)
