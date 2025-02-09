#------------------------------------------------------------------------------
# Sources
#------------------------------------------------------------------------------
set(ECMECH_HEADERS
    ${ECMECH_HEADER_INCLUDE_DIR}/ECMech_config.h
    cases/ECMech_cases_bcc_defs.h
    cases/ECMech_cases_fcc_defs.h
    cases/ECMech_cases_hcp_defs.h
    cases/ECMech_cases_util.h
    evptn/ECMech_base_classes.h
    evptn/ECMech_base_fcns.h
    evptn/ECMech_evptnSngl.h
    evptn/ECMech_evptn.h
    kinetics/ECMech_kinetics.h
    kinetics/ECMech_kinetics_KMBalD.h
    kinetics/ECMech_kinetics_OrowanD.h
    kinetics/ECMech_kinetics_VocePL.h
    kinetics/ECMech_kinetics_BCCMD.h
    slipgeom/ECMech_slipgeom_base.h
    slipgeom/ECMech_slipgeom_fcc.h
    slipgeom/ECMech_slipgeom_bcc.h
    slipgeom/ECMech_slipgeom_hcp.h
    util/mc_vars_set.h
    util/vad_vars_set.h
    util/vadl_vars_set.h
    util/vb_d_vars_set.h
    util/vw_vars_set.h
    ECMech_cases.h
    ECMech_const.h
    ECMech_core.h
    ECMech_gpu_portability.h
    ECMech_eosSimple.h
    ECMech_elastic.h
    ECMech_port.h
    ECMech_slipgeom.h
    ECMech_util.h
    ECMech_unused.h
    ECMech_matModelBase.h
    ECMech_evptnWrap.h)

set(ECMECH_SOURCES
    cases/ECMech_cases.cxx
    cases/ECMech_cases_fcc.cxx
    cases/ECMech_cases_bcc.cxx
    cases/ECMech_cases_bcc_oro.cxx
    cases/ECMech_cases_bcc_oro_big.cxx
    cases/ECMech_cases_hcp.cxx)

#------------------------------------------------------------------------------
# Dependencies
#------------------------------------------------------------------------------
set(ECMECH_DEPENDS)

ecmech_fill_depends_list(LIST_NAME  ECMECH_DEPENDS
                        DEPENDS_ON RAJA camp snls)

if (${BLT_VERSION} VERSION_GREATER_EQUAL 0.6.0)
    if(ENABLE_CUDA)
        list(APPEND ECMECH_DEPENDS blt::cuda_runtime blt::cuda) 
    endif()
    if(ENABLE_OPENMP)
        list(APPEND ECMECH_DEPENDS blt::openmp) 
    endif()
else()
    if(ENABLE_CUDA)
        list(APPEND ECMECH_DEPENDS cuda cuda_runtime) 
    endif()
    if(ENABLE_OPENMP)
        list(APPEND ECMECH_DEPENDS openmp) 
    endif()
endif()

if(ENABLE_HIP)
    list(APPEND ECMECH_DEPENDS blt::hip blt::hip_runtime)
endif()

if(NOT USE_BUILT_SNLS)
    list(APPEND ECMECH_DEPENDS snls)
    if (SNLS_USE_RAJA_PORT_SUITE)
        list(APPEND ECMECH_DEPENDS chai umpire)
    endif()
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

install(DIRECTORY ${PROJECT_SOURCE_DIR}/src/ DESTINATION include/ FILES_MATCHING PATTERN *.h)
install(
  FILES ${ECMECH_HEADER_INCLUDE_DIR}/ECMech_config.h
  DESTINATION include/ecmech/)


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
