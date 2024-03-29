set(ECMECH_MINIAPP_DEPENDS )

ecmech_fill_depends_list(LIST_NAME  ECMECH_MINIAPP_DEPENDS
                        DEPENDS_ON raja snls cuda openmp blt::hip blt::hip_runtime)

if(BUILD_STATIC_LIBS)
    list(APPEND ECMECH_MINIAPP_DEPENDS ecmech_static ${SNLS_STATIC_TARGET})
else()
    list(APPEND ECMECH_MINIAPP_DEPENDS ecmech_shared ${SNLS_SHARED_TARGET})
endif()

set(MINIAPP_SOURCES
    orientation_evolution.cxx
    material_kernels.cxx
    retrieve_kernels.cxx
    setup_kernels.cxx
    )

set(MINIAPP_HEADERS
    miniapp_util.h
    material_kernels.h
    retrieve_kernels.h
    setup_kernels.h
    )



set(ECMECH_MINIAPP_EXEC orientation_evolution)

blt_add_executable(NAME       ${ECMECH_MINIAPP_EXEC}
                   HEADERS    ${MINIAPP_HEADERS}
                   SOURCES    ${MINIAPP_SOURCES}
                   OUTPUT_DIR ${MINIAPP_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_MINIAPP_DEPENDS})

#------------------------------------------------------------------------------
# Install files
#------------------------------------------------------------------------------

install(TARGETS ${ECMECH_MINIAPP_EXEC} DESTINATION bin)
