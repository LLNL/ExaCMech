set(ECMECH_TEST_DEPENDS )

SET(CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}   -DDEBUG")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

if(BUILD_STATIC_LIBS)
    list(APPEND ECMECH_MINIAPP_DEPENDS ecmech_static ${SNLS_STATIC_TARGET})
else()
    list(APPEND ECMECH_MINIAPP_DEPENDS ecmech_shared ${SNLS_SHARED_TARGET})
endif()

if(ENABLE_OPENMP)
   list(APPEND ECMECH_MINIAPP_DEPENDS openmp)
endif()

if(ENABLE_CUDA)
   list(APPEND ECMECH_MINIAPP_DEPENDS cuda)
endif()

set(MINIAPP_SOURCES
    orientation_evolution.cpp
    material_kernels.cpp
    retrieve_kernels.cpp
    setup_kernels.cpp
    )

set(MINIAPP_HEADERS
    miniapp_util.h
    material_kernels.h
    retrieve_kernels.h
    setup_kernels.h
    )



blt_add_executable(NAME       orientation_evolution
                   HEADERS    ${MINIAPP_HEADERS}
                   SOURCES    ${MINIAPP_SOURCES}
                   OUTPUT_DIR ${MINIAPP_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_MINIAPP_DEPENDS})
