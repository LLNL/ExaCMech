set(ECMECH_TEST_DEPENDS )

SET(CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}   -DDEBUG")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

if(BUILD_STATIC_LIBS)
    list(APPEND ECMECH_TEST_DEPENDS ecmech_static snls)
else()
    list(APPEND ECMECH_TEST_DEPENDS ecmech_shared snls_shared)
endif()

blt_add_executable(NAME       test_evptn
                   SOURCES    test_evptn.cxx
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_TEST_DEPENDS})
blt_add_test(NAME    test_evptn
             COMMAND test_evptn)

