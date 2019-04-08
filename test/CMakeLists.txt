set(ECMECH_TEST_DEPENDS )

SET(CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}   -DDEBUG")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

if(BUILD_STATIC_LIBS)
    list(APPEND ECMECH_TEST_DEPENDS ecmech_static snls)
else()
    list(APPEND ECMECH_TEST_DEPENDS ecmech_shared snls_shared)
endif()

blt_add_executable(NAME       test_hardening
                   SOURCES    test_hardening.cxx
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_TEST_DEPENDS})
blt_add_test(NAME    test_hardening
             COMMAND test_hardening)

blt_add_executable(NAME       test_evptn_0
                   SOURCES    test_evptn.cxx
                   DEFINES    KIN_TYPE=0
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_TEST_DEPENDS})
blt_add_test(NAME    test_evptn_0
             COMMAND test_evptn_0)

blt_add_executable(NAME       test_evptn_1
                   SOURCES    test_evptn.cxx
                   DEFINES    KIN_TYPE=1
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_TEST_DEPENDS})
blt_add_test(NAME    test_evptn_1
             COMMAND test_evptn_1)

blt_add_executable(NAME       test_updst
                   SOURCES    test_updst.cxx
                   DEFINES    AS_DRIVER=0
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_TEST_DEPENDS})
blt_add_test(NAME    test_updst
             COMMAND test_updst)

blt_add_executable(NAME       test_driver
                   SOURCES    test_updst.cxx
                   DEFINES    AS_DRIVER=1
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_TEST_DEPENDS})
blt_add_test(NAME    test_driver
             COMMAND test_driver)

blt_add_executable(NAME       test_mtan
                   SOURCES    test_updst.cxx
                   DEFINES    AS_DRIVER=1 DO_FD_CHECK_MTAN=1
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_TEST_DEPENDS})
blt_add_test(NAME    test_mtan
             COMMAND test_mtan)

