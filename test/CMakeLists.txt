set(ECMECH_TEST_DEPENDS )

ecmech_fill_depends_list(LIST_NAME  ECMECH_TEST_DEPENDS
                        DEPENDS_ON raja snls cuda openmp blt::hip blt::hip_runtime)

if(BUILD_STATIC_LIBS)
    list(APPEND ECMECH_TEST_DEPENDS ecmech_static ${SNLS_STATIC_TARGET})
else()
    list(APPEND ECMECH_TEST_DEPENDS ecmech_shared ${SNLS_SHARED_TARGET})
endif()

blt_add_executable(NAME       test_hardening_0
                   SOURCES    test_hardening.cxx
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_TEST_DEPENDS} gtest)
blt_add_test(NAME    test_hardening_0
             COMMAND test_hardening_0)

blt_add_executable(NAME       test_hardening_1
                   SOURCES    test_hardening.cxx
                   DEFINES    KIN_NONLINEAR
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_TEST_DEPENDS} gtest)
blt_add_test(NAME    test_hardening_1
             COMMAND test_hardening_1)

blt_add_executable(NAME       test_aniso_hardening_0
                   SOURCES    test_aniso_hardening.cxx
                   DEFINES    ANISO_HARDENING
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_TEST_DEPENDS} gtest)
blt_add_test(NAME    test_aniso_hardening_0
             COMMAND test_aniso_hardening_0)

blt_add_executable(NAME       test_aniso_hardening_1
                   SOURCES    test_aniso_hardening.cxx
                   DEFINES    KIN_NONLINEAR ANISO_HARDENING
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_TEST_DEPENDS} gtest)
blt_add_test(NAME    test_aniso_hardening_1
             COMMAND test_aniso_hardening_1)

blt_add_executable(NAME       test_evptn_0
                   SOURCES    test_evptn.cxx
                   DEFINES    KIN_TYPE=0
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_TEST_DEPENDS} gtest)
blt_add_test(NAME    test_evptn_0
             COMMAND test_evptn_0)

blt_add_executable(NAME       test_evptn_1
                   SOURCES    test_evptn.cxx
                   DEFINES    KIN_TYPE=1
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_TEST_DEPENDS} gtest)
blt_add_test(NAME    test_evptn_1
             COMMAND test_evptn_1)

blt_add_executable(NAME       test_evptn_2
                   SOURCES    test_evptn.cxx
                   DEFINES    KIN_TYPE=2
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_TEST_DEPENDS} gtest)
blt_add_test(NAME    test_evptn_2
             COMMAND test_evptn_2)

blt_add_executable(NAME       test_evptn_3
                   SOURCES    test_evptn.cxx
                   DEFINES    KIN_TYPE=3
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_TEST_DEPENDS} gtest)
blt_add_test(NAME    test_evptn_3
             COMMAND test_evptn_3)

blt_add_executable(NAME       test_updst
                   SOURCES    test_updst.cxx
                   DEFINES    AS_DRIVER=0
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_TEST_DEPENDS} gtest)
blt_add_test(NAME    test_updst
             COMMAND test_updst)

blt_add_executable(NAME       test_mtan
                   SOURCES    test_updst.cxx
                   DEFINES    AS_DRIVER=1 DO_FD_CHECK_MTAN=1
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_TEST_DEPENDS} gtest)
blt_add_test(NAME    test_mtan
             COMMAND test_mtan)

blt_add_executable(NAME       test_px
                   SOURCES    test_px.cxx
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_TEST_DEPENDS} gtest)
blt_add_test(NAME    test_px
             COMMAND test_px)

