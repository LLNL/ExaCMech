set(ECMECH_TEST_DEPENDS )

ecmech_fill_depends_list(LIST_NAME  ECMECH_TEST_DEPENDS
                        DEPENDS_ON RAJA snls)

if(NOT USE_BUILT_SNLS)
    list(APPEND ECMECH_TEST_DEPENDS snls)
    if (SNLS_USE_RAJA_PORT_SUITE)
        list(APPEND ECMECH_TEST_DEPENDS chai umpire camp fmt::fmt)
    endif()
endif()

if (${BLT_VERSION} VERSION_GREATER_EQUAL 0.6.0)
    if(ENABLE_CUDA)
        list(APPEND ECMECH_TEST_DEPENDS blt::cuda_runtime blt::cuda) 
    endif()
    if(ENABLE_OPENMP)
        list(APPEND ECMECH_TEST_DEPENDS blt::openmp) 
    endif()
else()
    if(ENABLE_CUDA)
        list(APPEND ECMECH_TEST_DEPENDS cuda cuda_runtime) 
    endif()
    if(ENABLE_OPENMP)
        list(APPEND ECMECH_TEST_DEPENDS openmp) 
    endif()
endif()

if(ENABLE_HIP)
    list(APPEND ECMECH_TEST_DEPENDS blt::hip blt::hip_runtime)
endif()

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

blt_add_executable(NAME       test_orowan_hardening_0
                   SOURCES    test_orowan_hardening.cxx
                   DEFINES    LARGE_DD
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_TEST_DEPENDS} gtest)
blt_add_test(NAME    test_orowan_hardening_0
             COMMAND test_orowan_hardening_0)

blt_add_executable(NAME       test_orowan_hardening_1
                   SOURCES    test_orowan_hardening.cxx
                   DEFINES
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_TEST_DEPENDS} gtest)
blt_add_test(NAME    test_orowan_hardening_1
             COMMAND test_orowan_hardening_1)

blt_add_executable(NAME       test_orowan_kinetics_0
                   SOURCES    test_orowan_kinetics.cxx
                   DEFINES    LARGE_DD
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_TEST_DEPENDS} gtest)
blt_add_test(NAME    test_orowan_kinetics_0
             COMMAND test_orowan_kinetics_0)

blt_add_executable(NAME       test_orowan_px_0
                   SOURCES    test_orowan_px.cxx
                   DEFINES    KIN_BCC
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_TEST_DEPENDS} gtest)
blt_add_test(NAME    test_orowan_px_0
             COMMAND test_orowan_px_0)

blt_add_executable(NAME       test_orowan_px_0_ns
                   SOURCES    test_orowan_px.cxx
                   DEFINES    KIN_BCC_NS
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_TEST_DEPENDS} gtest)
blt_add_test(NAME    test_orowan_px_0_ns
             COMMAND test_orowan_px_0_ns)

blt_add_executable(NAME       test_orowan_px_1
                   SOURCES    test_orowan_px.cxx
                   DEFINES    KIN_FCC
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_TEST_DEPENDS} gtest)
blt_add_test(NAME    test_orowan_px_1
             COMMAND test_orowan_px_1)

blt_add_executable(NAME       test_orowan_kinetics_1
                   SOURCES    test_orowan_kinetics.cxx
                   DEFINES
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_TEST_DEPENDS} gtest)
blt_add_test(NAME    test_orowan_kinetics_1
             COMMAND test_orowan_kinetics_1)

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
                   DEFINES    AS_DRIVER=1 DO_FD_CHECK_MTAN=1 ECMECH_DEBUG
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

