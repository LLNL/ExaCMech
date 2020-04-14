set(ECMECH_TEST_DEPENDS )

SET(CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}   -DDEBUG")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

if(BUILD_STATIC_LIBS)
    list(APPEND ECMECH_TEST_DEPENDS ecmech_static ${SNLS_STATIC_TARGET})
else()
    list(APPEND ECMECH_TEST_DEPENDS ecmech_shared ${SNLS_SHARED_TARGET})
endif()

if(ENABLE_OPENMP)
   list(APPEND ECMECH_TEST_DEPENDS openmp)
endif()

if(ENABLE_CUDA)
   list(APPEND ECMECH_TEST_DEPENDS cuda)
endif()

blt_add_executable(NAME       test_hardening
                   SOURCES    test_hardening.cxx
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_TEST_DEPENDS} gtest)
blt_add_test(NAME    test_hardening
             COMMAND test_hardening)

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

blt_add_executable(NAME       test_evptn_hcp
                   SOURCES    test_evptn_hcp.cxx
                   OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                   DEPENDS_ON ${ECMECH_TEST_DEPENDS} gtest)
blt_add_test(NAME    test_evptn_hcp
             COMMAND test_evptn_hcp)

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

