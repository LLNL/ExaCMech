#include(${PYBIND_SOURCE_DIR}/CMakeLists.txt)

set(PYECMECH_DEPENDS )

ecmech_fill_depends_list(LIST_NAME  PYECMECH_DEPENDS
                        DEPENDS_ON  RAJA snls)

if (${BLT_VERSION} VERSION_GREATER_EQUAL 0.6.0)
    if(ENABLE_CUDA)
        list(APPEND PYECMECH_DEPENDS blt::cuda_runtime blt::cuda) 
    endif()
    if(ENABLE_OPENMP)
        list(APPEND PYECMECH_DEPENDS blt::openmp) 
    endif()
else()
    if(ENABLE_CUDA)
        list(APPEND PYECMECH_DEPENDS cuda cuda_runtime) 
    endif()
    if(ENABLE_OPENMP)
        list(APPEND PYECMECH_DEPENDS openmp) 
    endif()
endif()

if(ENABLE_HIP)
    list(APPEND PYECMECH_DEPENDS blt::hip blt::hip_runtime)
endif()

if(NOT USE_BUILT_SNLS)
    list(APPEND PYECMECH_DEPENDS snls)
endif()

list(APPEND PYECMECH_DEPENDS pybind11::pybind11 ecmech_static)

message("-- PYECMECH_DEPENDS: ${PYECMECH_DEPENDS}")

blt_add_library(NAME        pyecmech_lib_static
                OUTPUT_NAME pyecmech_lib
                DEFINES     ${ECMECH_DEFINES}
                DEPENDS_ON  ${PYECMECH_DEPENDS}
                HEADERS     ecmechpy.hpp
                SOURCES     ecmechpy.cpp
                SHARED      FALSE)

install(TARGETS pyecmech_lib_static DESTINATION python)

pybind11_add_module(pyecmech ./ecmech_pybind11.cpp)
target_include_directories(pyecmech SYSTEM PUBLIC ${ECMECH_HEADERS})

target_link_libraries(pyecmech PUBLIC ecmech_static pyecmech_lib_static PRIVATE )

if(ENABLE_OPENMP)
    target_link_libraries(pyecmech PRIVATE OpenMP::OpenMP_CXX)
endif()

install(TARGETS pyecmech DESTINATION python)