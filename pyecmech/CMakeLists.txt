#include(${PYBIND_SOURCE_DIR}/CMakeLists.txt)

pybind11_add_module(ecmech ./ecmech_pybind11.cpp ./ecmechpy.cpp)
target_include_directories(ecmech SYSTEM PUBLIC ${ECMECH_HEADERS})

target_link_libraries(ecmech PUBLIC ecmech_static PRIVATE )

if(ENABLE_OPENMP)
    target_link_libraries(ecmech PRIVATE OpenMP::OpenMP_CXX)
endif()