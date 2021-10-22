#include(${PYBIND_SOURCE_DIR}/CMakeLists.txt)

pybind11_add_module(pyecmech ./ecmech_pybind11.cpp ./ecmechpy.cpp)
target_include_directories(pyecmech SYSTEM PUBLIC ${ECMECH_HEADERS})

target_link_libraries(pyecmech PUBLIC ecmech_static PRIVATE )

if(ENABLE_OPENMP)
    target_link_libraries(pyecmech PRIVATE OpenMP::OpenMP_CXX)
endif()

install(TARGETS pyecmech DESTINATION python)