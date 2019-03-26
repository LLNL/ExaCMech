option(BUILD_STATIC_LIBS "Build static libraries" ON)
option(BUILD_SHARED_LIBS "Build shared libraries" ON)

option(ENABLE_THROW "Enable throwing of exceptions" ON)

option(ENABLE_TESTS "Enable tests" ON)

option(ENABLE_CUDA "Enable CUDA" OFF)

# Force atleast static if user turns off both
if(NOT BUILD_STATIC_LIBS AND NOT BUILD_SHARED_LIBS)
    message("Both static and shared libaries were disabled."
            "Building static libraries re-enabled.")
    set(BUILD_STATIC_LIBS ON CACHE BOOL "Build static libraries" FORCE)
endif(NOT BUILD_STATIC_LIBS AND NOT BUILD_SHARED_LIBS)
