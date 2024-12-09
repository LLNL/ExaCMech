#pragma once

#include "ECMech_config.h"
#include "ECMech_gpu_portability.h"

#if defined(UNUSED)
#elif defined(__GNUC__)
# define UNUSED(x) UNUSED_ ## x __attribute__((unused))
#else
# define UNUSED(x) x
#endif

#if defined(__ecmech_gpu_active__)
#define UNUSED_GPU(x) x
#else
#if defined(UNUSED)
# define UNUSED_GPU(x) UNUSED(x)
#elif defined(__GNUC__)
# define UNUSED_GPU(x) UNUSED_ ## x __attribute__((unused))
#else
# define UNUSED_GPU(x) x
#endif
#endif

#if defined(ECMECH_EXTRA_SOLVERS)
#define UNUSED_EXTRA(x) x
#else
#if defined(UNUSED)
# define UNUSED_EXTRA(x) UNUSED(x)
#elif defined(__GNUC__)
# define UNUSED_EXTRA(x) UNUSED_ ## x __attribute__((unused))
#else
# define UNUSED_EXTRA(x) x
#endif
#endif