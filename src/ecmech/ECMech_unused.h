#pragma once

#include "ECMech_config.h"

#if defined(UNUSED)
#elif defined(__GNUC__)
# define UNUSED(x) UNUSED_ ## x __attribute__((unused))
#else
# define UNUSED(x) x
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