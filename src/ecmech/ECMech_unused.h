/**
 * @file ECMech_unused.h
 * @brief Macros for silencing "unused parameter" compiler warnings that only apply
 * conditionally, based on the active build configuration.
 *
 * Some function parameters are only used when a particular build option is enabled
 * (e.g. GPU support, or the optional "extra solvers"); when that option is disabled, the
 * parameter goes unused and would otherwise trigger a compiler warning on GCC-family
 * compilers. Wrapping the parameter name in the appropriate macro below suppresses that
 * warning only in the configurations where it would actually fire, while leaving the
 * parameter name intact (and thus usable/readable) in configurations where it's needed.
 */

#pragma once

#include "ECMech_config.h"
#include "ECMech_gpu_portability.h"

/**
 * @def UNUSED(x)
 * @brief Marks parameter `x` as intentionally unused, suppressing GCC/Clang's
 * `-Wunused-parameter` warning for it. If a `UNUSED` macro is already defined by
 * something else that was included first, this leaves it untouched.
 */
#if defined(UNUSED)
#elif defined(__GNUC__)
# define UNUSED(x) UNUSED_ ## x __attribute__((unused))
#else
# define UNUSED(x) x
#endif

/**
 * @def UNUSED_GPU(x)
 * @brief Like #UNUSED, but only actually suppresses the warning when
 * `__ecmech_gpu_active__` is *not* defined — i.e. for a parameter that is unused on CPU
 * builds but becomes used once GPU support is compiled in.
 */
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

/**
 * @def UNUSED_EXTRA(x)
 * @brief Like #UNUSED, but only actually suppresses the warning when
 * `ECMECH_EXTRA_SOLVERS` is *not* defined — i.e. for a parameter that is unused unless
 * the optional extra solvers are compiled in.
 */
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
