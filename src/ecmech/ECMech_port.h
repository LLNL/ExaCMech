/**
 * @file ECMech_port.h
 * @brief Portable error-reporting macros used in place of directly throwing exceptions
 * or calling an external logging library.
 *
 * GPU device code cannot throw C++ exceptions, so ExaCMech routines that need to report
 * a fatal error or a warning do so through the `ECMECH_FAIL`/`ECMECH_WARN` macros defined
 * here rather than calling `throw` or a logger directly. This lets the same calling code
 * be compiled for either host or device:
 *
 * - On the host, `ECMECH_FAIL` throws a `std::runtime_error` (or calls into LLNL's MSLib
 *   logging, if `ECMECH_HAVE_MSLIB` is enabled) so the failure can be caught and handled
 *   by the calling application.
 * - On the device, exceptions are not available, so `ECMECH_FAIL` instead prints an
 *   error message; the caller is still responsible for checking return codes/status
 *   flags to detect the failure (see, e.g., the `nFEvals < 0` convention used by
 *   `updateH1`/`updateHN` in ECMech_kinetics.h).
 *
 * Both macros take a `loc` string (conventionally `__func__`) identifying where the
 * failure/warning occurred, and a `str` message describing it.
 *
 * @see ECMech_gpu_portability.h for the `__ecmech_host_only__` macro used to select
 * between the throwing and non-throwing implementations
 */

// -*-c++-*-

#ifndef ECMECH_port_h__
#define ECMECH_port_h__

#include "ECMech_gpu_portability.h"

#if ECMECH_HAVE_MSLIB

#include "MS_port.h"
#include "MS_Log.h"

/**
 * @def ECMECH_FAIL(loc, str)
 * @brief Report a fatal error. When `ECMECH_HAVE_MSLIB` is enabled, defers to LLNL's
 * MSLib `MS_Fail` on both host and device.
 * @param loc Location string identifying the call site (typically `__func__`).
 * @param str Message describing the failure.
 */
/**
 * @def ECMECH_WARN(loc, str)
 * @brief Report a non-fatal warning via MSLib's `MS_Warn`.
 * @param loc Location string identifying the call site (typically `__func__`).
 * @param str Message describing the warning.
 */
#ifdef __ecmech_host_only__
#define ECMECH_FAIL(loc, str) MS_Fail(loc, str);
#define ECMECH_WARN(loc, str) MS_Warn(loc, str);
#else
#define ECMECH_FAIL(loc, str) MS_Fail(loc, str);
#define ECMECH_WARN(loc, str) MS_Warn(loc, str);
#endif

#else
// ECMECH_HAVE_MSLIB


/**
 * @def ECMECH_FAIL(loc, str)
 * @brief Report a fatal error.
 *
 * On the host, throws a `std::runtime_error` combining `loc` and `str` so the caller can
 * catch it. On the device (where exceptions are unavailable), prints the error message
 * instead; the caller must detect the failure through some other means (a negative
 * return code, etc.).
 * @param loc Location string identifying the call site (typically `__func__`).
 * @param str Message describing the failure.
 */
/**
 * @def ECMECH_WARN(loc, str)
 * @brief Report a non-fatal warning by printing a formatted message to stdout, on both
 * host and device.
 * @param loc Location string identifying the call site (typically `__func__`).
 * @param str Message describing the warning.
 */
#ifdef __ecmech_host_only__
#include <stdio.h>
#include <exception>
#include <stdexcept>
#define ECMECH_FAIL(loc, str) throw std::runtime_error(std::string("at ") + std::string(loc) + std::string( \
                                                          " failure : ") + std::string(str));
#define ECMECH_WARN(loc, str) printf("WARNING : ECMECH warning in %s : %s\n", loc, str);
#else
#define ECMECH_FAIL(loc, str) printf("ERROR : ECMECH failure in %s : %s\n", loc, str);
#define ECMECH_WARN(loc, str) printf("WARNING : ECMECH warning in %s : %s\n", loc, str);
#endif

#endif
// ECMECH_HAVE_MSLIB

#endif
// ECMECH_port_h__
