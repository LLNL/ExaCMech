// -*-c++-*-

#ifndef ECMECH_port_h__
#define ECMECH_port_h__

#include "ECMech_gpu_portability.h"

#if ECMECH_HAVE_MSLIB

#include "MS_port.h"
#include "MS_Log.h"

#ifdef __ecmech_host_only__
#define ECMECH_FAIL(loc, str) MS_Fail(loc, str);
#define ECMECH_WARN(loc, str) MS_Warn(loc, str);
#else
#define ECMECH_FAIL(loc, str) MS_Fail(loc, str);
#define ECMECH_WARN(loc, str) MS_Warn(loc, str);
#endif

#else
// ECMECH_HAVE_MSLIB


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
