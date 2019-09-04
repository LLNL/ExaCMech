#pragma once

#ifndef __ECMECH_CUDA_PORTABILITY_H
#define __ECMECH_CUDA_PORTABILITY_H

#include <stdlib.h>

#ifdef __CUDACC__
#include <cuda_runtime_api.h>
#endif

// When compiling using the Nvidia/CUDA tools, nvcc defines the host, device, and global
// labels to identify the compilation target for a particular module. Routines that
// are intended for the host need to be declared with __host__.  Similarly, routines
// that are intended for the GPU need to be declared using __device__. Routines
// that are intended for both host and GPU need to be declared using both __host__ and
// __device__.
//
// For non-CUDA builds, we need to declare empty macros for portability.
//----------------------------------------------------------------------------------------

#ifdef __CUDACC__
#define __ecmech_host__   __host__
#define __ecmech_device__ __device__
#define __ecmech_global__ __global__
#define __ecmech_hdev__   __host__ __device__
#else
#define __ecmech_host__
#define __ecmech_device__
#define __ecmech_global__
#define __ecmech_hdev__
#endif

// declare a CUDA runtime error handler and macro...
//----------------------------------------------------------------------------------------

#ifdef __CUDACC__
extern void CUDART_Check(const cudaError_t err, const char *file, const char *func, const int ln);

#define CUDART_CHECK(err) CUDART_Check(err, __FILE__, __func__, __LINE__);
#else
#define CUDART_CHECK(err)
#endif

// __CUDA_ARCH__ is defined when compiling for the device, the macro below is used
// to filter code that cannot be compiled for the device.
//----------------------------------------------------------------------------------------

#ifdef __CUDA_ARCH__
#define __cuda_device_only__
#else
#define __cuda_host_only__
#endif

#endif  // __ECMECH_CUDA_PORTABILITY_H
