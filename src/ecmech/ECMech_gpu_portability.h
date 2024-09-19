#pragma once

#ifdef __CUDACC__
#include <cuda_runtime_api.h>
#endif

#if defined(__HIPCC__)
#include <hip/hip_runtime.h>
#endif

// When compiling using the Nvidia/CUDA tools, nvcc defines the host, device, and global
// labels to identify the compilation target for a particular module. Routines that
// are intended for the host need to be declared with __host__.  Similarly, routines
// that are intended for the GPU need to be declared using __device__. Routines
// that are intended for both host and GPU need to be declared using both __host__ and
// __device__.
//
// For non-CUDA builds, we need to declare empty macros for portability.
// ----------------------------------------------------------------------------------------

#if defined(__CUDACC__) || defined(__HIPCC__) 
#define __ecmech_gpu_active__
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

// Modify our number of threads as needed if we need to set it to something non-standard
#define ECMECH_GPU_THREADS 256

// __CUDA_ARCH__ is defined when compiling for the device, the macro below is used
// to filter code that cannot be compiled for the device.
// ----------------------------------------------------------------------------------------

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))  || defined(__HIP_DEVICE_COMPILE__)
#define __ecmech_device_only__
#else
#define __ecmech_host_only__
#endif
