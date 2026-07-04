/**
 * @file ECMech_gpu_portability.h
 * @brief Compiler-portability macros for host/device function annotations across CUDA,
 * HIP, and CPU-only builds.
 *
 * ExaCMech's kernels are written once and compiled either for the CPU or for a GPU back
 * end (CUDA or HIP) selected at build time. Rather than sprinkling `#ifdef __CUDACC__`
 * throughout the codebase, every function that needs a `__host__`/`__device__`
 * annotation uses one of the `__ecmech_*__` macros defined here, which expand to the
 * appropriate CUDA/HIP annotation when compiling for a GPU back end and to nothing when
 * compiling for the CPU.
 *
 * @see ECMech_core.h for the umbrella header that pulls this in
 */

#pragma once

#ifdef __CUDACC__
#include <cuda_runtime_api.h>
#endif

#if defined(__HIPCC__)
#include <hip/hip_runtime.h>
#endif

/**
 * @def __ecmech_gpu_active__
 * @brief Defined (with no value) when compiling for a CUDA or HIP GPU target; can be
 * used to conditionally compile GPU-only code paths. Not defined for CPU-only builds.
 */
/**
 * @def __ecmech_host__
 * @brief Marks a function as callable from host code only. Expands to `__host__` when
 * compiling for CUDA/HIP, or to nothing for CPU-only builds.
 */
/**
 * @def __ecmech_device__
 * @brief Marks a function as callable from device (GPU) code only. Expands to
 * `__device__` when compiling for CUDA/HIP, or to nothing for CPU-only builds.
 */
/**
 * @def __ecmech_global__
 * @brief Marks a function as a GPU kernel entry point. Expands to `__global__` when
 * compiling for CUDA/HIP, or to nothing for CPU-only builds.
 */
/**
 * @def __ecmech_hdev__
 * @brief Marks a function as callable from both host and device code — the annotation
 * used by the vast majority of ExaCMech's math/kinetics routines so the same
 * implementation runs unmodified on CPU or GPU. Expands to `__host__ __device__` when
 * compiling for CUDA/HIP, or to nothing for CPU-only builds.
 */
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

/**
 * @def ECMECH_GPU_THREADS
 * @brief Default number of threads per block used when launching ExaCMech's GPU
 * kernels. Adjust here if a non-standard block size is needed for a given target
 * architecture.
 */
#define ECMECH_GPU_THREADS 256

/**
 * @def __ecmech_device_only__
 * @brief Defined when the current compilation pass is generating device (GPU) code.
 * Use to guard code that can only run/compile on the device.
 *
 * nvcc/hipcc compile each `__host__ __device__` function twice — once for the host and
 * once for the device — so this and #__ecmech_host_only__ let code differentiate
 * between those two passes via `__CUDA_ARCH__` (only defined, and > 0, during the
 * device compilation pass) or `__HIP_DEVICE_COMPILE__`.
 */
/**
 * @def __ecmech_host_only__
 * @brief Defined when the current compilation pass is generating host (CPU) code only —
 * i.e. `__ecmech_device_only__` is not defined. Use to guard code (such as
 * exception-based error handling) that is only valid on the host.
 */
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))  || defined(__HIP_DEVICE_COMPILE__)
#define __ecmech_device_only__
#else
#define __ecmech_host_only__
#endif
