#ifndef MINIAPP_UTIL_H
#define MINIAPP_UTIL_H

#include "RAJA/RAJA.hpp"

#if defined(RAJA_ENABLE_CUDA)
#include "RAJA/policy/cuda/raja_cudaerrchk.hpp"
#endif

//We're going to use this to determine what RAJA code to run for our
//kernels.
//The HIP backend won't be able to run on AMD GPGPUs
//until device function pointers are supported.
enum class Accelerator { CPU, CUDA, OPENMP };

//The below is a simple memory manager taken directly from the RAJA repo and as such
//the necessary copyright/LICENSE information is provided towards the bottom of this file
//for it.
//The memory manager is also sufficient for our basic needs for this miniapp. The CUDA calls
//make use of unified memory which might not be sufficient for our needs later on if we need to start
//experimenting with using pinned memory type models for GPU runs.

/*
  As RAJA does not manage memory we include a general purpose memory
  manager which may be used to perform c++ style allocation/deallocation
  or allocate/deallocate CUDA unified memory. The type of memory allocated
  is dependent on how RAJA was configured.
*/
namespace memoryManager
{
   template <typename T>
   T *allocate(RAJA::Index_type size)
   {
      T *ptr;
#if defined(RAJA_ENABLE_CUDA)
      cudaErrchk(
         cudaMallocManaged((void * *) &ptr, sizeof(T) * size, cudaMemAttachGlobal));
#else
      ptr = new T[size];
#endif
      return ptr;
   }

   template <typename T>
   void deallocate(T *&ptr)
   {
      if (ptr) {
#if defined(RAJA_ENABLE_CUDA)
         cudaErrchk(cudaFree(ptr));
#else
         delete[] ptr;
#endif
         ptr = nullptr;
      }
   }
};  // namespace memoryManager

// Copyright (c) 2016-19, Lawrence Livermore National Security, LLC.
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the disclaimer below.

// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the disclaimer (as noted below) in the
//   documentation and/or other materials provided with the distribution.

// * Neither the name of the LLNS/LLNL nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
// THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
// OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
// EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// Additional BSD Notice

// 1. This notice is required to be provided under our contract with the U.S.
// Department of Energy (DOE). This work was produced at Lawrence Livermore
// National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.

// 2. Neither the United States Government nor Lawrence Livermore National
// Security, LLC nor any of their employees, makes any warranty, express or
// implied, or assumes any liability or responsibility for the accuracy,
// completeness, or usefulness of any information, apparatus, product, or
// process disclosed, or represents that its use would not infringe
// privately-owned rights.

// 3. Also, reference herein to any specific commercial products, process,
// or services by trade name, trademark, manufacturer or otherwise does not
// necessarily constitute or imply its endorsement, recommendation, or favoring
// by the United States Government or Lawrence Livermore National Security, LLC.
// The views and opinions of authors expressed herein do not necessarily state
// or reflect those of the United States Government or Lawrence Livermore
// National Security, LLC, and shall not be used for advertising or product
// endorsement purposes.

#endif