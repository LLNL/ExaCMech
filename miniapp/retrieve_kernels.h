#ifndef RETRIEVE_KERNELS_H
#define RETRIEVE_KERNELS_H

#include "ECMech_cases.h"
#include "ECMech_evptnWrap.h"
#include "RAJA/RAJA.hpp"
#include "miniapp_util.h"

#include <math.h>

using namespace ecmech;

// This will then be the final function/kernel to save off all the data at
// each time step. We could also show this kernel being used to compute the
// internal force (divergence of the Cauchy stress) for each point. However,
// we might make one last kernel to accomplish that goal.
void retrieve_data(ecmech::Accelerator accel, const int nqpts, const int nstatev,
                   const double* stress_svec_p_array, const double* vol_ratio_array,
                   const double* eng_int_array, double* state_vars_array,
                   double* stress_array);


#endif

