#ifndef RETRIEVE_KERNELS_H
#define RETRIEVE_KERNELS_H

#include "ECMech_cases.h"
#include "ECMech_evptnWrap.h"
#include "RAJA/RAJA.hpp"
#include "miniapp_util.h"

#include <math.h>

using namespace ecmech;

// This will then be the final function/kernel to save off all the data at
// each time step.
void retrieve_data(ecmech::ExecutionStrategy accel, const int nqpts, const int nstatev,
                   const double* stress_svec_p_array, const double* vol_ratio_array,
                   const double* eng_int_array, double* state_vars_array,
                   double* stress_array);


#endif

