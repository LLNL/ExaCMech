#ifndef MATERIAL_KERNELS_H
#define MATERIAL_KERNELS_H

#include "ECMech_cases.h"
#include "ECMech_evptnWrap.h"
#include "RAJA/RAJA.hpp"
#include "miniapp_util.h"

#include <math.h>

using namespace ecmech;

void mat_model_kernel(const ecmech::matModelBase* mat_model_base,
                      const int nqpts, const double dt,
                      const int nstatev, double* state_vars_array,
                      double* stress_svec_p_array, double* d_svec_p_array,
                      double* w_vec_array, double* ddsdde_array,
                      double* vol_ratio_array, double* eng_int_array);

#endif