#pragma once

namespace ecmech
{
    class matModelBase;
}

void mat_model_kernel(const ecmech::matModelBase* mat_model_base,
                      const int nqpts, const double dt, double* state_vars_array,
                      double* stress_svec_p_array, double* d_svec_p_array,
                      double* w_vec_array, double* ddsdde_array,
                      double* vol_ratio_array, double* eng_int_array,
                      double* temp_array, double* sdd_array);

