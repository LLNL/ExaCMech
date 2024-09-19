#pragma once

namespace ecmech
{
    class matModelBase;
}

// This will then be the final function/kernel to save off all the data at
// each time step.
void retrieve_data(const int nqpts, const int nstatev,
                   const double* stress_svec_p_array, const double* vol_ratio_array,
                   const double* eng_int_array, double* state_vars_array,
                   double* stress_array);

