#pragma once

namespace ecmech
{
    class matModelBase;
}

// This will then be the final function/kernel to save off all the data at
// each time step.
void retrieve_data(const int nqpts, const int nstatev,
                   const double* cauchy_stress_d6p_array, const double* rel_vol_ratios_array,
                   const double* internal_energy_array, double* state_vars_array,
                   double* cauchy_stress_array);

