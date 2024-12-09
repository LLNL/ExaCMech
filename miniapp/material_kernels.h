#pragma once

namespace ecmech
{
    class matModelBase;
}

void mat_model_kernel(const ecmech::matModelBase* mat_model_base,
                      const int nqpts, const double dt, double* state_vars_array,
                      double* cauchy_stress_d6p_array, double* def_rate_d6v_array,
                      double* spin_vec_array, double* ddsdde_array,
                      double* rel_vol_ratios_array, double* internal_energy_array,
                      double* tkelv_array, double* sdd_array);

