/**
 * @file material_kernels.cxx
 * @brief Implementation of `mat_model_kernel` (declared in `material_kernels.h`).
 */

#include "material_kernels.h"

#include "ECMech_evptnWrap.h"

// See material_kernels.h for the full @brief/@param doc. All of the parallelization
// (across nqpts material points, per the model's configured ExecutionStrategy) happens
// inside getResponseECM itself -- this function only forwards the arguments.
void mat_model_kernel(const ecmech::matModelBase* mat_model_base,
                      const int nqpts, const double dt, double* state_vars_array,
                      double* cauchy_stress_d6p_array, double* def_rate_d6v_array,
                      double* spin_vec_array, double* ddsdde_array,
                      double* rel_vol_ratios_array, double* internal_energy_array,
                      double* tkelv_array, double* sdd_array){
   mat_model_base->getResponseECM(dt, def_rate_d6v_array, spin_vec_array, rel_vol_ratios_array,
                                  internal_energy_array, cauchy_stress_d6p_array, state_vars_array,
                                  tkelv_array, sdd_array, ddsdde_array, nqpts);
} // end of mat_model_kernel

