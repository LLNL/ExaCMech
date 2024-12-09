#pragma once

#include <vector>

namespace ecmech
{
    class matModelBase;
}

// Here we're going to initialize all of the data that's going inside of
// of our material update function call.
// This function is used to initialize the data originally
void init_data(const std::vector<double>& ori_vec, const ecmech::matModelBase* mat_model_base,
               const int nqpts, const int num_hardness,
               const int num_slip, const int ind_gdot,
               const int state_var_vdim, double* state_vars);

// This sets the macroscopic velocity_grad to be purely deviatoric and behaving as a tension test in the
// z direction. More interesting velocity_grads could be created just as easily as well where we also have some
// spin terms as well. We could also create a case where there is some sort of spin term as well.
void setup_velocity_grad(const std::vector<double>& velocity_grad_input, double* const velocity_grad, const int nqpts);

// This function/kernel is used to set-up the problem at each time step
void setup_data(const int nqpts, const int nstatev,
                const double dt, const double* vel_grad_array,
                const double* cauchy_stress_array, const double* state_vars_array,
                double* cauchy_stress_d6p_array, double* def_rate_d6v_array,
                double* spin_vec_array, double* ddsdde_array,
                double* rel_vol_ratios_array, double* internal_energy_array,
                double* tkelv_array);
