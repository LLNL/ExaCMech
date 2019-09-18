#ifndef SETUP_KERNELS_H
#define SETUP_KERNELS_H

#include "ECMech_cases.h"
#include "ECMech_evptnWrap.h"
#include "RAJA/RAJA.hpp"
#include "miniapp_util.h"

#include <math.h>

using namespace ecmech;

//Here we're going to initialize all of the data that's going inside of
//of our material update function call.
//This function is used to initialize the data originally
void init_data(Accelerator accel, const double* ori, const ecmech::matModelBase* mat_model_base,
               const int nqpts, const int num_hardness,
               const int num_slip, const int ind_gdot,
               const int state_var_vdim, double* state_vars);

//This sets the macroscopic vgrad to be purely deviatoric and behaving as a tension test in the
//z direction. More interesting vgrads could be created just as easily as well where we also have some
//spin terms as well. We could also create a case where there is some sort of spin term as well.
void setup_vgrad(double* vgrad, const int nqpts);

//This function/kernel is used to set-up the problem at each time step
void setup_data(Accelerator accel, const int nqpts, const int nstatev,
                const double dt, const double* vel_grad_array,
                const double* stress_array, const double* state_vars_array,
                double* stress_svec_p_array, double* d_svec_p_array,
                double* w_vec_array, double* ddsdde_array,
                double* vol_ratio_array, double* eng_int_array);

#endif