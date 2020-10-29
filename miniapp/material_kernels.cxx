#include "material_kernels.h"

using namespace ecmech;

// All of the parallelization operations are within the getResponse function of the material class.
void mat_model_kernel(const ecmech::matModelBase* mat_model_base,
                      const int nqpts, const double dt, double* state_vars_array,
                      double* stress_svec_p_array, double* d_svec_p_array,
                      double* w_vec_array, double* ddsdde_array,
                      double* vol_ratio_array, double* eng_int_array,
                      double* temp_array, double* sdd_array){
   mat_model_base->getResponse(dt, d_svec_p_array, w_vec_array, vol_ratio_array,
                               eng_int_array, stress_svec_p_array, state_vars_array,
                               temp_array, sdd_array, ddsdde_array, nqpts);
} // end of mat_model_kernel

