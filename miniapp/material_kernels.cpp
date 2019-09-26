#include "material_kernels.h"

using namespace ecmech;

//All of the parallelization operations are within the getResponse function of the material class.
void mat_model_kernel(const ecmech::matModelBase* mat_model_base,
                      const int nqpts, const double dt,
                      const int nstatev, double* state_vars_array,
                      double* stress_svec_p_array, double* d_svec_p_array,
                      double* w_vec_array, double* ddsdde_array,
                      double* vol_ratio_array, double* eng_int_array){
   //This really isn't an efficient way to do this.
   double* temp_array = memoryManager::allocate<double>(nqpts);// 300.;
   double* sdd_array = memoryManager::allocate<double>(nqpts * 2);

   for (int i = 0; i < nqpts; i++) {
      temp_array[i] = 300.;
   }

   mat_model_base->getResponse(dt, d_svec_p_array, w_vec_array, vol_ratio_array,
                               eng_int_array, stress_svec_p_array, state_vars_array,
                               temp_array, sdd_array, ddsdde_array, nqpts);
   //We'll later just keep these around for good but for now let's just make sure things work.
   memoryManager::deallocate(temp_array);
   memoryManager::deallocate(sdd_array);
}//end of mat_model_kernel
