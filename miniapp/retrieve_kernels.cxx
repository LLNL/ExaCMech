#include "retrieve_kernels.h"

#include "ECMech_evptnWrap.h"
#include <math.h>

// This will then be the final function/kernel to save off all the data at
// each time step.
void retrieve_data(const int nqpts, const int nstatev,
                   const double* stress_svec_p_array, const double* vol_ratio_array,
                   const double* eng_int_array, double* state_vars_array,
                   double* stress_array) {
   const int ind_int_eng = nstatev - ecmech::ne;
   const int ind_vols = ind_int_eng - 1;

   snls::forall(0, nqpts, [=]
      __ecmech_hdev__
      (int i_qpts)
   {
      // These are our outputs
      double* state_vars = &(state_vars_array[i_qpts * nstatev]);
      double* stress = &(stress_array[i_qpts * ecmech::nsvec]);
      // Here is all of our ouputs
      const double* eng_int = &(eng_int_array[i_qpts * ecmech::ne]);
      const double* vol_ratio = &(vol_ratio_array[i_qpts * ecmech::nvr]);
      // A few variables are set up as the 6-vec deviatoric + tr(tens) values
      int ind_svecp = i_qpts * ecmech::nsvp;
      const double* stress_svec_p = &(stress_svec_p_array[ind_svecp]);

      // We need to update our state variables to include the volume ratio and
      // internal energy portions
      state_vars[ind_vols] = vol_ratio[1];
      for (int i = 0; i < ecmech::ne; i++) {
         state_vars[ind_int_eng + i] = eng_int[i];
      }

      // Here we're converting back from our deviatoric + pressure representation of our
      // Cauchy stress back to the Voigt notation of stress.
      double stress_mean = -stress_svec_p[ecmech::iSvecP];
      for (int i = 0; i < ecmech::nsvec; i++) {
         stress[i] = stress_svec_p[i];
      }
      stress[0] += stress_mean;
      stress[1] += stress_mean;
      stress[2] += stress_mean;
   });
} // end of retrieve_data

