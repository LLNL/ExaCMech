/**
 * @file retrieve_kernels.cxx
 * @brief Implementation of `retrieve_data` (declared in `retrieve_kernels.h`).
 */

#include "retrieve_kernels.h"

#include "ECMech_evptnWrap.h"
#include <math.h>

// See retrieve_kernels.h for the full @brief/@param doc.
void retrieve_data(const int nqpts, const int nstatev,
                   const double* cauchy_stress_d6p_array, const double* rel_vol_ratios_array,
                   const double* internal_energy_array, double* state_vars_array,
                   double* cauchy_stress_array) {
   // The miniapp's state_vars layout tacks the volume-ratio and internal-energy slots
   // onto the end, right after the model's own history block (see init_data in
   // setup_kernels.cxx) -- so they sit at the last (ind_int_eng..ind_int_eng+ne-1) and
   // (ind_vols) offsets of the nstatev-wide per-point record.
   const int ind_int_eng = nstatev - ecmech::ne;
   const int ind_vols = ind_int_eng - 1;

   // snls::forall is SNLS's CPU/OpenMP/GPU-portable parallel-for; the execution back end
   // was already selected when the model was configured (matModelBase::setExecutionStrategy),
   // and getResponseECM ran on the same points with the same back end just before this.
   snls::forall(0, nqpts, [=]
      __ecmech_hdev__
      (int i_qpts)
   {
      // These are our outputs
      double* state_vars = &(state_vars_array[i_qpts * nstatev]);
      double* cauchy_stress = &(cauchy_stress_array[i_qpts * ecmech::nsvec]);
      // Here is all of our ouputs
      const double* internal_energy = &(internal_energy_array[i_qpts * ecmech::ne]);
      const double* rel_vol_ratios = &(rel_vol_ratios_array[i_qpts * ecmech::nvr]);
      // A few variables are set up as the 6-vec deviatoric + tr(tens) values
      int ind_svecp = i_qpts * ecmech::nsvp;
      const double* cauchy_stress_d6p = &(cauchy_stress_d6p_array[ind_svecp]);

      // We need to update our state variables to include the volume ratio and
      // internal energy portions
      state_vars[ind_vols] = rel_vol_ratios[1];
      for (int i = 0; i < ecmech::ne; i++) {
         state_vars[ind_int_eng + i] = internal_energy[i];
      }

      // Here we're converting back from our deviatoric + pressure representation of our
      // Cauchy stress back to the Voigt notation of stress. cauchy_stress_d6p[iSvecP]
      // holds -mean_stress (see ECMech_util.h's svecp/pressure convention), so negating
      // it and adding it onto the three normal components reconstructs the full
      // symmetric Voigt tensor.
      double stress_mean = -cauchy_stress_d6p[ecmech::iSvecP];
      for (int i = 0; i < ecmech::nsvec; i++) {
         cauchy_stress[i] = cauchy_stress_d6p[i];
      }
      cauchy_stress[0] += stress_mean;
      cauchy_stress[1] += stress_mean;
      cauchy_stress[2] += stress_mean;
   });
} // end of retrieve_data

