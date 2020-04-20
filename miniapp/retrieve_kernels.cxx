#include "retrieve_kernels.h"


using namespace ecmech;

namespace {
   void retrieve_data_cpu(const int nqpts, const int nstatev,
                          const double* stress_svec_p_array, const double* vol_ratio_array,
                          const double* eng_int_array, double* state_vars_array,
                          double* stress_array){
      const int ind_int_eng = nstatev - ecmech::ne;
      const int ind_vols = ind_int_eng - 1;

      RAJA::RangeSegment default_range(0, nqpts);

      RAJA::forall<RAJA::loop_exec>(default_range, [ = ](int i_qpts) {
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
         std::copy(stress_svec_p, stress_svec_p + ecmech::nsvec, stress);
         stress[0] += stress_mean;
         stress[1] += stress_mean;
         stress[2] += stress_mean;
      }); // end of qpts loop
   } // end of retrieve_data_cpu

#if defined(RAJA_ENABLE_OPENMP)
   void retrieve_data_openmp(const int nqpts, const int nstatev,
                             const double* stress_svec_p_array, const double* vol_ratio_array,
                             const double* eng_int_array, double* state_vars_array,
                             double* stress_array){
      const int ind_int_eng = nstatev - ecmech::ne;
      const int ind_vols = ind_int_eng - 1;

      RAJA::RangeSegment default_range(0, nqpts);

      RAJA::forall<RAJA::omp_parallel_for_exec>(default_range, [ = ](int i_qpts) {
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
         std::copy(stress_svec_p, stress_svec_p + ecmech::nsvec, stress);
         stress[0] += stress_mean;
         stress[1] += stress_mean;
         stress[2] += stress_mean;
      }); // end of qpts loop
   } // end of retrieve_data_openmp

#endif

#if defined(RAJA_ENABLE_CUDA)
   void retrieve_data_cuda(const int nqpts, const int nstatev,
                           const double* stress_svec_p_array, const double* vol_ratio_array,
                           const double* eng_int_array, double* state_vars_array,
                           double* stress_array){
      const int ind_int_eng = nstatev - ecmech::ne;
      const int ind_vols = ind_int_eng - 1;

      RAJA::RangeSegment default_range(0, nqpts);

      RAJA::forall<RAJA::cuda_exec<384> >(default_range, [ = ] RAJA_DEVICE(int i_qpts) {
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
      }); // end of qpts loop
   } // end of retrieve_data_cuda

#endif
} // end of private namespace

// This will then be the final function/kernel to save off all the data at
// each time step.
void retrieve_data(ecmech::VectorizationStrategy accel, const int nqpts, const int nstatev,
                   const double* stress_svec_p_array, const double* vol_ratio_array,
                   const double* eng_int_array, double* state_vars_array,
                   double* stress_array){
   switch ( accel ) {
#if defined(RAJA_ENABLE_OPENMP)
   case ecmech::VectorizationStrategy::OPENMP :
   {
      retrieve_data_openmp(nqpts, nstatev, stress_svec_p_array, vol_ratio_array,
                           eng_int_array, state_vars_array, stress_array);
   }
   break;
#endif
#if defined(RAJA_ENABLE_CUDA)
   case ecmech::VectorizationStrategy::CUDA :
   {
      retrieve_data_cuda(nqpts, nstatev, stress_svec_p_array, vol_ratio_array,
                         eng_int_array, state_vars_array, stress_array);
   }
   break;
#endif
   case ecmech::VectorizationStrategy::CPU :
   default :
   {
      retrieve_data_cpu(nqpts, nstatev, stress_svec_p_array, vol_ratio_array,
                        eng_int_array, state_vars_array, stress_array);
   }
   break;
   }
   
} // end of retrieve_data

