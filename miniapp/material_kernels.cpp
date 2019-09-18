#include "material_kernels.h"

using namespace ecmech;

//We want to keep all of the functions in this namespace private and not viewable to outside files.
namespace {
   class matModelBase_test
   {
      public:
         __host__ __device__
         matModelBase_test() {}

         __host__ __device__
         virtual void getResponse(const double  & dt,
                                  const double  * defRateV,
                                  const double  * spinV,
                                  const double  * volRatioV,
                                  double  * eIntV,
                                  double  * stressSvecPV,
                                  double  * histV,
                                  double  * tkelvV,
                                  double  * sddV,
                                  double  * mtanSDV,
                                  const int    & nPassed) const = 0;

         __host__ __device__
         ~matModelBase_test() {}
   };

   class matModel_test : public matModelBase_test
   {
      private:
         double var;
      public:
         __host__ __device__
         matModel_test() {}

         __host__ __device__
         virtual void getResponse(const double  & dt,
                                  const double  * defRateV,
                                  const double  * spinV,
                                  const double  * volRatioV,
                                  double  * eIntV,
                                  double  * stressSvecPV,
                                  double  * histV,
                                  double  * tkelvV,
                                  double  * sddV,
                                  double  * mtanSDV,
                                  const int    & nPassed) const override final {
            eIntV[0] = 1.0;
            stressSvecPV[0] = var;
         }

         __host__ void init_data(double num){
            var = num;
         }

         __host__ __device__
         ~matModel_test() {}
   };

   void mat_model_kernel_cpu(const ecmech::matModelBase* mat_model_base,
                             const int nqpts, const double dt,
                             const int nstatev, double* state_vars_array,
                             double* stress_svec_p_array, double* d_svec_p_array,
                             double* w_vec_array, double* ddsdde_array,
                             double* vol_ratio_array, double* eng_int_array){
      //The current thought is to potentially have each thread do like a group of items
      //and then go onto the next item. However, it might be seen that actually it's more efficient to have
      //each thread do one item and then go onto the next.
      const int num_items = 1;
      const int num_passes = ceil(nqpts / num_items);
      const int remainder = nqpts % num_items;
      double temp = 300.;
      //This was necessary since RAJA kept capturing it as a constant value
      //which is not what we wanted...
      double* temp_ref = &(temp);

      RAJA::RangeSegment default_range(0, num_passes);

      RAJA::forall<RAJA::loop_exec>(default_range, [ = ](int ipass) {
         // for(int ipass = 0; ipass < num_passes; ipass++){
         int i_qpts = ipass * num_items;
         //Might want to eventually set these all up using RAJA views. It might simplify
         //things later on.
         //These are our inputs
         double* state_vars = &(state_vars_array[i_qpts * nstatev]);
         //Here is all of our ouputs
         double* ddsdde = &(ddsdde_array[i_qpts * ecmech::nsvec * ecmech::nsvec]);
         double* eng_int = &(eng_int_array[i_qpts * ecmech::ne]);
         double* w_vec = &(w_vec_array[i_qpts * ecmech::nwvec]);
         double* vol_ratio = &(vol_ratio_array[i_qpts * ecmech::nvr]);
         //A few variables are set up as the 6-vec deviatoric + tr(tens) values
         int ind_svecp = i_qpts * ecmech::nsvp;
         double* stress_svec_p = &(stress_svec_p_array[ind_svecp]);
         double* d_svec_p = &(d_svec_p_array[ind_svecp]);
         double sdd[2];

         if (remainder == 0 | ipass < (num_passes - 1)) {
            mat_model_base->getResponse(dt, d_svec_p, w_vec, vol_ratio,
                                        eng_int, stress_svec_p, state_vars,
                                        temp_ref, &sdd[0], ddsdde, num_items);
         }
         else {
            mat_model_base->getResponse(dt, d_svec_p, w_vec, vol_ratio,
                                        eng_int, stress_svec_p, state_vars,
                                        temp_ref, &sdd[0], ddsdde, remainder);
         }
      });//end of npass loop
   }

#if defined(RAJA_ENABLE_OPENMP)
   void mat_model_kernel_openmp(const ecmech::matModelBase* mat_model_base,
                                const int nqpts, const double dt,
                                const int nstatev, double* state_vars_array,
                                double* stress_svec_p_array, double* d_svec_p_array,
                                double* w_vec_array, double* ddsdde_array,
                                double* vol_ratio_array, double* eng_int_array){
      //The current thought is to potentially have each thread do like a group of items
      //and then go onto the next item. However, it might be seen that actually it's more efficient to have
      //each thread do one item and then go onto the next.
      const int num_items = 1;
      const int num_passes = ceil(nqpts / num_items);
      const int remainder = nqpts % num_items;
      double temp = 300.;
      //This was necessary since RAJA kept capturing it as a constant value
      //which is not what we wanted...
      double* temp_ref = &(temp);

      RAJA::RangeSegment default_range(0, num_passes);

      RAJA::forall<RAJA::omp_parallel_for_exec>(default_range, [ = ](int ipass) {
         // for(int ipass = 0; ipass < num_passes; ipass++){
         int i_qpts = ipass * num_items;
         //Might want to eventually set these all up using RAJA views. It might simplify
         //things later on.
         //These are our inputs
         double* state_vars = &(state_vars_array[i_qpts * nstatev]);
         //Here is all of our ouputs
         double* ddsdde = &(ddsdde_array[i_qpts * ecmech::nsvec * ecmech::nsvec]);
         double* eng_int = &(eng_int_array[i_qpts * ecmech::ne]);
         double* w_vec = &(w_vec_array[i_qpts * ecmech::nwvec]);
         double* vol_ratio = &(vol_ratio_array[i_qpts * ecmech::nvr]);
         //A few variables are set up as the 6-vec deviatoric + tr(tens) values
         int ind_svecp = i_qpts * ecmech::nsvp;
         double* stress_svec_p = &(stress_svec_p_array[ind_svecp]);
         double* d_svec_p = &(d_svec_p_array[ind_svecp]);
         double sdd[2];

         if (remainder == 0 | ipass < (num_passes - 1)) {
            mat_model_base->getResponse(dt, d_svec_p, w_vec, vol_ratio,
                                        eng_int, stress_svec_p, state_vars,
                                        temp_ref, &sdd[0], ddsdde, num_items);
         }
         else {
            mat_model_base->getResponse(dt, d_svec_p, w_vec, vol_ratio,
                                        eng_int, stress_svec_p, state_vars,
                                        temp_ref, &sdd[0], ddsdde, remainder);
         }
      });//end of npass loop
   }//end of mat_model_openmp

#endif

#if defined(RAJA_ENABLE_CUDA)
   void mat_model_kernel_cuda(const ecmech::matModelBase* mat_model_base,
                              const int nqpts, const double dt,
                              const int nstatev, double* state_vars_array,
                              double* stress_svec_p_array, double* d_svec_p_array,
                              double* w_vec_array, double* ddsdde_array,
                              double* vol_ratio_array, double* eng_int_array){
      //The current thought is to potentially have each thread do like a group of items
      //and then go onto the next item. However, it might be seen that actually it's more efficient to have
      //each thread do one item and then go onto the next.
      const int num_items = 1;
      const int num_passes = ceil(nqpts / num_items);
      const int remainder = nqpts % num_items;
      double temp = 300.;
      //This was necessary since RAJA kept capturing it as a constant value
      //which is not what we wanted...
      double* temp_ref = &(temp);

      RAJA::RangeSegment default_range(0, num_passes);

      //The block size could change during compilation if we decide that there's better ways of doing this
      //or if we want to test different sizes
      //static const int block_size = 128;

      const int nsvec = ecmech::nsvec;
      const int ne = ecmech::ne;
      const int nwvec = ecmech::nwvec;
      const int nvr = ecmech::nvr;
      const int nsvp = ecmech::nsvp;
      /*
      matModel_test* arg =  memoryManager::allocate<matModel_test>(1);

      RAJA::forall<RAJA::cuda_exec<1>>(RAJA::RangeSegment(0, 1), [ = ] RAJA_HOST_DEVICE (int) {
     new(arg) matModel_test();
      });

      arg->init_data(1.0);

      matModelBase_test* arg2 = dynamic_cast<matModelBase_test*>(arg);
      */
      RAJA::forall<RAJA::cuda_exec<128> >(default_range, [ = ] RAJA_DEVICE(int ipass) {
         // for(int ipass = 0; ipass < num_passes; ipass++){
         int i_qpts = ipass * num_items;
         //double temp = 300;
         //Might want to eventually set these all up using RAJA views. It might simplify
         //things later on.
         //These are our inputs
         double* state_vars = &(state_vars_array[i_qpts * nstatev]);
         //Here is all of our ouputs
         double* ddsdde = &(ddsdde_array[i_qpts * nsvec * nsvec]);
         double* eng_int = &(eng_int_array[i_qpts * ne]);
         double* w_vec = &(w_vec_array[i_qpts * nwvec]);
         double* vol_ratio = &(vol_ratio_array[i_qpts * nvr]);
         //A few variables are set up as the 6-vec deviatoric + tr(tens) values
         int ind_svecp = i_qpts * nsvp;
         double* stress_svec_p = &(stress_svec_p_array[ind_svecp]);
         double* d_svec_p = &(d_svec_p_array[ind_svecp]);
         double sdd[2];

         //sdd[0] = 2;
         /*
         arg2->getResponse(dt, d_svec_p, w_vec, vol_ratio,
               eng_int, stress_svec_p, state_vars,
               &temp, &sdd[0], ddsdde, num_items);
         */
         if (remainder == 0 | ipass < (num_passes - 1)) {
            mat_model_base->getResponse(dt, d_svec_p, w_vec, vol_ratio,
                                        eng_int, stress_svec_p, state_vars,
                                        temp_ref, &sdd[0], ddsdde, num_items);
         }
         else {
            mat_model_base->getResponse(dt, d_svec_p, w_vec, vol_ratio,
                                        eng_int, stress_svec_p, state_vars,
                                        temp_ref, &sdd[0], ddsdde, remainder);
         }
      });//end of npass loop
         /*
         RAJA::forall<RAJA::cuda_exec<1>>(RAJA::RangeSegment(0, 1), [ = ] RAJA_HOST_DEVICE (int) {
       arg->~matModel_test();
         });

         memoryManager::deallocate(arg); */
   }//end of mat_model_cuda

#endif
}

void mat_model_kernel(Accelerator accel, const ecmech::matModelBase* mat_model_base,
                      const int nqpts, const double dt,
                      const int nstatev, double* state_vars_array,
                      double* stress_svec_p_array, double* d_svec_p_array,
                      double* w_vec_array, double* ddsdde_array,
                      double* vol_ratio_array, double* eng_int_array){
   #if defined(RAJA_ENABLE_OPENMP)
   if (accel == Accelerator::OPENMP) {
      mat_model_kernel_openmp(mat_model_base, nqpts, dt, nstatev, state_vars_array,
                              stress_svec_p_array, d_svec_p_array, w_vec_array, ddsdde_array,
                              vol_ratio_array, eng_int_array);
   }
   #endif
   #if defined(RAJA_ENABLE_CUDA)
   if (accel == Accelerator::CUDA) {
      mat_model_kernel_cuda(mat_model_base, nqpts, dt, nstatev, state_vars_array,
                            stress_svec_p_array, d_svec_p_array, w_vec_array, ddsdde_array,
                            vol_ratio_array, eng_int_array);
   }
   #endif
   if (accel == Accelerator::CPU) {
      mat_model_kernel_cpu(mat_model_base, nqpts, dt, nstatev, state_vars_array,
                           stress_svec_p_array, d_svec_p_array, w_vec_array, ddsdde_array,
                           vol_ratio_array, eng_int_array);
   }
}//end of mat_model_kernel
