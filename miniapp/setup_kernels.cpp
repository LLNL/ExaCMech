#include "setup_kernels.h"

using namespace ecmech;

//We want to keep all of the functions in this namespace private and not viewable to outside files.
namespace {
   void init_data_cpu(const double* ori, const std::vector<double>& histInit_vec,
                      const int nqpts, const int num_hardness, const int ind_gdot,
                      const int num_slip, const int vdim, double* state_vars){
      //Probably going to replace this all with just a constant index array and
      //a constant num_item array. We can then use that to pass into other function handles
      //easily.
      const int ind_dp_eff = ecmech::evptn::iHistA_shrateEff;
      const int ind_eql_pl_strain = ecmech::evptn::iHistA_shrEff;
      const int ind_num_evals = ecmech::evptn::iHistA_nFEval;
      const int ind_dev_elas_strain = ecmech::evptn::iHistLbE;
      const int ind_quats = ecmech::evptn::iHistLbQ;
      const int ind_hardness = ecmech::evptn::iHistLbH;
      //The number of vols -> we actually only need to save the previous time step value
      //instead of all 4 values used in the evalModel. The rest can be calculated from
      //this value.
      const int num_vols = 1;
      const int ind_vols = ind_gdot + num_slip;
      const int ind_int_eng = ind_vols + num_vols;


      //We're going to use RAJA here to initialize everything all at once
      //We should be able to just use OpenMP here. Since, we can assume everything is on the
      //host originally. We'll later want to migrate everything over to the GPU if we're
      //running things on there.
      RAJA::RangeSegment default_range(0, nqpts);

      RAJA::forall<RAJA::loop_exec>(default_range, [ = ](int i) {
         // for (int i = 0; i < nqpts; i++) {
         int ind = i * vdim;
         int ind_ori = i * ecmech::qdim;

         state_vars[ind + ind_dp_eff] = histInit_vec[ind_dp_eff];
         state_vars[ind + ind_eql_pl_strain] = histInit_vec[ind_eql_pl_strain];
         state_vars[ind + ind_num_evals] = histInit_vec[ind_num_evals];
         //Normally, we might have this assigned as 0 but that may not always
         //be the case.
         state_vars[ind + ind_vols] = 1.0;
         //Due to the smallness of these loops, I would expect several of these to be
         //completely unrolled by the compiler.
         for (int j = 0; j < ecmech::qdim; j++) {
            state_vars[ind + ind_quats + j] = ori[ind_ori + j];
         }

         for (int j = 0; j < ecmech::ntvec; j++) {
            state_vars[ind + ind_dev_elas_strain + j] = histInit_vec[ind_dev_elas_strain + j];
         }

         for (int j = 0; j < num_slip; j++) {
            state_vars[ind + ind_gdot + j] = histInit_vec[ind_gdot + j];
         }

         for (int j = 0; j < num_hardness; j++) {
            state_vars[ind + ind_hardness] = histInit_vec[ind_hardness + j];
         }

         //Normally, we might have this assigned as 0 but that may not always
         //be the case.
         for (int j = 0; j < ecmech::ne; j++) {
            state_vars[ind + ind_int_eng + j] = 0.0;
         }
      });  //end of qpt loop
   } //end of init_data_cpu

   void setup_data_cpu(const int nqpts, const int nstatev,
                       const double dt, const double* vel_grad_array,
                       const double* stress_array, const double* state_vars_array,
                       double* stress_svec_p_array, double* d_svec_p_array,
                       double* w_vec_array, double* ddsdde_array,
                       double* vol_ratio_array, double* eng_int_array){
      //vgrad is kinda a pain to deal with as a raw 1d array, so we're
      //going to just use a RAJA view here. The data is taken to be in col. major format.
      //It might be nice to eventually create a type alias for the below or
      //maybe something like it.

      const int ind_int_eng = nstatev - ecmech::ne;
      const int ind_vols = ind_int_eng - 1;

      const int DIM = 3;
      std::array<RAJA::idx_t, DIM> perm {{ 2, 1, 0 } };
      RAJA::Layout<DIM> layout = RAJA::make_permuted_layout({{ ecmech::ndim, ecmech::ndim, nqpts } }, perm);
      RAJA::View<const double, RAJA::Layout<DIM, RAJA::Index_type, 0> > vgrad_view(vel_grad_array, layout);

      //All of the below we could setup in one big RAJA loop/kernel
      // RAJA::RangeSegment default_range(0, nqpts);
      RAJA::RangeSegment default_range(0, nqpts);

      RAJA::forall<RAJA::loop_exec>(default_range, [ = ](int i_qpts) {
         // for (int i_qpts = 0; i_qpts < nqpts; i_qpts++) {
         //Might want to eventually set these all up using RAJA views. It might simplify
         //things later on.
         //These are our inputs
         const double* state_vars = &(state_vars_array[i_qpts * nstatev]);
         const double* stress = &(stress_array[i_qpts * ecmech::nsvec]);
         //Here is all of our ouputs
         double* ddsdde = &(ddsdde_array[i_qpts * ecmech::nsvec * ecmech::nsvec]);
         double* eng_int = &(eng_int_array[i_qpts * ecmech::ne]);
         double* w_vec = &(w_vec_array[i_qpts * ecmech::nwvec]);
         double* vol_ratio = &(vol_ratio_array[i_qpts * ecmech::nvr]);
         //A few variables are set up as the 6-vec deviatoric + tr(tens) values
         int ind_svecp = i_qpts * ecmech::nsvp;
         double* stress_svec_p = &(stress_svec_p_array[ind_svecp]);
         double* d_svec_p = &(d_svec_p_array[ind_svecp]);

         // initialize 6x6 2d arrays all to 0
         for (int i = 0; i < ecmech::nsvec; i++) {
            for (int j = 0; j < ecmech::nsvec; j++) {
               ddsdde[(i * ecmech::nsvec) +j] = 0.0;
            }
         }

         for (int i = 0; i < ecmech::ne; i++) {
            eng_int[i] = state_vars[ind_int_eng + i];
         }

         //Here we have the skew portion of our velocity gradient as represented as an
         //axial vector.
         w_vec[0] = 0.5 * (vgrad_view(2, 1, i_qpts) - vgrad_view(1, 2, i_qpts));
         w_vec[1] = 0.5 * (vgrad_view(0, 2, i_qpts) - vgrad_view(2, 0, i_qpts));
         w_vec[2] = 0.5 * (vgrad_view(1, 0, i_qpts) - vgrad_view(0, 1, i_qpts));

         //Really we're looking at the negative of J but this will do...
         double d_mean = -ecmech::onethird * (vgrad_view(0, 0, i_qpts) + vgrad_view(1, 1, i_qpts) + vgrad_view(2, 2, i_qpts));
         //The 1st 6 components are the symmetric deviatoric portion of our velocitu gradient
         //The last value is simply the trace of the deformation rate
         d_svec_p[0] = vgrad_view(0, 0, i_qpts) + d_mean;
         d_svec_p[1] = vgrad_view(1, 1, i_qpts) + d_mean;
         d_svec_p[2] = vgrad_view(2, 2, i_qpts) + d_mean;
         d_svec_p[3] = 0.5 * (vgrad_view(2, 1, i_qpts) + vgrad_view(1, 2, i_qpts));
         d_svec_p[4] = 0.5 * (vgrad_view(2, 0, i_qpts) + vgrad_view(0, 2, i_qpts));
         d_svec_p[5] = 0.5 * (vgrad_view(1, 0, i_qpts) + vgrad_view(0, 1, i_qpts));
         d_svec_p[6] = -3 * d_mean;

         vol_ratio[0] = state_vars[ind_vols];
         vol_ratio[1] = vol_ratio[0] * exp(d_svec_p[ecmech::iSvecP] * dt);
         vol_ratio[3] = vol_ratio[1] - vol_ratio[0];
         vol_ratio[2] = vol_ratio[3] / (dt * 0.5 * (vol_ratio[0] + vol_ratio[1]));

         std::copy(stress, stress + ecmech::nsvec, stress_svec_p);

         double stress_mean = -ecmech::onethird * (stress[0] + stress[1] + stress[2]);
         stress_svec_p[0] += stress_mean;
         stress_svec_p[1] += stress_mean;
         stress_svec_p[2] += stress_mean;
         stress_svec_p[ecmech::iSvecP] = stress_mean;
      });//end of qpt loop
   }//end setup_data_cpu

#if defined(RAJA_ENABLE_OPENMP)
   void init_data_openmp(const double* ori, const std::vector<double>& histInit_vec,
                         const int nqpts, const int num_hardness, const int ind_gdot,
                         const int num_slip, const int vdim, double* state_vars){
      //Probably going to replace this all with just a constant index array and
      //a constant num_item array. We can then use that to pass into other function handles
      //easily.
      const int ind_dp_eff = ecmech::evptn::iHistA_shrateEff;
      const int ind_eql_pl_strain = ecmech::evptn::iHistA_shrEff;
      const int ind_num_evals = ecmech::evptn::iHistA_nFEval;
      const int ind_dev_elas_strain = ecmech::evptn::iHistLbE;
      const int ind_quats = ecmech::evptn::iHistLbQ;
      const int ind_hardness = ecmech::evptn::iHistLbH;
      //The number of vols -> we actually only need to save the previous time step value
      //instead of all 4 values used in the evalModel. The rest can be calculated from
      //this value.
      const int num_vols = 1;
      const int ind_vols = ind_gdot + num_slip;
      const int ind_int_eng = ind_vols + num_vols;


      //We're going to use RAJA here to initialize everything all at once
      //We should be able to just use OpenMP here. Since, we can assume everything is on the
      //host originally. We'll later want to migrate everything over to the GPU if we're
      //running things on there.
      RAJA::RangeSegment default_range(0, nqpts);

      RAJA::forall<RAJA::omp_parallel_for_exec>(default_range, [ = ](int i) {
         // for (int i = 0; i < nqpts; i++) {
         int ind = i * vdim;
         int ind_ori = i * ecmech::qdim;

         state_vars[ind + ind_dp_eff] = histInit_vec[ind_dp_eff];
         state_vars[ind + ind_eql_pl_strain] = histInit_vec[ind_eql_pl_strain];
         state_vars[ind + ind_num_evals] = histInit_vec[ind_num_evals];
         //Normally, we might have this assigned as 0 but that may not always
         //be the case.
         state_vars[ind + ind_vols] = 1.0;
         //Due to the smallness of these loops, I would expect several of these to be
         //completely unrolled by the compiler.
         for (int j = 0; j < ecmech::qdim; j++) {
            state_vars[ind + ind_quats + j] = ori[ind_ori + j];
         }

         for (int j = 0; j < ecmech::ntvec; j++) {
            state_vars[ind + ind_dev_elas_strain + j] = histInit_vec[ind_dev_elas_strain + j];
         }

         for (int j = 0; j < num_slip; j++) {
            state_vars[ind + ind_gdot + j] = histInit_vec[ind_gdot + j];
         }

         for (int j = 0; j < num_hardness; j++) {
            state_vars[ind + ind_hardness] = histInit_vec[ind_hardness + j];
         }

         //Normally, we might have this assigned as 0 but that may not always
         //be the case.
         for (int j = 0; j < ecmech::ne; j++) {
            state_vars[ind + ind_int_eng + j] = 0.0;
         }
      });  //end of qpt loop
   } //end of init_data_openmp

   void setup_data_openmp(const int nqpts, const int nstatev,
                          const double dt, const double* vel_grad_array,
                          const double* stress_array, const double* state_vars_array,
                          double* stress_svec_p_array, double* d_svec_p_array,
                          double* w_vec_array, double* ddsdde_array,
                          double* vol_ratio_array, double* eng_int_array){
      //vgrad is kinda a pain to deal with as a raw 1d array, so we're
      //going to just use a RAJA view here. The data is taken to be in col. major format.
      //It might be nice to eventually create a type alias for the below or
      //maybe something like it.

      const int ind_int_eng = nstatev - ecmech::ne;
      const int ind_vols = ind_int_eng - 1;

      const int DIM = 3;
      std::array<RAJA::idx_t, DIM> perm {{ 2, 1, 0 } };
      RAJA::Layout<DIM> layout = RAJA::make_permuted_layout({{ ecmech::ndim, ecmech::ndim, nqpts } }, perm);
      RAJA::View<const double, RAJA::Layout<DIM, RAJA::Index_type, 0> > vgrad_view(vel_grad_array, layout);

      //All of the below we could setup in one big RAJA loop/kernel
      // RAJA::RangeSegment default_range(0, nqpts);
      RAJA::RangeSegment default_range(0, nqpts);

      RAJA::forall<RAJA::omp_parallel_for_exec>(default_range, [ = ](int i_qpts) {
         // for (int i_qpts = 0; i_qpts < nqpts; i_qpts++) {
         //Might want to eventually set these all up using RAJA views. It might simplify
         //things later on.
         //These are our inputs
         const double* state_vars = &(state_vars_array[i_qpts * nstatev]);
         const double* stress = &(stress_array[i_qpts * ecmech::nsvec]);
         //Here is all of our ouputs
         double* ddsdde = &(ddsdde_array[i_qpts * ecmech::nsvec * ecmech::nsvec]);
         double* eng_int = &(eng_int_array[i_qpts * ecmech::ne]);
         double* w_vec = &(w_vec_array[i_qpts * ecmech::nwvec]);
         double* vol_ratio = &(vol_ratio_array[i_qpts * ecmech::nvr]);
         //A few variables are set up as the 6-vec deviatoric + tr(tens) values
         int ind_svecp = i_qpts * ecmech::nsvp;
         double* stress_svec_p = &(stress_svec_p_array[ind_svecp]);
         double* d_svec_p = &(d_svec_p_array[ind_svecp]);

         // initialize 6x6 2d arrays all to 0
         for (int i = 0; i < ecmech::nsvec; i++) {
            for (int j = 0; j < ecmech::nsvec; j++) {
               ddsdde[(i * ecmech::nsvec) +j] = 0.0;
            }
         }

         for (int i = 0; i < ecmech::ne; i++) {
            eng_int[i] = state_vars[ind_int_eng + i];
         }

         //Here we have the skew portion of our velocity gradient as represented as an
         //axial vector.
         w_vec[0] = 0.5 * (vgrad_view(2, 1, i_qpts) - vgrad_view(1, 2, i_qpts));
         w_vec[1] = 0.5 * (vgrad_view(0, 2, i_qpts) - vgrad_view(2, 0, i_qpts));
         w_vec[2] = 0.5 * (vgrad_view(1, 0, i_qpts) - vgrad_view(0, 1, i_qpts));

         //Really we're looking at the negative of J but this will do...
         double d_mean = -ecmech::onethird * (vgrad_view(0, 0, i_qpts) + vgrad_view(1, 1, i_qpts) + vgrad_view(2, 2, i_qpts));
         //The 1st 6 components are the symmetric deviatoric portion of our velocitu gradient
         //The last value is simply the trace of the deformation rate
         d_svec_p[0] = vgrad_view(0, 0, i_qpts) + d_mean;
         d_svec_p[1] = vgrad_view(1, 1, i_qpts) + d_mean;
         d_svec_p[2] = vgrad_view(2, 2, i_qpts) + d_mean;
         d_svec_p[3] = 0.5 * (vgrad_view(2, 1, i_qpts) + vgrad_view(1, 2, i_qpts));
         d_svec_p[4] = 0.5 * (vgrad_view(2, 0, i_qpts) + vgrad_view(0, 2, i_qpts));
         d_svec_p[5] = 0.5 * (vgrad_view(1, 0, i_qpts) + vgrad_view(0, 1, i_qpts));
         d_svec_p[6] = -3 * d_mean;

         vol_ratio[0] = state_vars[ind_vols];
         vol_ratio[1] = vol_ratio[0] * exp(d_svec_p[ecmech::iSvecP] * dt);
         vol_ratio[3] = vol_ratio[1] - vol_ratio[0];
         vol_ratio[2] = vol_ratio[3] / (dt * 0.5 * (vol_ratio[0] + vol_ratio[1]));

         std::copy(stress, stress + ecmech::nsvec, stress_svec_p);

         double stress_mean = -ecmech::onethird * (stress[0] + stress[1] + stress[2]);
         stress_svec_p[0] += stress_mean;
         stress_svec_p[1] += stress_mean;
         stress_svec_p[2] += stress_mean;
         stress_svec_p[ecmech::iSvecP] = stress_mean;
      });//end of qpt loop
   }//end setup_data_openmp

#endif

#if defined(RAJA_ENABLE_CUDA)

   void setup_data_cuda(const int nqpts, const int nstatev,
                        const double dt, const double* vel_grad_array,
                        const double* stress_array, const double* state_vars_array,
                        double* stress_svec_p_array, double* d_svec_p_array,
                        double* w_vec_array, double* ddsdde_array,
                        double* vol_ratio_array, double* eng_int_array){
      //vgrad is kinda a pain to deal with as a raw 1d array, so we're
      //going to just use a RAJA view here. The data is taken to be in col. major format.
      //It might be nice to eventually create a type alias for the below or
      //maybe something like it.

      const int ind_int_eng = nstatev - ecmech::ne;
      const int ind_vols = ind_int_eng - 1;

      const int DIM = 3;
      std::array<RAJA::idx_t, DIM> perm {{ 2, 1, 0 } };
      RAJA::Layout<DIM> layout = RAJA::make_permuted_layout({{ ecmech::ndim, ecmech::ndim, nqpts } }, perm);
      RAJA::View<const double, RAJA::Layout<DIM, RAJA::Index_type, 0> > vgrad_view(vel_grad_array, layout);

      //All of the below we could setup in one big RAJA loop/kernel
      // RAJA::RangeSegment default_range(0, nqpts);
      RAJA::RangeSegment default_range(0, nqpts);
      //The block size could change during compilation if we decide that there's better ways of doing this
      //or if we want to test different sizes
      const int block_size = 256;

      RAJA::forall<RAJA::cuda_exec<256> >(default_range, [ = ] RAJA_DEVICE(int i_qpts) {
         // for (int i_qpts = 0; i_qpts < nqpts; i_qpts++) {
         //Might want to eventually set these all up using RAJA views. It might simplify
         //things later on.
         //These are our inputs
         const double* state_vars = &(state_vars_array[i_qpts * nstatev]);
         const double* stress = &(stress_array[i_qpts * ecmech::nsvec]);
         //Here is all of our ouputs
         double* ddsdde = &(ddsdde_array[i_qpts * ecmech::nsvec * ecmech::nsvec]);
         double* eng_int = &(eng_int_array[i_qpts * ecmech::ne]);
         double* w_vec = &(w_vec_array[i_qpts * ecmech::nwvec]);
         double* vol_ratio = &(vol_ratio_array[i_qpts * ecmech::nvr]);
         //A few variables are set up as the 6-vec deviatoric + tr(tens) values
         int ind_svecp = i_qpts * ecmech::nsvp;
         double* stress_svec_p = &(stress_svec_p_array[ind_svecp]);
         double* d_svec_p = &(d_svec_p_array[ind_svecp]);

         // initialize 6x6 2d arrays all to 0
         for (int i = 0; i < ecmech::nsvec; i++) {
            for (int j = 0; j < ecmech::nsvec; j++) {
               ddsdde[(i * ecmech::nsvec) +j] = 0.0;
            }
         }

         for (int i = 0; i < ecmech::ne; i++) {
            eng_int[i] = state_vars[ind_int_eng + i];
         }

         //Here we have the skew portion of our velocity gradient as represented as an
         //axial vector.
         w_vec[0] = 0.5 * (vgrad_view(2, 1, i_qpts) - vgrad_view(1, 2, i_qpts));
         w_vec[1] = 0.5 * (vgrad_view(0, 2, i_qpts) - vgrad_view(2, 0, i_qpts));
         w_vec[2] = 0.5 * (vgrad_view(1, 0, i_qpts) - vgrad_view(0, 1, i_qpts));

         //Really we're looking at the negative of J but this will do...
         double d_mean = -ecmech::onethird * (vgrad_view(0, 0, i_qpts) + vgrad_view(1, 1, i_qpts) + vgrad_view(2, 2, i_qpts));
         //The 1st 6 components are the symmetric deviatoric portion of our velocitu gradient
         //The last value is simply the trace of the deformation rate
         d_svec_p[0] = vgrad_view(0, 0, i_qpts) + d_mean;
         d_svec_p[1] = vgrad_view(1, 1, i_qpts) + d_mean;
         d_svec_p[2] = vgrad_view(2, 2, i_qpts) + d_mean;
         d_svec_p[3] = 0.5 * (vgrad_view(2, 1, i_qpts) + vgrad_view(1, 2, i_qpts));
         d_svec_p[4] = 0.5 * (vgrad_view(2, 0, i_qpts) + vgrad_view(0, 2, i_qpts));
         d_svec_p[5] = 0.5 * (vgrad_view(1, 0, i_qpts) + vgrad_view(0, 1, i_qpts));
         d_svec_p[6] = -3 * d_mean;

         vol_ratio[0] = state_vars[ind_vols];
         vol_ratio[1] = vol_ratio[0] * exp(d_svec_p[ecmech::iSvecP] * dt);
         vol_ratio[3] = vol_ratio[1] - vol_ratio[0];
         vol_ratio[2] = vol_ratio[3] / (dt * 0.5 * (vol_ratio[0] + vol_ratio[1]));

         //std::copy(stress, stress + ecmech::nsvec, stress_svec_p);
         for (int i = 0; i < ecmech::nsvec; i++) {
            stress_svec_p[i] = stress[i];
         }

         double stress_mean = -ecmech::onethird * (stress[0] + stress[1] + stress[2]);
         stress_svec_p[0] += stress_mean;
         stress_svec_p[1] += stress_mean;
         stress_svec_p[2] += stress_mean;
         stress_svec_p[ecmech::iSvecP] = stress_mean;
      });//end of qpt loop
   }//end setup_data_cuda

#endif
}

//Here we're going to initialize all of the data that's going inside of
//of our material update function call.
//This function is used to initialize the data originally
void init_data(Accelerator accel, const double* ori, const ecmech::matModelBase* mat_model_base,
               const int nqpts, const int num_hardness,
               const int num_slip, const int ind_gdot,
               const int state_var_vdim, double* state_vars){
   //We want to store the initial history vector and use information from it to instantiate
   //everything else.
   //When we pass this to the for loop we'll want to use just the raw double array
   #ifdef __cuda_host_only__
   std::vector<double> histInit_vec;
   {
      std::vector<std::string> names;
      std::vector<bool>        plot;
      std::vector<bool>        state;
      mat_model_base->getHistInfo(names, histInit_vec, plot, state);
   }

   const int vdim = state_var_vdim;

   #if defined(RAJA_ENABLE_OPENMP)
   if (accel == Accelerator::OPENMP) {
      init_data_openmp(ori, histInit_vec, nqpts, num_hardness, ind_gdot, num_slip, vdim, state_vars);
   }
   #endif
   if (accel == Accelerator::CPU || accel == Accelerator::CUDA) {
      init_data_cpu(ori, histInit_vec, nqpts, num_hardness, ind_gdot, num_slip, vdim, state_vars);
   }
   #endif
}//end of init_data

//This sets the macroscopic vgrad to be purely deviatoric and behaving as a tension test in the
//z direction. More interesting vgrads could be created just as easily as well where we also have some
//spin terms as well. We could also create a case where there is some sort of spin term as well.
void setup_vgrad(double* vgrad, const int nqpts){
   //vgrad is kinda a pain to deal with as a raw 1d array, so we're
   //going to just use a RAJA view here. The data is taken to be in col. major format.
   //It might be nice to eventually create a type alias for the below or
   //maybe something like it.
   const int DIM = 3;
   std::array<RAJA::idx_t, DIM> perm {{ 2, 1, 0 } };
   RAJA::Layout<DIM> layout = RAJA::make_permuted_layout({{ ecmech::ndim, ecmech::ndim, nqpts } }, perm);
   RAJA::View<double, RAJA::Layout<DIM, RAJA::Index_type, 0> > vgrad_view(vgrad, layout);

   RAJA::RangeSegment default_range(0, nqpts);

   RAJA::forall<RAJA::loop_exec>(default_range, [ = ](int i) {
      // for (int i = 0; i < nqpts; i++) {
      vgrad_view(0, 0, i) = -0.5;
      vgrad_view(0, 1, i) = 0.0;
      vgrad_view(0, 2, i) = 0.0;

      vgrad_view(1, 0, i) = 0.0;
      vgrad_view(1, 1, i) = -0.5;
      vgrad_view(1, 2, i) = 0.0;

      vgrad_view(2, 0, i) = 0.0;
      vgrad_view(2, 1, i) = 0.0;
      vgrad_view(2, 2, i) = 1.0;
   });//end of qpt loop
}//end of setup_vgrad

//This function/kernel is used to set-up the problem at each time step
void setup_data(Accelerator accel, const int nqpts, const int nstatev,
                const double dt, const double* vel_grad_array,
                const double* stress_array, const double* state_vars_array,
                double* stress_svec_p_array, double* d_svec_p_array,
                double* w_vec_array, double* ddsdde_array,
                double* vol_ratio_array, double* eng_int_array){
   #if defined(RAJA_ENABLE_OPENMP)
   if (accel == Accelerator::OPENMP) {
      setup_data_openmp(nqpts, nstatev, dt, vel_grad_array, stress_array, state_vars_array,
                        stress_svec_p_array, d_svec_p_array, w_vec_array, ddsdde_array,
                        vol_ratio_array, eng_int_array);
   }
   #endif
   #if defined(RAJA_ENABLE_CUDA)
   if (accel == Accelerator::CUDA) {
      setup_data_cuda(nqpts, nstatev, dt, vel_grad_array, stress_array, state_vars_array,
                      stress_svec_p_array, d_svec_p_array, w_vec_array, ddsdde_array,
                      vol_ratio_array, eng_int_array);
   }
   #endif
   if (accel == Accelerator::CPU) {
      setup_data_cpu(nqpts, nstatev, dt, vel_grad_array, stress_array, state_vars_array,
                     stress_svec_p_array, d_svec_p_array, w_vec_array, ddsdde_array,
                     vol_ratio_array, eng_int_array);
   }
}//end setup_data
