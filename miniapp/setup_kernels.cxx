/**
 * @file setup_kernels.cxx
 * @brief Implementation of `init_data`, `setup_velocity_grad`, and `setup_data`
 * (declared in `setup_kernels.h`).
 */

#include "setup_kernels.h"

#include "ECMech_evptnWrap.h"
#include <math.h>
// See setup_kernels.h for the full @brief/@param doc.
void init_data(const std::vector<double>& ori_vec, const ecmech::matModelBase* mat_model_base,
               const int nqpts, const int num_hardness,
               const int num_slip, const int ind_gdot,
               const int state_var_vdim, double* state_vars){
   // We want to store the initial history vector and use information from it to instantiate
   // everything else.
   // When we pass this to the for loop we'll want to use just the raw double array

   std::vector<double> histInit_vec;
   {
      std::vector<std::string> names;
      std::vector<bool> plot;
      std::vector<bool> state;
      mat_model_base->getHistInfo(names, histInit_vec, plot, state);
   }

   std::cout << " hist_init values: " << std::endl;
   for (size_t i = 0; i < histInit_vec.size(); i++) {
      std::cout << histInit_vec[i] << " ";
   }
   std::cout << std::endl;

   // histInit_vec/ori_vec are plain host std::vectors. Without the CHAI-backed RAJA
   // port suite, a raw pointer into them is fine for the snls::forall lambda below
   // (host-only build). With it, GPU execution needs its own device-visible copies, so
   // we stage both into chai::ManagedArrays on the host and then fetch each array's
   // pointer in whatever execution space snls::Device is currently configured for.
#if !defined(SNLS_RAJA_PORT_SUITE)
   const auto histInit_data = histInit_vec.data();
   const auto ori = ori_vec.data();
#else
   auto mm = snls::memoryManager::getInstance();
   auto mvec = mm.allocManagedArray<double>(histInit_vec.size());
   auto mvec_data = mvec.data(chai::ExecutionSpace::CPU);
   for (size_t i = 0; i < histInit_vec.size(); i++ ) {
      mvec_data[i] = histInit_vec[i];
   }
   auto ovec = mm.allocManagedArray<double>(ori_vec.size());
   auto ovec_data = ovec.data(chai::ExecutionSpace::CPU);
   for (size_t i = 0; i < ori_vec.size(); i++ ) {
      ovec_data[i] = ori_vec[i];
   }
   const auto histInit_data = mvec.data(snls::Device::GetInstance().GetCHAIES());
   const auto ori = ovec.data(snls::Device::GetInstance().GetCHAIES());
#endif

   const int vdim = state_var_vdim;

   const int ind_dp_eff = ecmech::evptn::iHistA_shrateEff;
   const int ind_eql_pl_strain = ecmech::evptn::iHistA_shrEff;
   const int ind_flow_stress = ecmech::evptn::iHistA_flowStr;
   const int ind_num_evals = ecmech::evptn::iHistA_nFEval;
   const int ind_dev_elas_strain = ecmech::evptn::iHistLbE;
   const int ind_quats = ecmech::evptn::iHistLbQ;
   const int ind_hardness = ecmech::evptn::iHistLbH;
   // The number of vols -> we actually only need to save the previous time step value
   // instead of all 4 values used in the evalModel. The rest can be calculated from
   // this value.
   const int num_vols = 1;
   // These two slots are *not* part of the model's own history layout (nothing at or
   // past ind_gdot + num_slip comes from histInit_data/getHistInfo) -- they're this
   // miniapp's own bookkeeping, appended right after the model's slip-rate block. See
   // this file's @file doc / init_data's doc in setup_kernels.h.
   const int ind_vols = ind_gdot + num_slip;
   const int ind_int_eng = ind_vols + num_vols;

   // We're going to use RAJA here to initialize everything all at once
   // We should be able to just use OpenMP here. Since, we can assume everything is on the
   // host originally. We'll later want to migrate everything over to the GPU if we're
   // running things on there.

   snls::forall(0, nqpts, [=]
      __ecmech_hdev__
      (int i)
   {
      int ind = i * vdim;
      int ind_ori = i * ecmech::qdim;

      state_vars[ind + ind_dp_eff] = histInit_data[ind_dp_eff];
      state_vars[ind + ind_eql_pl_strain] = histInit_data[ind_eql_pl_strain];
      state_vars[ind + ind_flow_stress] = histInit_data[ind_flow_stress];
      state_vars[ind + ind_num_evals] = histInit_data[ind_num_evals];
      // Normally, we might have this assigned as 0 but that may not always
      // be the case.
      state_vars[ind + ind_vols] = 1.0;
      // Due to the smallness of these loops, I would expect several of these to be
      // completely unrolled by the compiler.
      for (int j = 0; j < ecmech::qdim; j++) {
         state_vars[ind + ind_quats + j] = ori[ind_ori + j];
      }

      for (int j = 0; j < ecmech::ntvec; j++) {
         state_vars[ind + ind_dev_elas_strain + j] = histInit_data[ind_dev_elas_strain + j];
      }

      for (int j = 0; j < num_slip; j++) {
         state_vars[ind + ind_gdot + j] = histInit_data[ind_gdot + j];
      }

      for (int j = 0; j < num_hardness; j++) {
         state_vars[ind + ind_hardness + j] = histInit_data[ind_hardness + j];
      }

      // Normally, we might have this assigned as 0 but that may not always
      // be the case.
      for (int j = 0; j < ecmech::ne; j++) {
         state_vars[ind + ind_int_eng + j] = 0.0;
      }
   });
} // end of init_data

// See setup_kernels.h for the full @brief/@param doc. This function itself just
// broadcasts whatever 3x3 matrix it's given to every point; it's the option file's
// velocity-gradient entry (see orientation_evolution.cxx) that determines whether the
// resulting run is a deviatoric uniaxial-tension test, includes spin, etc. -- the
// default option files use a purely deviatoric z-direction tension velocity gradient.
void setup_velocity_grad(const std::vector<double>& velocity_grad_input, double* const velocity_grad, const int nqpts){
   // velocity_grad is stored as a RAJA view with a permuted (2,1,0) layout so that,
   // despite indexing as (row, col, point), the fastest-varying dimension in memory is
   // the point index -- i.e. all nqpts copies of a given (row, col) entry are
   // contiguous, which is the layout setup_data's own velocity_grad_view expects.

#if !defined(SNLS_RAJA_PORT_SUITE)
   const auto velocity_grad_data = velocity_grad_input.data();
#else
   auto mm = snls::memoryManager::getInstance();
   auto mvec = mm.allocManagedArray<double>(velocity_grad_input.size());
   auto mvec_data = mvec.data(chai::ExecutionSpace::CPU);
   for (size_t i = 0; i < velocity_grad_input.size(); i++ ) {
      mvec_data[i] = velocity_grad_input[i];
   }
   const auto velocity_grad_data = mvec.data(snls::Device::GetInstance().GetCHAIES());
#endif

   const int DIM = 3;
   std::array<RAJA::idx_t, DIM> perm { { 2, 1, 0 } };
   RAJA::Layout<DIM> layout = RAJA::make_permuted_layout({ { ecmech::ndim, ecmech::ndim, nqpts } }, perm);
   RAJA::View<double, RAJA::Layout<DIM, RAJA::Index_type, 0> > velocity_grad_view(velocity_grad, layout);

   RAJA::RangeSegment default_range(0, nqpts);

   snls::forall(0, nqpts, [=]
      __ecmech_hdev__
      (int i)
   {
      velocity_grad_view(0, 0, i) = velocity_grad_data[0];
      velocity_grad_view(0, 1, i) = velocity_grad_data[1];
      velocity_grad_view(0, 2, i) = velocity_grad_data[2];

      velocity_grad_view(1, 0, i) = velocity_grad_data[3];
      velocity_grad_view(1, 1, i) = velocity_grad_data[4];
      velocity_grad_view(1, 2, i) = velocity_grad_data[5];

      velocity_grad_view(2, 0, i) = velocity_grad_data[6];
      velocity_grad_view(2, 1, i) = velocity_grad_data[7];
      velocity_grad_view(2, 2, i) = velocity_grad_data[8];
   }); // end of qpt loop
} // end of setup_velocity_grad

// See setup_kernels.h for the full @brief/@param doc.
void setup_data(const int nqpts, const int nstatev,
                const double dt, const double* vel_grad_array,
                const double* cauchy_stress_array, const double* state_vars_array,
                double* cauchy_stress_d6p_array, double* def_rate_d6v_array,
                double* spin_vec_array, double* ddsdde_array,
                double* rel_vol_ratios_array, double* internal_energy_array,
                double* tkelv_array){
   // Same permuted-layout view as setup_velocity_grad, just read-only here.

   // Mirror of init_data/retrieve_data's layout: the volume-ratio and internal-energy
   // bookkeeping slots are the last 1 + ecmech::ne entries of the nstatev-wide record.
   const int ind_int_eng = nstatev - ecmech::ne;
   const int ind_vols = ind_int_eng - 1;

   const int DIM = 3;
   std::array<RAJA::idx_t, DIM> perm { { 2, 1, 0 } };
   RAJA::Layout<DIM> layout = RAJA::make_permuted_layout({ { ecmech::ndim, ecmech::ndim, nqpts } }, perm);
   RAJA::View<const double, RAJA::Layout<DIM, RAJA::Index_type, 0> > velocity_grad_view(vel_grad_array, layout);

   snls::forall(0, nqpts, [=]
      __ecmech_hdev__
      (int i_qpts)
   {
      // Might want to eventually set these all up using RAJA views. It might simplify
      // things later on.
      // These are our inputs
      const double* state_vars = &(state_vars_array[i_qpts * nstatev]);
      const double* cauchy_stress = &(cauchy_stress_array[i_qpts * ecmech::nsvec]);
      // Here is all of our ouputs
      double* ddsdde = &(ddsdde_array[i_qpts * ecmech::nsvec * ecmech::nsvec]);
      double* internal_energy = &(internal_energy_array[i_qpts * ecmech::ne]);
      double* spin_vec = &(spin_vec_array[i_qpts * ecmech::nwvec]);
      double* rel_vol_ratios = &(rel_vol_ratios_array[i_qpts * ecmech::nvr]);
      // A few variables are set up as the 6-vec deviatoric + tr(tens) values
      int ind_svecp = i_qpts * ecmech::nsvp;
      double* cauchy_stress_d6p = &(cauchy_stress_d6p_array[ind_svecp]);
      double* def_rate_d6p = &(def_rate_d6v_array[ind_svecp]);

      tkelv_array[i_qpts] = 300.;

      // initialize 6x6 2d arrays all to 0
      for (int i = 0; i < ecmech::nsvec; i++) {
         for (int j = 0; j < ecmech::nsvec; j++) {
            ddsdde[(i * ecmech::nsvec) +j] = 0.0;
         }
      }

      for (int i = 0; i < ecmech::ne; i++) {
         internal_energy[i] = state_vars[ind_int_eng + i];
      }

      // Here we have the skew portion of our velocity gradient as represented as an
      // axial vector.
      spin_vec[0] = 0.5 * (velocity_grad_view(2, 1, i_qpts) - velocity_grad_view(1, 2, i_qpts));
      spin_vec[1] = 0.5 * (velocity_grad_view(0, 2, i_qpts) - velocity_grad_view(2, 0, i_qpts));
      spin_vec[2] = 0.5 * (velocity_grad_view(1, 0, i_qpts) - velocity_grad_view(0, 1, i_qpts));

      // def_rate_mean = -trace(L)/3 (the sign matches the svecp/pressure convention in
      // ECMech_util.h: adding it onto the diagonal below removes the trace, and
      // negating it back out at the end restores it as the volumetric-rate slot).
      double def_rate_mean = -ecmech::onethird * (velocity_grad_view(0, 0, i_qpts) + velocity_grad_view(1, 1, i_qpts) + velocity_grad_view(2, 2, i_qpts));
      // The 1st 6 components are the symmetric, trace-removed (deviatoric) portion of
      // the velocity gradient -- diagonal entries directly (already symmetric), shear
      // entries symmetrized as 0.5*(L_ij + L_ji). The 7th (index ecmech::iSvecP) is
      // trace(L) = trace(D) itself (the volumetric rate), matching the def_rate_d6vV
      // convention documented on matModelBase::getResponseECM.
      def_rate_d6p[0] = velocity_grad_view(0, 0, i_qpts) + def_rate_mean;
      def_rate_d6p[1] = velocity_grad_view(1, 1, i_qpts) + def_rate_mean;
      def_rate_d6p[2] = velocity_grad_view(2, 2, i_qpts) + def_rate_mean;
      def_rate_d6p[3] = 0.5 * (velocity_grad_view(2, 1, i_qpts) + velocity_grad_view(1, 2, i_qpts));
      def_rate_d6p[4] = 0.5 * (velocity_grad_view(2, 0, i_qpts) + velocity_grad_view(0, 2, i_qpts));
      def_rate_d6p[5] = 0.5 * (velocity_grad_view(1, 0, i_qpts) + velocity_grad_view(0, 1, i_qpts));
      def_rate_d6p[6] = -3 * def_rate_mean;
      // Integrate the volumetric rate over dt (J_{n+1} = J_n * exp(trace(D)*dt)) to get
      // this step's relative volume, then derive the rate/delta slots -- see
      // ECMech_const.h's `nvr` doc for the [rel_vol_n, rel_vol_n+1, rate, delta] layout.
      rel_vol_ratios[0] = state_vars[ind_vols];
      rel_vol_ratios[1] = rel_vol_ratios[0] * exp(def_rate_d6p[ecmech::iSvecP] * dt);
      rel_vol_ratios[3] = rel_vol_ratios[1] - rel_vol_ratios[0];
      rel_vol_ratios[2] = rel_vol_ratios[3] / (dt * 0.5 * (rel_vol_ratios[0] + rel_vol_ratios[1]));

      // Convert the previous step's Cauchy stress from plain Voigt form to ExaCMech's
      // deviatoric + pressure (svecp) form -- the inverse of what retrieve_data
      // (retrieve_kernels.cxx) does going back out: subtract the mean normal stress
      // from the three normal components to make them traceless, and stash it (negated,
      // per the svecp/pressure convention in ECMech_util.h) in the iSvecP slot.
      for (int i = 0; i < ecmech::nsvec; i++) {
         cauchy_stress_d6p[i] = cauchy_stress[i];
      }

      double stress_mean = -ecmech::onethird * (cauchy_stress[0] + cauchy_stress[1] + cauchy_stress[2]);
      cauchy_stress_d6p[0] += stress_mean;
      cauchy_stress_d6p[1] += stress_mean;
      cauchy_stress_d6p[2] += stress_mean;
      cauchy_stress_d6p[ecmech::iSvecP] = stress_mean;
   }); // end of qpt loop
} // end setup_data

