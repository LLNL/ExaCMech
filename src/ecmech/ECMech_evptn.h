// -*-c++-*-

#ifndef ECMECH_EVPTN_H
#define ECMECH_EVPTN_H

#include <cassert>

#include "ECMech_core.h"
#include "ECMech_util.h"
#include "ECMech_elastic.h"
#include "ECMech_eosSimple.h"

#include "SNLS_lup_solve.h"
#include "SNLS_TrDLDenseG.h"
#include "SNLS_HybrdTrDLDenseG.h"

#include "RAJA/RAJA.hpp"

namespace ecmech {
   namespace evptn {
      const int numHistAux = 4; // effective shearing rate, accumulated shear, flow strength, nFEval
      //
      const int iHistLbA = 0;
      const int iHistA_shrateEff = iHistLbA + 0;
      const int iHistA_shrEff = iHistLbA + 1;
      const int iHistA_flowStr = iHistLbA + 2;
      const int iHistA_nFEval = iHistLbA + 3;
      const int iHistLbE = numHistAux;
      const int iHistLbQ = numHistAux + ecmech::ntvec;
      const int iHistLbH = numHistAux + ecmech::ntvec + ecmech::qdim;

      /*
       * just a container for a traits
       */
      template<class SlipGeom, class Kinetics, class ThermoElastN, class EosModel>
      class NumHist
      {
         public:
            // see n_rsv_matmod in F90 code
            static const int iHistLbGdot = iHistLbH + Kinetics::nH;
            static const int numHist = iHistLbH + Kinetics::nH + SlipGeom::nslip;
      }; // NumHist

      template<class ThermoElastN>
      class EvptnLatticeStrainProblem
      {
         public:
         static const int nDimSys = ecmech::ntvec;

         __ecmech_hdev__
         EvptnLatticeStrainProblem(const ThermoElastN& thermoElastN,
                                 const double dt,
                                 const double detV, 
                                 const double eVref, 
                                 const double p_EOS, 
                                 const double tK,
                                 const double* const e_vecd_n)
            : m_thermo_elast_n(thermoElastN),
            m_dt(dt), m_det_vol(detV), m_elast_vol_ref(eVref),
            m_pressure_eos(p_EOS), m_temp_k(tK),
            m_elast_dev_vec_n(e_vecd_n),
            m_inv_dt(1.0 / dt),
            m_inv_det_vol(1.0 / m_det_vol),
            m_a_vol(pow(m_det_vol, onethird)),
            m_inv_a_vol(1.0 / m_a_vol)
            {};

         __ecmech_hdev__
         ~EvptnLatticeStrainProblem() = default;

         // used to be elastNEtoT
         __ecmech_hdev__
         inline
         void elas_strain_to_kirchoff_stress(double* const kirchoff_stress, // nsvec
                                             const double* const elas_dev_vec // ntvec
                                          ) const
         {
            //// do not need to use elaw_T_BT here as T and BT are the same
            //
            // specialize to cem%l_lin_lnsd
            // CALL elawn_T(s_meas, e_vecd_f, crys%elas, tK, .TRUE., a_V, &
            // & p_EOS, eVref, crys%i_eos_model, crys%eos_const &
            // &)
            double elas_dev_vol_vec[ecmech::nsvec];
            vecsVxa<ntvec>(elas_dev_vol_vec, m_inv_a_vol, elas_dev_vec);
            //// tr_Ee = three * DLOG(a_V%r)
            //// CALL trace_to_vecds_s(s_meas%Ee_vecds(SVEC), tr_Ee)
            elas_dev_vol_vec[iSvecS] = sqr3 * log(m_a_vol); // could go into constructor
            //
            //// Kirchhoff stress from elas_dev_vol_vec
            // CALL elawn_lin_op(s_meas%T_vecds, s_meas%Ee_vecds, cem, tK, &
            // & p_EOS, eVref, i_eos_model, eos_const)
            m_thermo_elast_n.eval(kirchoff_stress, elas_dev_vol_vec, m_temp_k, m_pressure_eos, m_elast_vol_ref);
         }

         // used to be elastNEtoC
         __ecmech_hdev__
         inline
         void elas_strain_to_cauchy_stress(double* const cauchy, // nsvec
                                          const double* const e_vecd_f // ntvec
                                          ) const
         {
            double kirchoff[ecmech::nsvec];
            this->elas_strain_to_kirchoff_stress(kirchoff, e_vecd_f);
            m_thermo_elast_n.getCauchy(cauchy, kirchoff, m_inv_det_vol);
         }

         __ecmech_hdev__
         template<bool calc_strain_rate = false>
         inline
         void get_elas_strain_state(double* const elas_delta_dev_vec,
                                    double* const elas_dt_dev_vec,
                                    const double* const x) const
         {
            //////////////////////////////
            // PULL VALUES out of x, with scalings
            //
            // double edot_vecd[ecmech::ntvec];
            vecsVxa<ntvec>(elas_dt_dev_vec, ecmech::e_scale, x); // elas_dt_dev_vec is now the delta, _not_ yet edot_vecd
            // e_vecd_f is end-of-step
            // double e_vecd_f[ntvec];
            vecsVapb<ntvec>(elas_delta_dev_vec, elas_dt_dev_vec, m_elast_dev_vec_n);
            if constexpr(calc_strain_rate) {
               vecsVsa<ntvec>(elas_dt_dev_vec, m_inv_dt); // _now_ elas_dt_dev_vec has dt contributions
            }
         }

         /*
            * NOTES :
            * () should be equivalent to what happens in get_elas_strain_state<false>
            * () not necessarily safe if e_vecd is the same memory as _e_vecd_n or quat is the same as _Cn_quat
            */
         __ecmech_hdev__
         inline
         void stateFromX(double* const elast_dev_vec,
                        const double* const x) const
         {
            double elast_dev_vec_delta[ecmech::ntvec] = {};
            this->get_elas_strain_state(elast_dev_vec, elast_dev_vec_delta, x);
         }

         __ecmech_hdev__
         inline
         void get_elas_strain_residual(double* const residual,
                                       const double epsdot_scale_inv,
                                       const double* const elas_dt_dev_vec,
                                       const double* const plastic_def_rate_dev_vec,
                                       const double* const def_rate_dev_vec_lattice) const
         {
            for (size_t iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
               residual[m_ind_sub_elas + iTvec] = epsdot_scale_inv * ( // SCALING
                  m_inv_a_vol * elas_dt_dev_vec[iTvec] + plastic_def_rate_dev_vec[iTvec] - def_rate_dev_vec_lattice[iTvec]);
            }
         }

         __ecmech_hdev__
         template<size_t JAC_SIZE>
         inline
         void get_deriv_elast_strain_wrt_elast_strain(double* const jacobian,
                                                      const double* const dDp_hat_delast_strain) const
         {
            RAJA::View<double, RAJA::Layout<2>> jacob_ee(jacobian, JAC_SIZE, JAC_SIZE);
            // dislocation plasticity;
            // first contribution; overwrite
            for (size_t jTvec = 0; jTvec < ecmech::ntvec; ++jTvec) {
               for (size_t iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                  jacob_ee(iTvec, jTvec) = dDp_hat_delast_strain[ECMECH_NN_INDX(iTvec, jTvec, ecmech::ntvec)];
               }
            }
            // elastic rate
            //
            {
               const double adti = m_inv_a_vol * m_inv_dt;
               for (size_t iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                  jacob_ee(iTvec, iTvec) += adti;
               }
            }
         } 

         __ecmech_hdev__
         template<size_t JAC_SIZE, size_t ind_sub_r>
         inline
         void get_deriv_omega_wrt_elast_strain(double* const jacobian,
                                               const double* const elas_dt_dev_vec,
                                               const double elast_elast_factor,
                                               const double* const dWp_hat_delast_strain,
                                               const double* const A_e_M35) const
         {
            // d(B_xi)/d(e_vecds_f)
            //
            RAJA::View<double, RAJA::Layout<2>> jacob_re(jacobian, JAC_SIZE, JAC_SIZE);

            double A_edot_M35[ecmech::nwvec * ecmech::ntvec];
            M35_d_AAoB_dA(A_edot_M35, elas_dt_dev_vec);

            double dt_ee_fac = m_dt * elast_elast_factor;

            for (size_t iWvec = 0; iWvec < ecmech::nwvec; ++iWvec) {
               for (size_t jTvec = 0; jTvec < ecmech::ntvec; ++jTvec) {
                  int ijWT = ECMECH_NM_INDX(iWvec, jTvec, ecmech::nwvec, ecmech::ntvec);
                  jacob_re(iWvec + ind_sub_r, jTvec) =
                     m_dt * dWp_hat_delast_strain[ijWT] - dt_ee_fac * (A_e_M35[ijWT] * m_inv_dt - A_edot_M35[ijWT]);
               }
            }         
         }

         __ecmech_hdev__
         template<size_t JAC_SIZE, size_t ind_sub_h, size_t num_hard, size_t num_slip>
         inline
         void get_deriv_hardening_wrt_elast_strain(double* const jacobian,
                                                   const double* const dhard_dgdot,
                                                   const double* const dgdot_delast_strain) const
         {
            // d(B_h) / d(e)
            // jacob_he = dt * (dhdot/dgdot)(dgdot/de)
            // nh x ntvec matrix
            double dhdot_delas_strain[num_hard * ecmech::ntvec];
            vecsMABT<num_hard, ecmech::ntvec, num_slip>(dhdot_delas_strain, dhard_dgdot, dgdot_delast_strain);
            RAJA::View<double, RAJA::Layout<2>> jacob_he(jacobian, JAC_SIZE, JAC_SIZE);
            for (size_t iH = 0; iH < num_hard; ++iH) {
               for (size_t jE = 0; jE < ecmech::ntvec; ++jE) {
                  // could also make dhdot_delas_strain into a RAJA view, but not really needed
                  jacob_he(iH + ind_sub_h, jE) = -m_dt * dhdot_delas_strain[ECMECH_NM_INDX(iH, jE, num_hard, ecmech::ntvec) ];
               }
            }       
         }

         public:
            static constexpr size_t m_ind_sub_elas = 0; // ntvec end_point
            const ThermoElastN& m_thermo_elast_n;
            const double m_dt, m_det_vol, m_elast_vol_ref;
            const double m_pressure_eos, m_temp_k;
            const double* const m_elast_dev_vec_n;
            const double m_inv_dt, m_inv_det_vol, m_a_vol, m_inv_a_vol;
      };

      template <size_t ind_sub_r=ecmech::ntvec>
      class EvptnLatticeRotationProblem {
         public:
         static constexpr size_t nDimSys = ecmech::nwvec;

         public:
         __ecmech_hdev__
         EvptnLatticeRotationProblem(const double dt,
                                    const double* const Cn_quat)
            : m_dt(dt), m_xtal_ori_quat_n(Cn_quat) {};
         
         __ecmech_hdev__
         ~EvptnLatticeRotationProblem() = default;

         __ecmech_hdev__
         inline
         void get_rotation_state(double* const delta_omega,
                                 double* const xtal_rmat,
                                 double* const xtal_rot_mat5,
                                 const double* const x) const
         {
            vecsVxa<ecmech::nwvec>(delta_omega, ecmech::r_scale, &(x[ind_sub_r]));
            //
            // not done in EvpC :
            // CALL exp_map_cpvec(A, xi_f)
            // CALL get_c(c, A, C_n)
            //
            double xtal_ori_quat_delta[ecmech::qdim];
            double xtal_ori_quat_n1[ecmech::qdim];

            emap_to_quat(xtal_ori_quat_delta, delta_omega);
            get_c_quat(xtal_ori_quat_n1, xtal_ori_quat_delta, m_xtal_ori_quat_n);
            quat_to_tensor(xtal_rmat, xtal_ori_quat_n1);
            get_rot_mat_vecd(xtal_rot_mat5, xtal_rmat);
         }

         /*
         * NOTES :
         * () should be equivalent to what happens in get_elas_strain_state<false>
         * () not necessarily safe if e_vecd is the same memory as _e_vecd_n or quat is the same as _Cn_quat
         */
         __ecmech_hdev__
         inline
         void stateFromX(double* const xtal_ori_quat,
                        const double* const x) const
         {
            double delta_omega[ecmech::nwvec];
            double xtal_ori_quat_delta[ecmech::qdim];
            vecsVxa<ecmech::nwvec>(delta_omega, ecmech::r_scale, &(x[ind_sub_r]));
            emap_to_quat(xtal_ori_quat_delta, delta_omega);
            get_c_quat(xtal_ori_quat, xtal_ori_quat_delta, m_xtal_ori_quat_n);
         }

         __ecmech_hdev__
         inline
         void get_omega_residual(double* const residual,
                                 const double rot_incr_scale_inv,
                                 const double ee_fac,
                                 const double* const delta_omega,
                                 const double* const spin_vec_lat,
                                 const double* const plastic_spin_vec,
                                 const double* const ee_spin_vec) const
         {
            // RESIDUAL B_omega
            for (int iWvec = 0; iWvec < ecmech::nwvec; ++iWvec) {
               residual[ind_sub_r + iWvec] = rot_incr_scale_inv * // SCALING
                                          (delta_omega[iWvec] - m_dt * (spin_vec_lat[iWvec]
                                          - plastic_spin_vec[iWvec] + ee_fac * ee_spin_vec[iWvec]));
            }
         }

         __ecmech_hdev__
         template<size_t JAC_SIZE>
         inline
         void get_deriv_omega_wrt_omega(double* const jacobian,
                                        const double* const dspin_samp_domega) const
         {
            // d(B_xi)/d(xi_f)
            RAJA::View<double, RAJA::Layout<2>> jacob_rr(jacobian, JAC_SIZE, JAC_SIZE);
            for (int iWvec = 0; iWvec < ecmech::nwvec; ++iWvec) {
               for (int jWvec = 0; jWvec < ecmech::nwvec; ++jWvec) {
                  int ijWW = ECMECH_NN_INDX(iWvec, jWvec, ecmech::nwvec);
                  jacob_rr(iWvec + ind_sub_r, jWvec + ind_sub_r) = -m_dt * dspin_samp_domega[ijWW];
               }
               jacob_rr(iWvec + ind_sub_r, iWvec + ind_sub_r) += one;
            }
         }
         __ecmech_hdev__
         template<size_t JAC_SIZE>
         inline
         void get_deriv_elast_strain_wrt_omega(double* const jacobian,
                                             const double* const ddef_rate_samp_domega) const
         {
            // d(B_S)/d(xi_f)
            //
            // jacob_er = -dDsm_dxi(:,:)
            RAJA::View<double, RAJA::Layout<2> > jacob_er(jacobian, JAC_SIZE, JAC_SIZE);
            for (int jWvec = 0; jWvec<ecmech::nwvec; ++jWvec) {
               for (int iTvec = 0; iTvec<ecmech::ntvec; ++iTvec) {
                  // could also make dDsm_dxi into a RAJA view, but not really needed
                  const size_t ind = ECMECH_NM_INDX(iTvec, jWvec, ecmech::ntvec, ecmech::nwvec);
                  jacob_er(iTvec, jWvec + ind_sub_r) = -ddef_rate_samp_domega[ind];
               }
            }
         }

         /**
            This method does nothing as most crystal plasticity hardening models have no dependency directly or indirectly with the
            lattice rotation.
         */
         __ecmech_hdev__
         inline
         void get_deriv_hardening_wrt_omega(double* const /* jacobian */){}

         private:
         const double m_dt;
         const double* const m_xtal_ori_quat_n;
      };

      __ecmech_hdev__
      inline
      void get_xtal_frame_vel_grad_terms(double* const def_rate_dev_vec_xtal,
                                         double* const spin_vec_xtal,
                                         const double* const def_rate_dev_vec_sample,
                                         const double* const spin_vec_sample,
                                         const double* const xtal_rmat,
                                         const double* const xtal_rot_mat5)
      {
         vecsVMTa<ecmech::ntvec>(def_rate_dev_vec_xtal, xtal_rot_mat5, def_rate_dev_vec_sample);
         vecsVMTa<ecmech::ndim>(spin_vec_xtal, xtal_rmat, spin_vec_sample);
      }

      __ecmech_hdev__
      template<class SlipGeom, class SlipKinetics>
      inline
      void get_slip_rate_terms(double* const dgdot_dtau,
                               double* const plastic_def_rate,
                               double* const plastic_spin_vec,
                               const double* const kirchoff,
                               const double* const kinetic_values,
                               const SlipGeom& slip_geom,
                               const SlipKinetics& slip_kinetics
                               )
      {        
         if constexpr (SlipGeom::nslip > 0) {
            // default initialize everything to 0.0 
            double abs_resolved_shear_stress[SlipGeom::nslip] = {};
            double gdot[SlipGeom::nslip] = {};
            // resolve stress onto slip systems
            // CALL resolve_tau_a_n(crys%tmp4_slp, s_meas%T_vecds, crys)
            //vecsVaTM<ntvec, SlipGeom::nslip>(taua, T_vecds, slipP);
            slip_geom.evalRSS(abs_resolved_shear_stress, kirchoff, slip_geom.getP());
            //
            // CALL plaw_eval(pl_vecd, pl_wvec, gss, crys, tK, ierr)
            // chi values are passed within extended taua array
            slip_kinetics.evalGdots(gdot, dgdot_dtau, nullptr, abs_resolved_shear_stress, kinetic_values);
            
            //
            // CALL sum_slip_def(pl_vecd, pl_wvec, crys%tmp1_slp, crys) ;
            vecsVMa<ntvec, SlipGeom::nslip>(plastic_def_rate, slip_geom.getP(), gdot);
            vecsVMa<nwvec, SlipGeom::nslip>(plastic_spin_vec, slip_geom.getQ(), gdot);
         }
      }

      // This function performs the necessary chain rules to go from the:
      // dgammadot_dRSS -> dDp_hat_dElas_strain
      // dgammadot_dRSS -> dWp_hat_dElas_strain
      // terms used typically in either the Jacobian or material tangent stiffness matrix
      __ecmech_hdev__
      template<class SlipGeom, class ThermoElastN>
      inline
      void get_slip_rate_deriv_terms(double* const dDp_hat_delast_strain,
                                     double* const dWp_hat_delast_strain,
                                     const double* const dgdot_dtau,
                                     const double inv_a_vol,
                                     const SlipGeom& slip_geom,
                                     const ThermoElastN& thermoElastN
                                    )
      {
         if constexpr (SlipGeom::nslip > 0) {
            double dtaua_deps[ ecmech::ntvec * SlipGeom::nslip ];
            thermoElastN.multDTDepsT(dtaua_deps, slip_geom.getP(), inv_a_vol, SlipGeom::nslip);

            double dgdot_deps[ ecmech::ntvec * SlipGeom::nslip ];
            for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
               for (int iSlip = 0; iSlip < SlipGeom::nslip; ++iSlip) {
                  int ijThis = ECMECH_NM_INDX(iTvec, iSlip, ecmech::ntvec, SlipGeom::nslip);
                  dgdot_deps[ijThis] = dgdot_dtau[iSlip] * dtaua_deps[ijThis];
               }
            }
            vecsMABT<ntvec, SlipGeom::nslip>(dDp_hat_delast_strain, slip_geom.getP(), dgdot_deps);
            vecsMABT<nwvec, ntvec, SlipGeom::nslip>(dWp_hat_delast_strain, slip_geom.getQ(), dgdot_deps);
         }
      }

      // The values returned here are usually useful for post-processing and might have application in
      // application codes. However, they are not really state variables in that everything can be
      // calculated post-state variable update. 
      __ecmech_hdev__
      template<class SlipGeom, class SlipKinetics, class Elasticty>
      inline
      void get_slip_contributions(double& pl_disipation_rate,
                                  double& effective_shear_rate,
                                  double* const gdot,
                                  const double inv_det_vol,
                                  const double* const elas_strain,
                                  const double* const kinetic_values,
                                  const SlipGeom& slip_geom,
                                  const SlipKinetics& slip_kinetics,
                                  const Elasticty& elasticity
                                 )
      {
         pl_disipation_rate = 0.0;
         effective_shear_rate = 0.0;
         if constexpr (SlipGeom::nslip > 0) {
            // default initialize everything to 0.0 
            double abs_resolved_shear_stress[SlipGeom::nslip] = {};
            double junk[SlipGeom::nslip] = {};
            double kirchoff[ecmech::nsvec] = {};
            elasticity.elas_strain_to_kirchoff_stress(kirchoff, elas_strain);
            // resolve stress onto slip systems
            slip_geom.evalRSS(abs_resolved_shear_stress, kirchoff, slip_geom.getP());
            slip_kinetics.evalGdots(gdot, junk, junk, abs_resolved_shear_stress, kinetic_values);
            inv_det_vol * vecsyadotb<SlipGeom::nslip>(abs_resolved_shear_stress, gdot);
#if defined(ECMECH_USE_DPEFF)
            vecsVMa<ntvec, SlipGeom::nslip>(plastic_def_rate, slip_geom.getP(), gdot);
            effective_shear_rate = vecd_Deff(plastic_def_rate);
#else
            effective_shear_rate = vecsssumabs<SlipGeom::nslip>(gdot);
#endif
            pl_disipation_rate = inv_det_vol * vecsyadotb<SlipGeom::nslip>(abs_resolved_shear_stress, gdot);
         }
      }

      // These terms are commonly used in the residual and jacobian terms related to omega
      // We should probably try to create better names for these outputs...
      __ecmech_hdev__
      inline
      void elasticity_higher_order_terms(double* const A_e_M35,
                                         double* const ee_spin_vec,
                                         double& ee_fac,
                                         const double inv_a_vol,
                                         const double* const elast_dev_vec,
                                         const double* const elas_dt_dev_vec
                                        )
      {
         // from e edot product term in spin (formerly neglected)
         M35_d_AAoB_dA(A_e_M35, elast_dev_vec);
         vecsVMa<nwvec, ntvec>(ee_spin_vec, A_e_M35, elas_dt_dev_vec);
         ee_fac = onehalf * inv_a_vol * inv_a_vol;
      }

      // Quite a few of these variables could be calculated on the fly
      // However, we already calculate most of them as part of the computeRJ portion of things
      // so just do it again here...
      // Might be able to rework this in a better way at some point...
      __ecmech_hdev__
      template<class ThermoElastN, size_t JAC_SIZE, size_t ind_sub_omega>
      inline
      void get_material_tangent_stiffness(double* const material_tangent,
                                          const double* const jacobian,
                                          const double* const dquat_domega_t, 
                                          const double* const qr5x5_ls,
                                          const double* const quat,
                                          const double* const rmat,
                                          const double* const cauchy_stress,
                                          const double inv_det_vol,
                                          const double inv_a_vol,
                                          const ThermoElastN& thermo_elast_n
                                          )
      {
         // Only concerned with the dCauchy/dDefRate at this point in time
         // and more specifically only considering the deviatoric portion 
         constexpr int nRHS = ecmech::ntvec;

         // mtan_sample_frame
         // = d/dDefRate_sample(Cauchy) = d/dDefRate(Q * alpha * (C_{elas} : lattice_strain)) 
         // Q is the 5x5 rotation oper from crystal to sample 
         // alpha is necessary scaling from Kirchoff to Cauchy
         // C_{elas} is the elasticity tensor (deviatoric contributions)
         // = d/dDefRate_s (Cauchy) = d(Q Cauchy) / dOmega_c * dOmega / dDefRate_s 
         //   + Q * d(Cauchy)/delas_strain * delas_strain / dDefRate_s
         // From our Jacobian we can calculate the dOmega / dDefRate_s and delas_strain / dDefRate_s terms
         // The other ones are either simple to calculate or require some math...
         //
         // dstrainomega_ddef_rate_t => [dlat_strain_ddef_rate_sample; domega_ddef_rate_sample];
         // Initially we set it to be our RHS
         double dstrainomega_ddef_rate_t[ nRHS * JAC_SIZE ] = {}; // transpose for use in SNLS_LUP_SolveX !
         {
            // RHS calculations
            //
            // dstrainomega_ddef_rate_t => dResidual_ddef_rate_sample
            {
               // negatives cancel
               // dstrainomega_ddef_rate_t[0:ind_omega_vec,:] = qr5x5_c2s
               for (int jE = 0; jE < nRHS; ++jE) {
                  for (int iE = 0; iE < ntvec; ++iE) { // ntvec, _not_ nDimSys // iE is same as index 
                     dstrainomega_ddef_rate_t[ECMECH_NM_INDX(jE, iE, nRHS, JAC_SIZE)] = qr5x5_ls[ECMECH_NN_INDX(jE, iE, ntvec)];
                  }
               }
               // If we had hardening terms then we'd add those here as well...
               // dRhard_ddef_rate but typically we can but for most problems safe to treat that
               // set of terms as being = 0
            }
            // Now solve for our dstrain_ddefrate and domega_ddefrate terms
            int err = SNLS_LUP_SolveX<JAC_SIZE>(const_cast<double* const>(jacobian), dstrainomega_ddef_rate_t, nRHS);
            if (err != 0) {
               ECMECH_FAIL(__func__, "error from SNLS_LUP_SolveX");
            }
         }

         {
            double temp_M6[ ecmech::nsvec2 ];
            for (int iTvec = 0; iTvec<ecmech::ntvec; ++iTvec) {
               for (int jTvec = iTvec; jTvec<nRHS; ++jTvec) {
                  const int offset1 = ECMECH_NM_INDX(iTvec, jTvec, nRHS, JAC_SIZE);
                  const int offset2 = ECMECH_NM_INDX(jTvec, iTvec, nRHS, JAC_SIZE);
                  const double tmp = dstrainomega_ddef_rate_t[offset1];
                  dstrainomega_ddef_rate_t[offset1] = dstrainomega_ddef_rate_t[offset2];
                  dstrainomega_ddef_rate_t[offset2] = tmp;
               }
            }
            thermo_elast_n.template multCauchyDif<nRHS, JAC_SIZE>(temp_M6, dstrainomega_ddef_rate_t, inv_det_vol, inv_a_vol);
            // Apply final rotation
            qr6x6_pre_mul<ecmech::nsvec, false>(material_tangent, temp_M6, qr5x5_ls);
         }

         // Calculate the d(QCauchy) / dDefRate_s term now and add that to
         {
            double dcauchy_dquat[ ecmech::ntvec * ecmech::qdim ];
            {
               double drmat_dquat[ ecmech::ndim * ecmech::ndim * ecmech::qdim ];
               d_quat_to_tensor(drmat_dquat, quat);
               double dcauchy_drmat[ ecmech::ntvec * ecmech::ndim * ecmech::ndim ];
               d_rot_mat_vecd_smop(dcauchy_drmat, rmat, cauchy_stress);
               vecsMAB<ntvec, qdim, ndim*ndim>(dcauchy_dquat, dcauchy_drmat, drmat_dquat);
            }

            // We now need to be able to go from our domega_ddef_rate to dquat_ddef_rate
            // dquat_ddef_rate = dquat_domega_t * domega_ddef_rate_t
            double dquat_ddef_rate[ ecmech::qdim * nRHS ];
            for (int ii_I = 0; ii_I < nRHS; ++ii_I) {
               for (int ii_Q = 0; ii_Q < ecmech::qdim; ++ii_Q) {
                  int iiQI = ECMECH_NM_INDX(ii_Q, ii_I, ecmech::qdim, nRHS);
                  dquat_ddef_rate[iiQI] = 0.0;
                  for (int ii_W = 0; ii_W < ecmech::nwvec; ++ii_W) {
                     dquat_ddef_rate[iiQI] +=
                        dquat_domega_t[ECMECH_NM_INDX(ii_W, ii_Q, ecmech::nwvec, ecmech::qdim)] *
                        dstrainomega_ddef_rate_t[ECMECH_NM_INDX(ii_I, ind_sub_omega + ii_W, nRHS, JAC_SIZE)];
                  }
               }
            }

            // Now get the dcauchy_lattice_dI terms by doing ->
            // dcauchy_lattice_dI = d_cauchy_lattice_dquat * dRmat_quat_dI
            double dqcauchy_ddefrate[ ecmech::ntvec * nRHS ];
            vecsMAB<ecmech::ntvec, nRHS, ecmech::qdim>(dcauchy_dquat, dcauchy_dquat, dquat_ddef_rate);
            for (int ii_T = 0; ii_T < ecmech::ntvec; ++ii_T) {
               for (int ii_I = 0; ii_I < nRHS; ++ii_I) {
                  // NOTE : only looping over ntvec, but mtan_sI is nsvec in the first dimension
                  material_tangent[ECMECH_NN_INDX(ii_T, ii_I, ecmech::nsvec)] += dqcauchy_ddefrate[ECMECH_NM_INDX(ii_T, ii_I, ecmech::ntvec, nRHS)];
               }
            }
         }
      }                                         

      template<class SlipGeom, class Kinetics, class ThermoElastN>
      class EvptnUpdstProblem
      {
         public:

            static constexpr int nDimSys = ecmech::ntvec + ecmech::nwvec;

            // constructor
            __ecmech_hdev__
            EvptnUpdstProblem(const SlipGeom& slipGeom,
                              const Kinetics& kinetics,
                              const ThermoElastN& thermoElastN,
                              double dt,
                              double detV, double eVref, double p_EOS, double tK,
                              const double* const h_state,
                              const double* const e_vecd_n,
                              const double* const Cn_quat,
                              const double* const d_vecd_sm, // okay to pass d_vecds_sm, but d_vecd_sm[iSvecS] is not used
                              const double* const w_veccp_sm
                              )
               : _slipGeom(slipGeom),
               _kinetics(kinetics),
               _thermoElastN(thermoElastN),
               _lattice_strain_prob(thermoElastN, dt, detV, eVref, p_EOS, tK, e_vecd_n),
               _lattice_rot_prob(dt, Cn_quat),
               _dt(dt),
               _detV(detV),
               _eVref(eVref),
               _p_EOS(p_EOS),
               _tK(tK),
               _h_state(h_state),
               _e_vecd_n(e_vecd_n),
               _Cn_quat(Cn_quat),
               _d_vecd_sm(d_vecd_sm), // vel_grad_sm%d_vecds
               _w_veccp_sm(w_veccp_sm), // vel_grad_sm%w_veccp
               _mtan_sI(nullptr)
            {
               _dt_ri = 1.0 / _dt;
               _detV_ri = 1.0 / _detV;
               _a_V = pow(detV, onethird);
               _a_V_ri = 1.0 / _a_V;

               _hdn_scale = _kinetics.getVals(_kin_vals, _p_EOS, _tK, _h_state);

               double adots_ref = _kinetics.getFixedRefRate(_kin_vals);
               double eff = vecNorm<ntvec>(_d_vecd_sm); // do not worry about factor of sqrt(twothird)
               if (eff < epsdot_scl_nzeff * adots_ref) {
                  _epsdot_scale_inv = one / adots_ref;
               }
               else {
                  _epsdot_scale_inv = fmin(one / eff, 1e6 * _dt);
               }
               //
               _rotincr_scale_inv = _dt_ri * _epsdot_scale_inv;
            }

            // deconstructor
            __ecmech_hdev__
            ~EvptnUpdstProblem() {}

            __ecmech_hdev__
            inline
            void provideMTan(double* mtan_sI) { _mtan_sI = mtan_sI; }

            __ecmech_hdev__
            inline
            void clearMTan( ) { _mtan_sI = nullptr; }

            __ecmech_hdev__
            inline
            double getDtRi() const { return _dt_ri; }

            __ecmech_hdev__
            inline
            double getShrateEff() const { return 1.0; }

            __ecmech_hdev__
            inline
            double getDisRate() const { return 1.0; }

            __ecmech_hdev__
            inline
            double getHdnScale() const { return _hdn_scale; }

            /*
             * NOTES :
             * () should be equivalent to what happens in computeRJ
             * () not necessarily safe if e_vecd is the same memory as _e_vecd_n or quat is the same as _Cn_quat
             */
            __ecmech_hdev__
            inline
            void stateFromX(double* const e_vecd,
                            double* const quat,
                            const double* const x) {
               _lattice_strain_prob.stateFromX(e_vecd,  &(x[_i_sub_e]));
               _lattice_rot_prob.stateFromX(quat,  x);
            }

            __ecmech_hdev__
            inline
            void elastNEtoC(double* const C_vecds, // nsvec
                            const double* const e_vecd_f // ntvec
                            ) const {
               _lattice_strain_prob.elas_strain_to_cauchy_stress(C_vecds, e_vecd_f);
            }

            __ecmech_hdev__
            inline
            bool computeRJ(double* const resid,
                           double* const Jacobian,
                           const double* const x) {
               bool doComputeJ = (Jacobian != nullptr);

               if (doComputeJ) {
                  // zero the Jacobian so that do not need to worry about zero
                  // entries in the midst of other things later
                  //
                  for (int ijJ = 0; ijJ<_nXnDim; ++ijJ) {
                     Jacobian[ijJ] = 0.0;
                  }
               }
               //
               for (int iR = 0; iR<nDimSys; ++iR) {
                  resid[iR] = 0.0;
               }

               //////////////////////////////
               // PULL VALUES out of x, with scalings
               //
               double edot_vecd[ecmech::ntvec];
               vecsVxa<ntvec>(edot_vecd, ecmech::e_scale, &(x[_i_sub_e]) ); // edot_vecd is now the delta, _not_ yet edot_vecd
               // e_vecd_f is end-of-step
               double e_vecd_f[ntvec];
               vecsVapb<ntvec>(e_vecd_f, edot_vecd, _e_vecd_n);
               vecsVsa<ntvec>(edot_vecd, _dt_ri); // _now_ edot_vecd has edot_vecd
               //
               double xi_f[nwvec];
               vecsVxa<nwvec>(xi_f, ecmech::r_scale, &(x[_i_sub_r]) );
               //
               // not done in EvpC :
               // CALL exp_map_cpvec(A, xi_f)
               // CALL get_c(c, A, C_n)
               //
               double A_quat[ecmech::qdim];
               emap_to_quat(A_quat, xi_f);
               //
               double C_quat[ecmech::qdim];
               get_c_quat(C_quat, A_quat, _Cn_quat);
               //
               double C_matx[ecmech::ndim * ecmech::ndim];
               quat_to_tensor(C_matx, C_quat);
               //
               double qr5x5_ls[ecmech::ntvec * ecmech::ntvec];
               get_rot_mat_vecd(qr5x5_ls, C_matx);
               //
               // CALL matt_x_vec_5(qr5x5_ls, vel_grad_sm%d_vecds(1:TVEC), d_vecd_lat)
               double d_vecd_lat[ecmech::ntvec];
               // vecsVMTa<ntvec>(d_vecd_lat, qr5x5_ls, _d_vecd_sm);
               // d_vecds_lat(SVEC) = vel_grad_sm%d_vecds(SVEC)
               //
               //// CALL rot_mat_vecd(A, qr5x5_A)
               //
               // CALL rot_mat_wveccp(C_matx, qr3x3_ls) // amounts to qr3x3_ls = C_matx
               // CALL matt_x_vec_3(qr3x3_ls, vel_grad_sm%w_veccp, w_vec_lat)
               double w_vec_lat[ecmech::nwvec]; // assumes nwvec = ndim
               // vecsVMTa<ndim>(w_vec_lat, C_matx, _w_veccp_sm);

               get_xtal_frame_vel_grad_terms(d_vecd_lat, w_vec_lat, _d_vecd_sm,  _w_veccp_sm, C_matx, qr5x5_ls);

               //////////////////////////////
               // CALCULATIONS

               double T_vecds[ecmech::nsvec];
               _lattice_strain_prob.elas_strain_to_kirchoff_stress(T_vecds, e_vecd_f);

               double dgdot_dtau[SlipGeom::nslip] = { 0.0 }; // crys%tmp2_slp
               double dgdot_dg[SlipGeom::nslip] = { 0.0 }; // crys%tmp3_slp
               double pl_vecd[ecmech::ntvec] = { 0.0 };
               double pl_wvec[ecmech::nwvec] = { 0.0 }; // \pcDhat

               get_slip_rate_terms(dgdot_dtau, pl_vecd, pl_wvec, T_vecds, _kin_vals, _slipGeom, _kinetics);

               // Higher-order terms related to the elasticity stuff that's used in the residuals and jacobian calculation
               double A_e_M35[ecmech::nwvec * ecmech::ntvec];
               double ee_wvec[ecmech::nwvec];
               double ee_fac;
               elasticity_higher_order_terms(A_e_M35, ee_wvec, ee_fac, _a_V_ri, e_vecd_f, edot_vecd);

               // Residual Calculations
               _lattice_strain_prob.get_elas_strain_residual(resid, _epsdot_scale_inv, edot_vecd, pl_vecd, d_vecd_lat);
               _lattice_rot_prob.get_omega_residual(resid, _rotincr_scale_inv, ee_fac, xi_f, w_vec_lat, pl_wvec, ee_wvec);


               //////////////////////////////////////////////////////////////////////
               // JACOBIAN, fixed hardness and temperature
               //
               if (doComputeJ) {
                  // use RAJA::View machinery to simplify indexing for blocks in the Jacobian matrix ;
                  // can always swap this out later if it ends up being too heavyweight ;
                  // RAJA defaults to "row-major" -- final dimension indexing the fastest
                  //
                  // preliminaries
                  //
                  double dpl_deps_symm[ ecmech::ntvec * ecmech::ntvec ] = { 0.0 };
                  double dpl_deps_skew[ ecmech::nwvec * ecmech::ntvec ] = { 0.0 };
                  get_slip_rate_deriv_terms(dpl_deps_symm, dpl_deps_skew, dgdot_dtau, _a_V_ri, _slipGeom, _thermoElastN);

                  //
                  //
                  // derivatives with respect to lattice orientation changes
                  double dC_quat_dxi_T[ ecmech::nwvec * ecmech::qdim ];
                  double dDsm_dxi[ ecmech::ntvec * ecmech::nwvec ];
                  double dWsm_dxi[ ecmech::nwvec * ecmech::nwvec ];
                  eval_d_dxi_impl_quat(dC_quat_dxi_T, dDsm_dxi, dWsm_dxi,
                                       _d_vecd_sm, _w_veccp_sm,
                                       xi_f, _Cn_quat, C_matx, C_quat);

                  // d(B_S)/d(e_vecd_f)
                  //
                  _lattice_strain_prob.template get_deriv_elast_strain_wrt_elast_strain<nDimSys>(Jacobian, dpl_deps_symm);
                  // d(B_S)/d(xi_f)
                  //
                  // jacob_er = -dDsm_dxi(:,:)
                  _lattice_rot_prob.template get_deriv_elast_strain_wrt_omega<nDimSys>(Jacobian, dDsm_dxi);
                  // d(B_xi)/d(e_vecds_f)
                  //
                  _lattice_strain_prob.template get_deriv_omega_wrt_elast_strain<nDimSys, _i_sub_r>(Jacobian, edot_vecd, ee_fac, dpl_deps_skew, A_e_M35); 
                  // d(B_xi)/d(xi_f)
                  //
                  _lattice_rot_prob.template get_deriv_omega_wrt_omega<nDimSys>(Jacobian, dWsm_dxi);

                  if (_mtan_sI) {
                     double cauchy_stress_lattice[ ecmech::nsvec ];
                     _thermoElastN.getCauchy(cauchy_stress_lattice, T_vecds, _detV_ri);
                     get_material_tangent_stiffness<ThermoElastN, nDimSys, _i_sub_r>(_mtan_sI, Jacobian,
                                                                                     dC_quat_dxi_T, qr5x5_ls,
                                                                                     C_quat, C_matx,
                                                                                     cauchy_stress_lattice, _detV_ri,
                                                                                     _a_V_ri, _thermoElastN);
                  }

                  // SCALING
                  {
                     double scaleFactorJ;
                     for (int iJ = 0; iJ<_i_sub_r; ++iJ) {
                        // Jacobian(i_sub_e:i_sup_e,i_sub_e:i_sup_e) = jacob_ee * epsdot_scale_inv  * e_scale ! resid, x
                        scaleFactorJ = _epsdot_scale_inv * ecmech::e_scale;
                        for (int jJ = 0; jJ<_i_sub_r; ++jJ) { // <=_i_sup_e
                           int ijJ = ECMECH_NN_INDX(iJ, jJ, nDimSys);
                           Jacobian[ ijJ ] *= scaleFactorJ;
                        }

                        // Jacobian(i_sub_e:i_sup_e,i_sub_r:i_sup_r) = jacob_er * epsdot_scale_inv  * r_scale
                        scaleFactorJ = _epsdot_scale_inv * ecmech::r_scale;
                        for (int jJ = _i_sub_r; jJ<nDimSys; ++jJ) { // <_i_sup_r
                           int ijJ = ECMECH_NN_INDX(iJ, jJ, nDimSys);
                           Jacobian[ ijJ ] *= scaleFactorJ;
                        }
                     }

                     for (int iJ = _i_sub_r; iJ<nDimSys; ++iJ) {
                        // Jacobian(i_sub_r:i_sup_r,i_sub_e:i_sup_e) = jacob_re * rotincr_scale_inv * e_scale
                        scaleFactorJ = _rotincr_scale_inv * ecmech::e_scale;
                        for (int jJ = 0; jJ<_i_sub_r; ++jJ) { // <=_i_sup_e
                           int ijJ = ECMECH_NN_INDX(iJ, jJ, nDimSys);
                           Jacobian[ ijJ ] *= scaleFactorJ;
                        }

                        // Jacobian(i_sub_r:i_sup_r,i_sub_r:i_sup_r) = jacob_rr * rotincr_scale_inv * r_scale
                        scaleFactorJ = _rotincr_scale_inv * ecmech::r_scale;
                        for (int jJ = _i_sub_r; jJ<nDimSys; ++jJ) { // <_i_sup_r
                           int ijJ = ECMECH_NN_INDX(iJ, jJ, nDimSys);
                           Jacobian[ ijJ ] *= scaleFactorJ;
                        }
                     }
                  } // SCALING
               }
               return true;
            } // computeRJ

            __ecmech_hdev__
            inline
            const double* getGdot() const { return _gdot; };

            __ecmech_hdev__
            inline
            void get_slip_contribution(double& pl_disipation_rate,
                                       double& effective_shear_rate,
                                       double* const gdot,
                                       const double* const elas_strain
                                      )
            {
               get_slip_contributions(pl_disipation_rate, effective_shear_rate, gdot,
                                      _detV_ri, elas_strain, _kin_vals,
                                      _slipGeom, _kinetics, _lattice_strain_prob);
            }
                              

         private:

            const SlipGeom &_slipGeom;
            const Kinetics &_kinetics;
            const ThermoElastN &_thermoElastN;
            const EvptnLatticeStrainProblem<ThermoElastN> _lattice_strain_prob;
            const EvptnLatticeRotationProblem<ecmech::ntvec> _lattice_rot_prob;

            double _dt, _detV, _eVref, _p_EOS, _tK, _a_V;
            double _dt_ri, _a_V_ri, _detV_ri;

            double _hdn_scale;
            double _epsdot_scale_inv, _rotincr_scale_inv;

            double _gdot[SlipGeom::nslip]; // crys%tmp1_slp

            double _kin_vals[Kinetics::nVals];

            const double* const _h_state;
            const double* const _e_vecd_n;
            const double* const _Cn_quat;
            const double* const _d_vecd_sm; // d_vecds_sm would be fine too -- but do not use _d_vecd_sm[iSvecS];
            const double* const _w_veccp_sm;

            static constexpr size_t _nXnDim = nDimSys * nDimSys;
            static constexpr size_t _i_sub_e = 0; // ntvec
            static constexpr size_t _i_sub_r = ecmech::ntvec; // nwvec

            // for mtan (material tangent stiffnes)
            double* _mtan_sI; // null if not wanting tangent evaluation
      }; // class EvptnUpdstProblem

      /*
       * for steady-flow capability, might want to check out Dlsmm_getEnabled() stuff in EvpC.c
       *
       * convention for spin coming in should be consistent with w_veccp_sm convention
       */
      template<class SlipGeom, class Kinetics, class ThermoElastN, class EosModel>
      __ecmech_hdev__
      inline
      bool getResponseSngl(const SlipGeom& slipGeom,
                           const Kinetics& kinetics,
                           const ThermoElastN& elastN,
                           const EosModel& eos,
                           double    dt,
                           double    tolerance,
                           const double  * d_svec_kk_sm, // defRate,
                           const double  * w_veccp_sm, // spin
                           const double  * volRatio,
                           double  * eInt,
                           double  * stressSvecP,
                           double  * hist,
                           double  & tkelv,
                           double  * sdd,
                           double  * mtanSD,

                           int outputLevel = 0)
      {
         static const int iHistLbGdot = NumHist<SlipGeom, Kinetics, ThermoElastN, EosModel>::iHistLbGdot;

         // NOTE : mtanSD can be nullptr
         //
         const bool haveMtan = (mtanSD != nullptr);

         // convert deformation rate convention
         //
         double d_vecd_sm[ecmech::ntvec];
         svecToVecd(d_vecd_sm, d_svec_kk_sm);

         // pointers to state
         //
         double* h_state = &(hist[iHistLbH]);
         double* gdot = &(hist[iHistLbGdot]);
         //
         // copies, to keep beginning-of-step state safe
         //
         double e_vecd_n[ecmech::ntvec];

         for (int i_hist = 0; i_hist < ecmech::ntvec; i_hist++) {
            e_vecd_n[i_hist] = hist[iHistLbE + i_hist];
         }

         double quat_n[ecmech::qdim];
         for (int i_hist = 0; i_hist < ecmech::qdim; i_hist++) {
            quat_n[i_hist] = hist[iHistLbQ + i_hist];
         }

         //
         // normalize quat just in case
         vecsVNormalize<qdim>(quat_n);

         // total increment in the deviatoric part of the strain energy using
         // trapezoidal rule integration
         //
         // just beginning-of-step stress part so far
         //
         double halfVMidDt = oneqrtr * (volRatio[0] + volRatio[1]) * dt;
         double eDevTot = halfVMidDt * vecsInnerSvecDev(stressSvecP, d_svec_kk_sm);

         // EOS
         //
         double eOld = eInt[ecmech::i_ne_total];
         double pOld = stressSvecP[6];
         double pEOS, eNew, bulkNew;
         //
         // get tkelv from beginning-of-step to avoid tangent stiffness contributions
         {
            double pBOS;
            double vOld = volRatio[0];
            eos.evalPT(pBOS, tkelv, vOld, eOld);
         }
         //
         double tkelvNew;
         {
            double dpde, dpdv, dtde;
            updateSimple<EosModel>(eos, pEOS, tkelvNew, eNew, bulkNew,
                                   dpde, dpdv, dtde,
                                   volRatio[1], volRatio[3],
                                   eOld, pOld);
         }

         // update hardness state to the end of the step
         // gdot is still at beginning-of-step
         //
         double h_state_u[Kinetics::nH];
         double hvals[SlipGeom::nslip] = { 0.0 }; // additional values needed to update the hardening state
         if (SlipGeom::dynamic) {
            // For dynamic slip systems we need the chi angle
            double P[ecmech::ntvec * SlipGeom::nslip];
            double Q[ecmech::nwvec * SlipGeom::nslip];
            slipGeom.getPQ(hvals, P, Q, stressSvecP);
         }
         kinetics.updateH(h_state_u, h_state, dt, gdot, hvals, tkelv);

         double Cstr_vecds_lat[ecmech::nsvec];
         //
         double* e_vecd_u = &(hist[iHistLbE]);
         double* quat_u = &(hist[iHistLbQ]);
         double vNew = volRatio[1];
         {

            EvptnUpdstProblem prob(slipGeom, kinetics, elastN,
                                   dt,
                                   vNew, eNew, pEOS, tkelv,
                                   h_state_u, e_vecd_n, quat_n,
                                   d_vecd_sm, w_veccp_sm);

            snls::SNLSTrDlDenseG<decltype(prob)> solver(prob);

            snls::TrDeltaControl deltaControl;
            deltaControl._deltaInit = 1e0;
            {
               static const int maxIter = 200;
               solver.setupSolver(maxIter, tolerance, &deltaControl, outputLevel);
            }

            // set initial guess
            //
            for (int iX = 0; iX < prob.nDimSys; ++iX) {
               solver._x[iX] = 0e0;
            }

            snls::SNLSStatus_t status = solver.solve( );
            snls::SNLSStatus_t status2 = status;
            if (status != snls::converged) {
               ECMECH_WARN(__func__, "Trust Region Dogleg Solver failed to converge -- trying again with a Hybrid Nonlinear Solver");
               snls::SNLSHybrdTrDLDenseG<EvptnUpdstProblem<SlipGeom, Kinetics, ThermoElastN> > solver2(prob);
               static const int maxIter = 200;
               deltaControl._xiDecDelta = 0.6;
               solver2.setupSolver(maxIter, tolerance, &deltaControl, outputLevel);
               for (int iX = 0; iX < prob.nDimSys; ++iX) {
                  solver2.m_x[iX] = (status2 == snls::converged) ? solver2.m_x[iX] : 0.0;
               }
               status2 = solver2.solve( );
               for (int iX = 0; iX < prob.nDimSys; ++iX) {
                  solver._x[iX] = (status2 == snls::converged) ? solver2.m_x[iX] : 0.0;
               }
#ifdef __ecmech_host_only__
               if (status2 != snls::converged) {
                  std::cout << "trust region solver residual " << solver.getRes() << " exit status " << status << std::endl;
                  std::cout << "hybrid solver residual " << solver2.getRes() << " exit status " << status2 << std::endl;
               }
#endif
            }
            //
            if (status != snls::converged && status2 != snls::converged) {
#ifdef __ecmech_host_only__
               ECMECH_WARN(__func__, "Both solvers failed to converge -- will try again with implicit elastic strain solve only");
#endif
               return false;
               // ECMECH_FAIL(__func__, "Solver failed to converge!");
            }

            if (haveMtan) {
               double mtanSD_vecds[ ecmech::nsvec2 ];
               prob.provideMTan(mtanSD_vecds);

               {
                  double residual[prob.nDimSys], Jacobian[prob.nDimSys * prob.nDimSys];
                  solver.computeRJ(&residual[0], &Jacobian[0]);
               }

               prob.clearMTan();

               // currently have derivative with-respsect-to deformation rate ;
               // to get derivative with-respsect-to strain increment,
               // multiply by 1/dt
               //
               double dt_ri = prob.getDtRi();
               for (int i = 0; i<ecmech::nsvec2; ++i) {
                  mtanSD_vecds[i] = mtanSD_vecds[i] * dt_ri;
               }

               // contribution to stiffness from EOS
               // this is a bit crude, but should do the trick for now;
               // neglects effect of pEOS and vNew on workings of evptn
               //
               mtanSD_vecds[ECMECH_NN_INDX(iSvecS, iSvecS, ecmech::nsvec)] = three * bulkNew;

               // convert from vecds notation to svec notation
               //
               mtan_conv_sd_svec<true>(mtanSD, mtanSD_vecds);
            }

            // store updated state
            //
            prob.stateFromX(e_vecd_u, quat_u, solver._x);
            for (int i_hstate = 0; i_hstate < Kinetics::nH; i_hstate++) {
               h_state[i_hstate] = h_state_u[i_hstate];
            }

            double pl_disipation_rate = 0.0;
            double effective_shear_rate = 0.0;

            prob.get_slip_contribution(pl_disipation_rate, effective_shear_rate,
                                        gdot, e_vecd_u);

            //
            // {
            //    const double* gdot_u = prob.getGdot();
            //    for (int i_gdot = 0; i_gdot < SlipGeom::nslip; i_gdot++) {
            //       gdot[i_gdot] = gdot_u[i_gdot];
            //    }
            // }
            //
            hist[iHistA_shrateEff] = effective_shear_rate;
            hist[iHistA_shrEff] += hist[iHistA_shrateEff] * dt;
            //
            {
               double dEff = vecd_Deff(d_vecd_sm);
               double flow_strength = prob.getHdnScale();
               if (dEff > idp_tiny_sqrt) {
                  flow_strength = pl_disipation_rate / dEff;
               }
               hist[iHistA_flowStr] = flow_strength;
            }
            //
            hist[iHistA_nFEval] = solver.getNFEvals(); // does _not_ include updateH iterations

            // get Cauchy stress
            //
            prob.elastNEtoC(Cstr_vecds_lat, e_vecd_u);
         }

         double C_matx[ecmech::ndim * ecmech::ndim];
         quat_to_tensor(C_matx, quat_u);
         //
         double qr5x5_ls[ecmech::ntvec * ecmech::ntvec];
         get_rot_mat_vecd(qr5x5_ls, C_matx);
         //
         double Cstr_vecds_sm[ecmech::nsvec];
         vecsVMa<ntvec>(Cstr_vecds_sm, qr5x5_ls, Cstr_vecds_lat);
         Cstr_vecds_sm[iSvecS] = Cstr_vecds_lat[iSvecS];
         //
         // put end-of-step stress in stressSvecP
         vecdsToSvecP(stressSvecP, Cstr_vecds_sm);
         //
         // and now the second half of the trapezoidal integration
         //
         eDevTot += halfVMidDt * vecsInnerSvecDev(stressSvecP, d_svec_kk_sm);

         // adjust sign on quat so that as close as possible to quat_o;
         // more likely to keep orientations clustered this way;
         // this flip through the origin is equivalent under antipodal symmetry
         //
         if (vecsyadotb<qdim>(quat_u, quat_n) < zero) {
            for (int iQ = 0; iQ<ecmech::qdim; ++iQ) {
               quat_u[iQ] = -quat_u[iQ];
            }
         }

         {
            double gmod = elastN.getGmod(tkelv, pEOS, eNew);
            sdd[i_sdd_bulk] = bulkNew;
            sdd[i_sdd_gmod] = gmod;
         }
#ifdef ECMECH_DEBUG
         assert(ecmech::nsdd == 2);
#endif

         eNew = eNew + eDevTot;
         //
         // could update pressure and temperature again, but do not bother

         eInt[ecmech::i_ne_total] = eNew;
#ifdef ECMECH_DEBUG
         assert(ecmech::ne == 1);
#endif
         return true;
      } // getResponseSngl
   } // namespace evptn
} // namespace ecmech

#endif // ECMECH_EVPTN_H
