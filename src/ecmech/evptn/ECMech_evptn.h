// -*-c++-*-

#ifndef ECMECH_EVPTN_H
#define ECMECH_EVPTN_H

#include <cassert>

#include "ECMech_core.h"
#include "ECMech_util.h"
#include "evptn/ECMech_base_classes.h"
#include "evptn/ECMech_base_fcns.h"

namespace ecmech {
   namespace evptn {

      template<class SlipGeom, class Kinetics, class ThermoElastN, class ProblemState>
      class EvptnUpdstProblem
      {
         public:

            static constexpr int nDimSys = ecmech::ntvec + ecmech::nwvec;

            // constructor
            __ecmech_hdev__
            EvptnUpdstProblem(const SlipGeom& slipGeom,
                              const Kinetics& kinetics,
                              const ThermoElastN& thermoElastN,
                              ProblemState& prob_state
                              )
               : m_slipGeom(slipGeom),
               m_kinetics(kinetics),
               m_lattice_strain_prob(thermoElastN, prob_state.dt, prob_state.rel_vol_new, prob_state.energy_new, prob_state.pressure_EOS, prob_state.tkelv, prob_state.elast_d5_n),
               m_lattice_rot_prob(prob_state.dt, prob_state.quat_n),
               m_def_rate_d5_sample(prob_state.def_rate_d5_sample), // vel_grad_sm%d_vecds
               m_spin_vec_sample(prob_state.spin_vec_sample), // vel_grad_sm%w_veccp
               m_mtan_sI(nullptr)
            {
               m_hdn_scale = m_kinetics.getVals(m_kin_vals, prob_state.pressure_EOS, prob_state.tkelv, prob_state.h_state_u);

               double adots_ref = m_kinetics.getFixedRefRate(m_kin_vals);
               double eff = vecNorm<ntvec>(m_def_rate_d5_sample); // do not worry about factor of sqrt(twothird)
               if (eff < epsdot_scl_nzeff * adots_ref) {
                  m_epsdot_scale_inv = one / adots_ref;
               }
               else {
                  m_epsdot_scale_inv = fmin(one / eff, 1e6 * m_lattice_strain_prob.m_dt);
               }
               //
               m_rotincr_scale_inv = m_lattice_strain_prob.m_inv_dt * m_epsdot_scale_inv;
            }

            // deconstructor
            __ecmech_hdev__
            ~EvptnUpdstProblem() {}

            __ecmech_hdev__
            inline
            void provideMTan(double* mtan_sI) { m_mtan_sI = mtan_sI; }

            __ecmech_hdev__
            inline
            void clearMTan( ) { m_mtan_sI = nullptr; }

            __ecmech_hdev__
            inline
            double getDtRi() const { return m_lattice_strain_prob.m_inv_dt; }

            __ecmech_hdev__
            inline
            double getHdnScale() const { return m_hdn_scale; }

            /*
             * NOTES :
             * () should be equivalent to what happens in computeRJ
             * () not necessarily safe if elast_dev_press_vec is the same memory as m_elast_d5_n or quat is the same as _xtal_ori_quat_n
             */
            __ecmech_hdev__
            inline
            void stateFromX(double* const elast_dev_press_vec,
                            double* const quat,
                            const double* const x) {
               m_lattice_strain_prob.stateFromX(elast_dev_press_vec,  &(x[m_i_sub_e]));
               m_lattice_rot_prob.stateFromX(quat,  &(x[m_i_sub_r]));
            }

            __ecmech_hdev__
            inline
            void elastNEtoC(double* const cauchy_xtal, // nsvec
                            const double* const elast_d5_f // ntvec
                            ) const {
               m_lattice_strain_prob.elast_strain_to_cauchy_stress(cauchy_xtal, elast_d5_f);
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
                  for (size_t ijJ = 0; ijJ<m_nXnDim; ++ijJ) {
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
               double elast_dt_d5[ecmech::ntvec];
               vecsVxa<ntvec>(elast_dt_d5, ecmech::e_scale, &(x[m_i_sub_e]) ); // elast_dt_d5 is now the delta, _not_ yet elast_dt_d5
               // elast_d5_f is end-of-step
               double elast_d5_f[ntvec];
               vecsVapb<ntvec>(elast_d5_f, elast_dt_d5, m_lattice_strain_prob.m_elast_d5_n);
               vecsVsa<ntvec>(elast_dt_d5, m_lattice_strain_prob.m_inv_dt); // _now_ elast_dt_d5 has elast_dt_d5
               //
               double xi_f[nwvec];
               vecsVxa<nwvec>(xi_f, ecmech::r_scale, &(x[m_i_sub_r]) );

               double A_quat[ecmech::qdim];
               emap_to_quat(A_quat, xi_f);
               //
               double xtal_ori_quat[ecmech::qdim];
               get_c_quat(xtal_ori_quat, A_quat, m_lattice_rot_prob.m_xtal_ori_quat_n);
               //
               double xtal_rmat[ecmech::ndim * ecmech::ndim];
               quat_to_tensor(xtal_rmat, xtal_ori_quat);
               //
               double rmat_5x5_sample2xtal[ecmech::ntvec * ecmech::ntvec];
               get_rot_mat_vecd(rmat_5x5_sample2xtal, xtal_rmat);

               double def_rate_d5_xtal[ecmech::ntvec];
               double spin_vec_xtal[ecmech::nwvec]; // assumes nwvec = ndim

               get_xtal_frame_vel_grad_terms(def_rate_d5_xtal, spin_vec_xtal, m_def_rate_d5_sample,  m_spin_vec_sample, xtal_rmat, rmat_5x5_sample2xtal);

               //////////////////////////////
               // CALCULATIONS

               double kirchoff[ecmech::nsvec];
               m_lattice_strain_prob.elast_strain_to_kirchoff_stress(kirchoff, elast_d5_f);

               double dgdot_dtau[SlipGeom::nslip] = { 0.0 }; // crys%tmp2_slp
               double plastic_def_rate_d5[ecmech::ntvec] = { 0.0 };
               double plastic_spin_vec[ecmech::nwvec] = { 0.0 }; // \pcDhat

               get_slip_rate_terms(dgdot_dtau, plastic_def_rate_d5, plastic_spin_vec, kirchoff, m_kin_vals, m_slipGeom, m_kinetics);

               // Higher-order terms related to the elasticity stuff that's used in the residuals and jacobian calculation
               double A_e_M35[ecmech::nwvec * ecmech::ntvec];
               double ee_wvec[ecmech::nwvec];
               double ee_fac;
               elasticity_higher_order_terms(A_e_M35, ee_wvec, ee_fac, m_lattice_strain_prob.m_inv_a_vol, elast_d5_f, elast_dt_d5);

               // Residual Calculations
               m_lattice_strain_prob.get_elast_strain_residual(resid, m_epsdot_scale_inv, elast_dt_d5, plastic_def_rate_d5, def_rate_d5_xtal);
               m_lattice_rot_prob.get_omega_residual(resid, m_rotincr_scale_inv, ee_fac, xi_f, spin_vec_xtal, plastic_spin_vec, ee_wvec);

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
                  get_slip_rate_deriv_terms(dpl_deps_symm, dpl_deps_skew, dgdot_dtau, m_lattice_strain_prob.m_inv_a_vol, m_slipGeom, m_lattice_strain_prob.m_thermo_elast_n);

                  // derivatives with respect to lattice orientation changes
                  double dxtal_ori_quat_dxi_T[ ecmech::nwvec * ecmech::qdim ];
                  double dDsm_dxi[ ecmech::ntvec * ecmech::nwvec ];
                  double dWsm_dxi[ ecmech::nwvec * ecmech::nwvec ];
                  eval_d_dxi_impl_quat(dxtal_ori_quat_dxi_T, dDsm_dxi, dWsm_dxi,
                                       m_def_rate_d5_sample, m_spin_vec_sample,
                                       xi_f, 
                                       m_lattice_rot_prob.m_xtal_ori_quat_n,
                                       xtal_rmat, xtal_ori_quat);

                  // d(B_S)/d(elast_d5_f)
                  //
                  m_lattice_strain_prob.template get_deriv_elast_strain_wrt_elast_strain<nDimSys>(Jacobian, dpl_deps_symm);
                  // d(B_S)/d(xi_f)
                  //
                  // jacob_er = -dDsm_dxi(:,:)
                  m_lattice_rot_prob.template get_deriv_elast_strain_wrt_omega<nDimSys>(Jacobian, dDsm_dxi);
                  // d(B_xi)/d(elast_dev_press_vecs_f)
                  //
                  m_lattice_strain_prob.template get_deriv_omega_wrt_elast_strain<nDimSys, m_i_sub_r>(Jacobian, elast_dt_d5, ee_fac, dpl_deps_skew, A_e_M35);
                  // d(B_xi)/d(xi_f)
                  //
                  m_lattice_rot_prob.template get_deriv_omega_wrt_omega<nDimSys>(Jacobian, dWsm_dxi);

                  if (m_mtan_sI) {
                     double cauchy_stress_lattice[ ecmech::nsvec ];
                     m_lattice_strain_prob.m_thermo_elast_n.getCauchy(cauchy_stress_lattice, kirchoff, m_lattice_strain_prob.m_inv_det_v_e);
                     get_material_tangent_stiffness<ThermoElastN, nDimSys, m_i_sub_r>
                     (m_mtan_sI, Jacobian,
                     dxtal_ori_quat_dxi_T, rmat_5x5_sample2xtal,
                     xtal_ori_quat, xtal_rmat,
                     cauchy_stress_lattice,
                     m_lattice_strain_prob.m_inv_det_v_e,
                     m_lattice_strain_prob.m_inv_a_vol,
                     m_lattice_strain_prob.m_thermo_elast_n);
                  }

                  // SCALING
                  {
                     double scaleFactorJ;
                     for (size_t iJ = 0; iJ<m_i_sub_r; ++iJ) {
                        // Jacobian(i_sub_e:i_sup_e,i_sub_e:i_sup_e) = jacob_ee * epsdot_scale_inv  * e_scale ! resid, x
                        scaleFactorJ = m_epsdot_scale_inv * ecmech::e_scale;
                        for (size_t jJ = 0; jJ<m_i_sub_r; ++jJ) { // <=_i_sup_e
                           int ijJ = ECMECH_NN_INDX(iJ, jJ, nDimSys);
                           Jacobian[ ijJ ] *= scaleFactorJ;
                        }

                        // Jacobian(i_sub_e:i_sup_e,i_sub_r:i_sup_r) = jacob_er * epsdot_scale_inv  * r_scale
                        scaleFactorJ = m_epsdot_scale_inv * ecmech::r_scale;
                        for (int jJ = m_i_sub_r; jJ<nDimSys; ++jJ) { // <_i_sup_r
                           int ijJ = ECMECH_NN_INDX(iJ, jJ, nDimSys);
                           Jacobian[ ijJ ] *= scaleFactorJ;
                        }
                     }

                     for (int iJ = m_i_sub_r; iJ<nDimSys; ++iJ) {
                        // Jacobian(i_sub_r:i_sup_r,i_sub_e:i_sup_e) = jacob_re * rotincr_scale_inv * e_scale
                        scaleFactorJ = m_rotincr_scale_inv * ecmech::e_scale;
                        for (size_t jJ = 0; jJ<m_i_sub_r; ++jJ) { // <=_i_sup_e
                           int ijJ = ECMECH_NN_INDX(iJ, jJ, nDimSys);
                           Jacobian[ ijJ ] *= scaleFactorJ;
                        }

                        // Jacobian(i_sub_r:i_sup_r,i_sub_r:i_sup_r) = jacob_rr * rotincr_scale_inv * r_scale
                        scaleFactorJ = m_rotincr_scale_inv * ecmech::r_scale;
                        for (size_t jJ = m_i_sub_r; jJ<nDimSys; ++jJ) { // <_i_sup_r
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
            void get_slip_contribution(double& pl_disipation_rate,
                                       double& effective_shear_rate,
                                       double* const gdot,
                                       const double* const elast_strain
                                      )
            {
               get_slip_contributions(pl_disipation_rate, effective_shear_rate, gdot,
                                      m_lattice_strain_prob.m_inv_det_v_e, elast_strain, m_kin_vals,
                                      m_slipGeom, m_kinetics, m_lattice_strain_prob);
            }
                              

         private:

            const SlipGeom &m_slipGeom;
            const Kinetics &m_kinetics;
            const EvptnLatticeStrainProblem<ThermoElastN> m_lattice_strain_prob;
            const EvptnLatticeRotationProblem<ecmech::ntvec> m_lattice_rot_prob;

            double m_hdn_scale;
            double m_epsdot_scale_inv, m_rotincr_scale_inv;

            double m_kin_vals[Kinetics::nVals];

            const double* const m_def_rate_d5_sample; // d_vecds_sm would be fine too -- but do not use m_def_rate_d5_sample[iSvecS];
            const double* const m_spin_vec_sample;

            static constexpr size_t m_nXnDim = nDimSys * nDimSys;
            static constexpr size_t m_i_sub_e = 0; // ntvec
            static constexpr size_t m_i_sub_r = ecmech::ntvec; // nwvec

            // for mtan (material tangent stiffnes)
            double* m_mtan_sI; // null if not wanting tangent evaluation
      }; // class EvptnUpdstProblem

#if defined(ECMECH_EXTRA_SOLVERS)

      template<class SlipGeom, class ThermoElastN, class ProblemState>
      class RotUpdProblem
      {
         public:
         static constexpr int nDimSys = ecmech::nwvec;

         __ecmech_hdev__
         RotUpdProblem(const SlipGeom& slipGeom,
                        const ThermoElastN& thermoElastN,
                        ProblemState& prob_state
                        ) :
            m_slipGeom(slipGeom),
            m_lattice_strain_prob(thermoElastN, prob_state.dt, prob_state.rel_vol_new, prob_state.energy_new, prob_state.pressure_EOS, prob_state.tkelv, prob_state.elast_d5_n),
            m_lattice_rot_prob(prob_state.dt, prob_state.quat_n),
            m_def_rate_d5_sample(prob_state.def_rate_d5_sample), // vel_grad_sm%d_vecds
            m_spin_vec_sample(prob_state.spin_vec_sample) // vel_grad_sm%w_veccp
         {

            double adots_ref = vecNorm<SlipGeom::nslip>(prob_state.gdot);

            double eff = vecNorm<ecmech::ntvec>(m_def_rate_d5_sample); // do not worry about factor of sqrt(twothird)
            if (eff < epsdot_scl_nzeff * adots_ref) {
                  m_epsdot_scale_inv = one / adots_ref;
            }
            else {
                  m_epsdot_scale_inv = fmin(one / eff, 1e6 * m_lattice_strain_prob.m_dt);
            }
            //
            m_rotincr_scale_inv = m_lattice_strain_prob.m_inv_dt * m_epsdot_scale_inv;

            vecsVMa<ntvec, SlipGeom::nslip>(m_plastic_def_rate_d5, slipGeom.getP(), prob_state.gdot);
            vecsVMa<nwvec, SlipGeom::nslip>(m_plastic_spin_vec, slipGeom.getQ(), prob_state.gdot);

         }

         // deconstructor
         __ecmech_hdev__
         ~RotUpdProblem() {}

         __ecmech_hdev__
         inline
         void stateFromX(double* const quat,
                        const double* const x) {
            m_lattice_rot_prob.stateFromX(quat, x);
         }

         __ecmech_hdev__
         bool computeRJ(double* const resid,
                        double* const Jacobian,
                        const double* const x) {
            bool doComputeJ = (Jacobian != nullptr);

            if (doComputeJ) {
                  // zero the Jacobian so that do not need to worry about zero
                  // entries in the midst of other things later
                  //
                  for (int ijJ = 0; ijJ< m_nXnDim; ++ijJ) {
                     Jacobian[ijJ] = 0.0;
                  }
            }
            //
            for (int iR = 0; iR<nDimSys; ++iR) {
                  resid[iR] = 0.0;
            }

            double xi_f[nwvec];
            vecsVxa<nwvec>(xi_f, ecmech::r_scale, x);

            double A_quat[ecmech::qdim];
            emap_to_quat(A_quat, xi_f);
            //
            double xtal_ori_quat[ecmech::qdim];
            get_c_quat(xtal_ori_quat, A_quat, m_lattice_rot_prob.m_xtal_ori_quat_n);
            //
            double xtal_rmat[ecmech::ndim * ecmech::ndim];
            quat_to_tensor(xtal_rmat, xtal_ori_quat);
            //
            double rmat_5x5_sample2xtal[ecmech::ntvec * ecmech::ntvec];
            get_rot_mat_vecd(rmat_5x5_sample2xtal, xtal_rmat);

            double def_rate_d5_xtal[ecmech::ntvec];
            double spin_vec_xtal[ecmech::nwvec]; // assumes nwvec = ndim
            get_xtal_frame_vel_grad_terms(def_rate_d5_xtal, spin_vec_xtal, m_def_rate_d5_sample,  m_spin_vec_sample, xtal_rmat, rmat_5x5_sample2xtal);

            double elast_dt_d5[ecmech::ntvec];
            // Calculate what this elast_dt_d5 term should be given the current
            // state information.
            for (int i = 0; i < ecmech::ntvec; i++)
            {
                  elast_dt_d5[i] = m_lattice_strain_prob.m_inv_a_vol * (def_rate_d5_xtal[i] - m_plastic_def_rate_d5[i]);
            }

            // Higher-order terms related to the elasticity stuff that's used in the residuals and jacobian calculation
            double A_e_M35[ecmech::nwvec * ecmech::ntvec];
            double ee_wvec[ecmech::nwvec];
            double ee_fac;
            elasticity_higher_order_terms(A_e_M35, ee_wvec, ee_fac, m_lattice_strain_prob.m_inv_a_vol, m_lattice_strain_prob.m_elast_d5_n, elast_dt_d5);

            // Residual Calculations
            m_lattice_rot_prob.get_omega_residual(resid, m_rotincr_scale_inv, ee_fac, xi_f, spin_vec_xtal, m_plastic_spin_vec, ee_wvec);

            //////////////////////////////////////////////////////////////////////
            // JACOBIAN, fixed hardness and temperature
            //
            if (doComputeJ) {

                  //
                  //
                  // derivatives with respect to lattice orientation changes
                  double dxtal_ori_quat_dxi_T[ ecmech::nwvec * ecmech::qdim ];
                  double dDsm_dxi[ ecmech::ntvec * ecmech::nwvec ];
                  double dWsm_dxi[ ecmech::nwvec * ecmech::nwvec ];
                  eval_d_dxi_impl_quat(dxtal_ori_quat_dxi_T, dDsm_dxi, dWsm_dxi,
                                    m_def_rate_d5_sample, m_spin_vec_sample,
                                    xi_f,
                                    m_lattice_rot_prob.m_xtal_ori_quat_n,
                                    xtal_rmat, xtal_ori_quat);

                  // d(B_xi)/d(xi_f)
                  //
                  m_lattice_rot_prob.template get_deriv_omega_wrt_omega<nDimSys>(Jacobian, dWsm_dxi);

                  const double scaleFactorJ = m_rotincr_scale_inv * ecmech::r_scale;
                  for (int iJ = 0; iJ<nDimSys; ++iJ) {
                     // Jacobian(i_sub_r:i_sup_r,i_sub_r:i_sup_r) = jacob_rr * rotincr_scale_inv * r_scale
                     for (int jJ = 0; jJ<nDimSys; ++jJ) { // <_i_sup_r
                        int ijJ = ECMECH_NN_INDX(iJ, jJ, nDimSys);
                        Jacobian[ ijJ ] *= scaleFactorJ;
                     }
                  }
            }
            return true;
         }

         private:

         const SlipGeom &m_slipGeom;
         const EvptnLatticeStrainProblem<ThermoElastN> m_lattice_strain_prob;
         const EvptnLatticeRotationProblem<0> m_lattice_rot_prob;

         double m_epsdot_scale_inv, m_rotincr_scale_inv;

         double m_plastic_def_rate_d5[ecmech::ntvec];
         double m_plastic_spin_vec[ecmech::nwvec]; // \pcDhat

         const double* const m_def_rate_d5_sample; // d_vecds_sm would be fine too -- but do not use m_def_rate_d5_sample[iSvecS];
         const double* const m_spin_vec_sample;

         static constexpr int m_nXnDim = nDimSys * nDimSys;
      };

      template<class SlipGeom, class Kinetics, class ThermoElastN, class ProblemState>
      class EvptnNRUpdstProblem
      {
         public:

            static constexpr int nDimSys = ecmech::ntvec;

            // constructor
            __ecmech_hdev__
            EvptnNRUpdstProblem(const SlipGeom& slipGeom,
                                const Kinetics& kinetics,
                                const ThermoElastN& thermoElastN,
                                ProblemState& prob_state
                               )
                        :
                        m_slipGeom(slipGeom),
                        m_kinetics(kinetics),
                        m_thermoElastN(thermoElastN),
                        m_lattice_strain_prob(thermoElastN, prob_state.dt, prob_state.rel_vol_new, prob_state.energy_new, prob_state.pressure_EOS, prob_state.tkelv, prob_state.elast_d5_n),
                        m_lattice_rot_prob(prob_state.dt, prob_state.quat_n),
                        m_elast_d5_n(prob_state.elast_d5_n),
                        m_xtal_ori_quat(prob_state.quat_u),
                        m_def_rate_d5_sample(prob_state.def_rate_d5_sample), // vel_grad_sm%d_vecds
                        m_spin_vec_sample(prob_state.spin_vec_sample), // vel_grad_sm%w_veccp
                        m_mtan_sI(nullptr)
            {
               m_hdn_scale = m_kinetics.getVals(m_kin_vals, prob_state.pressure_EOS, prob_state.tkelv, prob_state.h_state_u);

               double adots_ref = m_kinetics.getFixedRefRate(m_kin_vals);
               double eff = vecNorm<ntvec>(m_def_rate_d5_sample); // do not worry about factor of sqrt(twothird)
               if (eff < (epsdot_scl_nzeff * adots_ref)) {
                  m_epsdot_scale_inv = one / adots_ref;
               }
               else {
                  m_epsdot_scale_inv = fmin(one / eff, 1e6 * m_lattice_strain_prob.m_dt);
               }
               //
               m_rotincr_scale_inv = m_lattice_strain_prob.m_inv_dt * m_epsdot_scale_inv;
            }

            // deconstructor
            __ecmech_hdev__
            ~EvptnNRUpdstProblem() {}

            __ecmech_hdev__
            inline
            void provideMTan(double* mtan_sI) { m_mtan_sI = mtan_sI; }

            __ecmech_hdev__
            inline
            void clearMTan( ) { m_mtan_sI = nullptr; }

            __ecmech_hdev__
            inline
            double getDtRi() const { return m_lattice_strain_prob.m_inv_dt; }

            __ecmech_hdev__
            inline
            double getHdnScale() const { return m_hdn_scale; }

            /*
             * NOTES :
             * () should be equivalent to what happens in computeRJ
             * () not necessarily safe if elast_dev_press_vec is the same memory as m_elast_d5_n or quat is the same as _xtal_ori_quat_n
             */
            __ecmech_hdev__
            inline
            void stateFromX(double* const elast_dev_press_vec,
                            const double* const x) {
               m_lattice_strain_prob.stateFromX(elast_dev_press_vec,  &(x[m_i_sub_e]));
            }

            __ecmech_hdev__
            inline
            void elastNEtoC(double* const cauchy_xtal, // nsvec
                            const double* const elast_d5_f // ntvec
                            ) const {
               m_lattice_strain_prob.elast_strain_to_cauchy_stress(cauchy_xtal, elast_d5_f);
            }

            __ecmech_hdev__
            bool computeRJ(double* const resid,
                           double* const Jacobian,
                           const double* const x) {

               bool doComputeJ = (Jacobian != nullptr);

               if (doComputeJ) {
                  // zero the Jacobian so that do not need to worry about zero
                  // entries in the midst of other things later
                  //
                  for (size_t ijJ = 0; ijJ<m_nXnDim; ++ijJ) {
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
               double elast_dt_d5[ecmech::ntvec];
               vecsVxa<ntvec>(elast_dt_d5, ecmech::e_scale, &(x[m_i_sub_e]) ); // elast_dt_d5 is now the delta, _not_ yet elast_dt_d5
               // elast_d5_f is end-of-step
               double elast_d5_f[ntvec];
               vecsVapb<ntvec>(elast_d5_f, elast_dt_d5, m_lattice_strain_prob.m_elast_d5_n);
               vecsVsa<ntvec>(elast_dt_d5, m_lattice_strain_prob.m_inv_dt); // _now_ elast_dt_d5 has elast_dt_d5
               //
               double xtal_rmat[ecmech::ndim * ecmech::ndim];
               quat_to_tensor(xtal_rmat, m_xtal_ori_quat);
               //
               double rmat_5x5_sample2xtal[ecmech::ntvec * ecmech::ntvec];
               get_rot_mat_vecd(rmat_5x5_sample2xtal, xtal_rmat);
               double def_rate_d5_xtal[ecmech::ntvec];
               double spin_vec_xtal[ecmech::nwvec]; // assumes nwvec = ndim

               get_xtal_frame_vel_grad_terms(def_rate_d5_xtal, spin_vec_xtal, m_def_rate_d5_sample,  m_spin_vec_sample, xtal_rmat, rmat_5x5_sample2xtal);

               //////////////////////////////
               // CALCULATIONS

               double kirchoff[ecmech::nsvec];
               m_lattice_strain_prob.elast_strain_to_kirchoff_stress(kirchoff, elast_d5_f);

               double dgdot_dtau[SlipGeom::nslip] = { 0.0 }; // crys%tmp2_slp
               double plastic_def_rate_d5[ecmech::ntvec] = { 0.0 };
               double plastic_spin_vec[ecmech::nwvec] = { 0.0 }; // \pcDhat

               get_slip_rate_terms(dgdot_dtau, plastic_def_rate_d5, plastic_spin_vec, kirchoff, m_kin_vals, m_slipGeom, m_kinetics);

               // Residual Calculations
               m_lattice_strain_prob.get_elast_strain_residual(resid, m_epsdot_scale_inv, elast_dt_d5, plastic_def_rate_d5, def_rate_d5_xtal);

               //////////////////////////////////////////////////////////////////////
               // JACOBIAN, fixed hardness and temperature
               //
               if (doComputeJ) {
                  //
                  // preliminaries
                  //
                  double dpl_deps_symm[ ecmech::ntvec * ecmech::ntvec ] = { 0.0 };
                  double dpl_deps_skew[ ecmech::nwvec * ecmech::ntvec ] = { 0.0 };
                  get_slip_rate_deriv_terms(dpl_deps_symm, dpl_deps_skew, dgdot_dtau, m_lattice_strain_prob.m_inv_a_vol, m_slipGeom, m_lattice_strain_prob.m_thermo_elast_n);

                  // d(B_S)/d(elast_d5_f)
                  //
                  m_lattice_strain_prob.template get_deriv_elast_strain_wrt_elast_strain<nDimSys>(Jacobian, dpl_deps_symm);

                  if (m_mtan_sI) {

                     // Higher-order terms related to the elasticity stuff that's used in the residuals and jacobian calculation
                     double A_e_M35[ecmech::nwvec * ecmech::ntvec];
                     double ee_wvec[ecmech::nwvec];
                     double ee_fac;
                     double xi_f[ecmech::nwvec] = {};

                     m_lattice_rot_prob.deltaOmegaFromState(xi_f, m_xtal_ori_quat, m_lattice_rot_prob.m_xtal_ori_quat_n);

                     elasticity_higher_order_terms(A_e_M35, ee_wvec, ee_fac, m_lattice_strain_prob.m_inv_a_vol, elast_d5_f, elast_dt_d5);

                     // derivatives with respect to lattice orientation changes
                     double dxtal_ori_quat_dxi_T[ ecmech::nwvec * ecmech::qdim ];
                     double dDsm_dxi[ ecmech::ntvec * ecmech::nwvec ];
                     double dWsm_dxi[ ecmech::nwvec * ecmech::nwvec ];
                     eval_d_dxi_impl_quat(dxtal_ori_quat_dxi_T, dDsm_dxi, dWsm_dxi,
                                          m_def_rate_d5_sample, m_spin_vec_sample,
                                          xi_f,
                                          m_lattice_rot_prob.m_xtal_ori_quat_n,
                                          xtal_rmat, m_xtal_ori_quat);

                     static constexpr int nDimSolve = nDimSys + ecmech::nwvec;
                     static constexpr int nDimSolve2 = nDimSolve * nDimSolve;

                     double Jacobian2[nDimSolve2] = {};

                     RAJA::View<double, RAJA::Layout<2> > pfrac_ee(Jacobian2, nDimSolve, nDimSolve);
                     RAJA::View<double, RAJA::Layout<2> > jacob_ee(Jacobian, nDimSys, nDimSys);
                     for (int i_jac = 0; i_jac < nDimSys; i_jac++) {
                        for (int j_jac = 0; j_jac < nDimSys; j_jac++)
                        pfrac_ee(i_jac, j_jac) = jacob_ee(i_jac, j_jac);
                     }

                     // d(B_S)/d(xi_f)
                     //
                     // jacob_er = -dDsm_dxi(:,:)
                     m_lattice_rot_prob.template get_deriv_elast_strain_wrt_omega<nDimSolve>(Jacobian2, dDsm_dxi);
                     // d(B_xi)/d(elast_dev_press_vecs_f)
                     //
                     m_lattice_strain_prob.template get_deriv_omega_wrt_elast_strain<nDimSolve, m_i_sub_r>(Jacobian2, elast_dt_d5, ee_fac, dpl_deps_skew, A_e_M35);
                     // d(B_xi)/d(xi_f)
                     //
                     m_lattice_rot_prob.template get_deriv_omega_wrt_omega<nDimSolve>(Jacobian2, dWsm_dxi);

                     double cauchy_stress_lattice[ ecmech::nsvec ];
                     m_lattice_strain_prob.m_thermo_elast_n.getCauchy(cauchy_stress_lattice, kirchoff, m_lattice_strain_prob.m_inv_det_v_e);
                     get_material_tangent_stiffness<ThermoElastN, nDimSolve, m_i_sub_r>(m_mtan_sI, Jacobian2,
                     dxtal_ori_quat_dxi_T, rmat_5x5_sample2xtal,
                     m_xtal_ori_quat, xtal_rmat,
                     cauchy_stress_lattice,
                     m_lattice_strain_prob.m_inv_det_v_e,
                     m_lattice_strain_prob.m_inv_a_vol,
                     m_lattice_strain_prob.m_thermo_elast_n);
                  }

                  // SCALING
                  {
                     double scaleFactorJ;
                     for (size_t iJ = 0; iJ<m_i_sub_r; ++iJ) {
                        // Jacobian(i_sub_e:i_sup_e,i_sub_e:i_sup_e) = jacob_ee * epsdot_scale_inv  * e_scale ! resid, x
                        scaleFactorJ = m_epsdot_scale_inv * ecmech::e_scale;
                        for (size_t jJ = 0; jJ<m_i_sub_r; ++jJ) { // <=_i_sup_e
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
            void get_slip_contribution(double& pl_disipation_rate,
                                       double& effective_shear_rate,
                                       double* const gdot,
                                       const double* const elast_strain
                                      )
            {
               get_slip_contributions(pl_disipation_rate, effective_shear_rate, gdot,
                                      m_lattice_strain_prob.m_inv_det_v_e, elast_strain, m_kin_vals,
                                      m_slipGeom, m_kinetics, m_lattice_strain_prob);
            }

         private:

            const SlipGeom &m_slipGeom;
            const Kinetics &m_kinetics;
            const ThermoElastN &m_thermoElastN;
            const EvptnLatticeStrainProblem<ThermoElastN> m_lattice_strain_prob;
            const EvptnLatticeRotationProblem<ecmech::ntvec> m_lattice_rot_prob;

            double m_hdn_scale;
            double m_epsdot_scale_inv, m_rotincr_scale_inv;

            double m_kin_vals[Kinetics::nVals];

            const double* const m_elast_d5_n;
            const double* const m_xtal_ori_quat;
            const double* const m_def_rate_d5_sample; // d_vecds_sm would be fine too -- but do not use m_def_rate_d5_sample[iSvecS];
            const double* const m_spin_vec_sample;

            static constexpr int m_nXnDim = nDimSys * nDimSys;
            static constexpr int m_i_sub_e = 0; // ntvec
            static constexpr int m_i_sub_r = ecmech::ntvec; // nwvec
            // for mtan (material tangent stiffnes)
            double* m_mtan_sI; // null if not wanting tangent evaluation
      }; // class EvptnNRUpdstProblem

#endif

   } // namespace evptn
} // namespace ecmech

#endif // ECMECH_EVPTN_H
