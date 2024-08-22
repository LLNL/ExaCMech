// -*-c++-*-

#ifndef ECMECH_EVPTN_H
#define ECMECH_EVPTN_H

#include <cassert>

#include "ECMech_core.h"
#include "ECMech_util.h"
#include "ECMech_base_classes.h"
#include "ECMech_base_fcns.h"

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
               m_lattice_strain_prob(thermoElastN, prob_state.dt, prob_state.vNew, prob_state.eNew, prob_state.pEOS, prob_state.tkelv, prob_state.e_vecd_n),
               m_lattice_rot_prob(prob_state.dt, prob_state.quat_n),
               m_d_vecd_sm(prob_state.d_vecd_sm), // vel_grad_sm%d_vecds
               m_w_veccp_sm(prob_state.w_veccp_sm), // vel_grad_sm%w_veccp
               m_mtan_sI(nullptr)
            {
               m_hdn_scale = m_kinetics.getVals(m_kin_vals, prob_state.pEOS, prob_state.tkelv, prob_state.h_state_u);

               double adots_ref = m_kinetics.getFixedRefRate(m_kin_vals);
               double eff = vecNorm<ntvec>(m_d_vecd_sm); // do not worry about factor of sqrt(twothird)
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
             * () not necessarily safe if e_vecd is the same memory as m_e_vecd_n or quat is the same as _Cn_quat
             */
            __ecmech_hdev__
            inline
            void stateFromX(double* const e_vecd,
                            double* const quat,
                            const double* const x) {
               m_lattice_strain_prob.stateFromX(e_vecd,  &(x[m_i_sub_e]));
               m_lattice_rot_prob.stateFromX(quat,  &(x[m_i_sub_r]));
            }

            __ecmech_hdev__
            inline
            void elastNEtoC(double* const C_vecds, // nsvec
                            const double* const e_vecd_f // ntvec
                            ) const {
               m_lattice_strain_prob.elas_strain_to_cauchy_stress(C_vecds, e_vecd_f);
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
               double edot_vecd[ecmech::ntvec];
               vecsVxa<ntvec>(edot_vecd, ecmech::e_scale, &(x[m_i_sub_e]) ); // edot_vecd is now the delta, _not_ yet edot_vecd
               // e_vecd_f is end-of-step
               double e_vecd_f[ntvec];
               vecsVapb<ntvec>(e_vecd_f, edot_vecd, m_lattice_strain_prob.m_elast_dev_vec_n);
               vecsVsa<ntvec>(edot_vecd, m_lattice_strain_prob.m_inv_dt); // _now_ edot_vecd has edot_vecd
               //
               double xi_f[nwvec];
               vecsVxa<nwvec>(xi_f, ecmech::r_scale, &(x[m_i_sub_r]) );

               double A_quat[ecmech::qdim];
               emap_to_quat(A_quat, xi_f);
               //
               double C_quat[ecmech::qdim];
               get_c_quat(C_quat, A_quat, m_lattice_rot_prob.m_xtal_ori_quat_n);
               //
               double C_matx[ecmech::ndim * ecmech::ndim];
               quat_to_tensor(C_matx, C_quat);
               //
               double qr5x5_ls[ecmech::ntvec * ecmech::ntvec];
               get_rot_mat_vecd(qr5x5_ls, C_matx);

               double d_vecd_lat[ecmech::ntvec];
               double w_vec_lat[ecmech::nwvec]; // assumes nwvec = ndim

               get_xtal_frame_vel_grad_terms(d_vecd_lat, w_vec_lat, m_d_vecd_sm,  m_w_veccp_sm, C_matx, qr5x5_ls);

               //////////////////////////////
               // CALCULATIONS

               double T_vecds[ecmech::nsvec];
               m_lattice_strain_prob.elas_strain_to_kirchoff_stress(T_vecds, e_vecd_f);

               double dgdot_dtau[SlipGeom::nslip] = { 0.0 }; // crys%tmp2_slp
               double pl_vecd[ecmech::ntvec] = { 0.0 };
               double pl_wvec[ecmech::nwvec] = { 0.0 }; // \pcDhat

               get_slip_rate_terms(dgdot_dtau, pl_vecd, pl_wvec, T_vecds, m_kin_vals, m_slipGeom, m_kinetics);

               // Higher-order terms related to the elasticity stuff that's used in the residuals and jacobian calculation
               double A_e_M35[ecmech::nwvec * ecmech::ntvec];
               double ee_wvec[ecmech::nwvec];
               double ee_fac;
               elasticity_higher_order_terms(A_e_M35, ee_wvec, ee_fac, m_lattice_strain_prob.m_inv_a_vol, e_vecd_f, edot_vecd);

               // Residual Calculations
               m_lattice_strain_prob.get_elas_strain_residual(resid, m_epsdot_scale_inv, edot_vecd, pl_vecd, d_vecd_lat);
               m_lattice_rot_prob.get_omega_residual(resid, m_rotincr_scale_inv, ee_fac, xi_f, w_vec_lat, pl_wvec, ee_wvec);

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
                  double dC_quat_dxi_T[ ecmech::nwvec * ecmech::qdim ];
                  double dDsm_dxi[ ecmech::ntvec * ecmech::nwvec ];
                  double dWsm_dxi[ ecmech::nwvec * ecmech::nwvec ];
                  eval_d_dxi_impl_quat(dC_quat_dxi_T, dDsm_dxi, dWsm_dxi,
                                       m_d_vecd_sm, m_w_veccp_sm,
                                       xi_f, 
                                       m_lattice_rot_prob.m_xtal_ori_quat_n,
                                       C_matx, C_quat);

                  // d(B_S)/d(e_vecd_f)
                  //
                  m_lattice_strain_prob.template get_deriv_elast_strain_wrt_elast_strain<nDimSys>(Jacobian, dpl_deps_symm);
                  // d(B_S)/d(xi_f)
                  //
                  // jacob_er = -dDsm_dxi(:,:)
                  m_lattice_rot_prob.template get_deriv_elast_strain_wrt_omega<nDimSys>(Jacobian, dDsm_dxi);
                  // d(B_xi)/d(e_vecds_f)
                  //
                  m_lattice_strain_prob.template get_deriv_omega_wrt_elast_strain<nDimSys, m_i_sub_r>(Jacobian, edot_vecd, ee_fac, dpl_deps_skew, A_e_M35);
                  // d(B_xi)/d(xi_f)
                  //
                  m_lattice_rot_prob.template get_deriv_omega_wrt_omega<nDimSys>(Jacobian, dWsm_dxi);

                  if (m_mtan_sI) {
                     double cauchy_stress_lattice[ ecmech::nsvec ];
                     m_lattice_strain_prob.m_thermo_elast_n.getCauchy(cauchy_stress_lattice, T_vecds, m_lattice_strain_prob.m_inv_det_vol);
                     get_material_tangent_stiffness<ThermoElastN, nDimSys, m_i_sub_r>
                     (m_mtan_sI, Jacobian,
                     dC_quat_dxi_T, qr5x5_ls,
                     C_quat, C_matx,
                     cauchy_stress_lattice,
                     m_lattice_strain_prob.m_inv_det_vol,
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
                                       const double* const elas_strain
                                      )
            {
               get_slip_contributions(pl_disipation_rate, effective_shear_rate, gdot,
                                      m_lattice_strain_prob.m_inv_det_vol, elas_strain, m_kin_vals,
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

            const double* const m_d_vecd_sm; // d_vecds_sm would be fine too -- but do not use m_d_vecd_sm[iSvecS];
            const double* const m_w_veccp_sm;

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
         static const int nDimSys = ecmech::nwvec;

         __ecmech_hdev__
         RotUpdProblem(const SlipGeom& slipGeom,
                        const ThermoElastN& thermoElastN,
                        ProblemState& prob_state
                        ) :
            m_slipGeom(slipGeom),
            m_lattice_strain_prob(thermoElastN, prob_state.dt, prob_state.vNew, prob_state.eNew, prob_state.pEOS, prob_state.tkelv, prob_state.e_vecd_n),
            m_lattice_rot_prob(prob_state.dt, prob_state.quat_n),
            m_d_vecd_sm(prob_state.d_vecd_sm), // vel_grad_sm%d_vecds
            m_w_veccp_sm(prob_state.w_veccp_sm) // vel_grad_sm%w_veccp
         {

            double adots_ref = vecNorm<SlipGeom::nslip>(prob_state.gdot);

            double eff = vecNorm<ecmech::ntvec>(m_d_vecd_sm); // do not worry about factor of sqrt(twothird)
            if (eff < epsdot_scl_nzeff * adots_ref) {
                  m_epsdot_scale_inv = one / adots_ref;
            }
            else {
                  m_epsdot_scale_inv = fmin(one / eff, 1e6 * m_lattice_strain_prob.m_dt);
            }
            //
            m_rotincr_scale_inv = m_lattice_strain_prob.m_inv_dt * m_epsdot_scale_inv;

            vecsVMa<ntvec, SlipGeom::nslip>(m_pl_vecd, slipGeom.getP(), prob_state.gdot);
            vecsVMa<nwvec, SlipGeom::nslip>(m_pl_wvec, slipGeom.getQ(), prob_state.gdot);

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
            double C_quat[ecmech::qdim];
            get_c_quat(C_quat, A_quat, m_lattice_rot_prob.m_xtal_ori_quat_n);
            //
            double C_matx[ecmech::ndim * ecmech::ndim];
            quat_to_tensor(C_matx, C_quat);
            //
            double qr5x5_ls[ecmech::ntvec * ecmech::ntvec];
            get_rot_mat_vecd(qr5x5_ls, C_matx);

            double d_vecd_lat[ecmech::ntvec];
            double w_vec_lat[ecmech::nwvec]; // assumes nwvec = ndim
            get_xtal_frame_vel_grad_terms(d_vecd_lat, w_vec_lat, m_d_vecd_sm,  m_w_veccp_sm, C_matx, qr5x5_ls);

            double edot_vecd[ecmech::ntvec];
            // Calculate what this edot_vecd term should be given the current
            // state information.
            for (int i = 0; i < ecmech::ntvec; i++)
            {
                  edot_vecd[i] = m_lattice_strain_prob.m_inv_a_vol * (d_vecd_lat[i] - m_pl_vecd[i]);
            }

            // Higher-order terms related to the elasticity stuff that's used in the residuals and jacobian calculation
            double A_e_M35[ecmech::nwvec * ecmech::ntvec];
            double ee_wvec[ecmech::nwvec];
            double ee_fac;
            elasticity_higher_order_terms(A_e_M35, ee_wvec, ee_fac, m_lattice_strain_prob.m_inv_a_vol, m_lattice_strain_prob.m_elast_dev_vec_n, edot_vecd);

            // Residual Calculations
            m_lattice_rot_prob.get_omega_residual(resid, m_rotincr_scale_inv, ee_fac, xi_f, w_vec_lat, m_pl_wvec, ee_wvec);

            //////////////////////////////////////////////////////////////////////
            // JACOBIAN, fixed hardness and temperature
            //
            if (doComputeJ) {

                  //
                  //
                  // derivatives with respect to lattice orientation changes
                  double dC_quat_dxi_T[ ecmech::nwvec * ecmech::qdim ];
                  double dDsm_dxi[ ecmech::ntvec * ecmech::nwvec ];
                  double dWsm_dxi[ ecmech::nwvec * ecmech::nwvec ];
                  eval_d_dxi_impl_quat(dC_quat_dxi_T, dDsm_dxi, dWsm_dxi,
                                    m_d_vecd_sm, m_w_veccp_sm,
                                    xi_f,
                                    m_lattice_rot_prob.m_xtal_ori_quat_n,
                                    C_matx, C_quat);

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

         double m_pl_vecd[ecmech::ntvec];
         double m_pl_wvec[ecmech::nwvec]; // \pcDhat

         const double* const m_d_vecd_sm; // d_vecds_sm would be fine too -- but do not use m_d_vecd_sm[iSvecS];
         const double* const m_w_veccp_sm;

         static const int m_nXnDim = nDimSys * nDimSys;
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
                        m_lattice_strain_prob(thermoElastN, prob_state.dt, prob_state.vNew, prob_state.eNew, prob_state.pEOS, prob_state.tkelv, prob_state.e_vecd_n),
                        m_lattice_rot_prob(prob_state.dt, prob_state.quat_n),
                        m_e_vecd_n(prob_state.e_vecd_n),
                        m_C_quat(prob_state.quat_u),
                        m_d_vecd_sm(prob_state.d_vecd_sm), // vel_grad_sm%d_vecds
                        m_w_veccp_sm(prob_state.w_veccp_sm), // vel_grad_sm%w_veccp
                        m_mtan_sI(nullptr)
            {
               m_hdn_scale = m_kinetics.getVals(m_kin_vals, prob_state.pEOS, prob_state.tkelv, prob_state.h_state_u);

               double adots_ref = m_kinetics.getFixedRefRate(m_kin_vals);
               double eff = vecNorm<ntvec>(m_d_vecd_sm); // do not worry about factor of sqrt(twothird)
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
             * () not necessarily safe if e_vecd is the same memory as m_e_vecd_n or quat is the same as _Cn_quat
             */
            __ecmech_hdev__
            inline
            void stateFromX(double* const e_vecd,
                            const double* const x) {
               m_lattice_strain_prob.stateFromX(e_vecd,  &(x[m_i_sub_e]));
            }

            __ecmech_hdev__
            inline
            void elastNEtoC(double* const C_vecds, // nsvec
                            const double* const e_vecd_f // ntvec
                            ) const {
               m_lattice_strain_prob.elas_strain_to_cauchy_stress(C_vecds, e_vecd_f);
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
               double edot_vecd[ecmech::ntvec];
               vecsVxa<ntvec>(edot_vecd, ecmech::e_scale, &(x[m_i_sub_e]) ); // edot_vecd is now the delta, _not_ yet edot_vecd
               // e_vecd_f is end-of-step
               double e_vecd_f[ntvec];
               vecsVapb<ntvec>(e_vecd_f, edot_vecd, m_lattice_strain_prob.m_elast_dev_vec_n);
               vecsVsa<ntvec>(edot_vecd, m_lattice_strain_prob.m_inv_dt); // _now_ edot_vecd has edot_vecd
               //
               double C_matx[ecmech::ndim * ecmech::ndim];
               quat_to_tensor(C_matx, m_C_quat);
               //
               double qr5x5_ls[ecmech::ntvec * ecmech::ntvec];
               get_rot_mat_vecd(qr5x5_ls, C_matx);
               double d_vecd_lat[ecmech::ntvec];
               double w_vec_lat[ecmech::nwvec]; // assumes nwvec = ndim

               get_xtal_frame_vel_grad_terms(d_vecd_lat, w_vec_lat, m_d_vecd_sm,  m_w_veccp_sm, C_matx, qr5x5_ls);

               //////////////////////////////
               // CALCULATIONS

               double T_vecds[ecmech::nsvec];
               m_lattice_strain_prob.elas_strain_to_kirchoff_stress(T_vecds, e_vecd_f);

               double dgdot_dtau[SlipGeom::nslip] = { 0.0 }; // crys%tmp2_slp
               double pl_vecd[ecmech::ntvec] = { 0.0 };
               double pl_wvec[ecmech::nwvec] = { 0.0 }; // \pcDhat

               get_slip_rate_terms(dgdot_dtau, pl_vecd, pl_wvec, T_vecds, m_kin_vals, m_slipGeom, m_kinetics);

               // Residual Calculations
               m_lattice_strain_prob.get_elas_strain_residual(resid, m_epsdot_scale_inv, edot_vecd, pl_vecd, d_vecd_lat);

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

                  // d(B_S)/d(e_vecd_f)
                  //
                  m_lattice_strain_prob.template get_deriv_elast_strain_wrt_elast_strain<nDimSys>(Jacobian, dpl_deps_symm);

                  if (m_mtan_sI) {

                     // Higher-order terms related to the elasticity stuff that's used in the residuals and jacobian calculation
                     double A_e_M35[ecmech::nwvec * ecmech::ntvec];
                     double ee_wvec[ecmech::nwvec];
                     double ee_fac;
                     double xi_f[ecmech::nwvec] = {};

                     m_lattice_rot_prob.deltaOmegaFromState(xi_f, m_C_quat, m_lattice_rot_prob.m_xtal_ori_quat_n);

                     elasticity_higher_order_terms(A_e_M35, ee_wvec, ee_fac, m_lattice_strain_prob.m_inv_a_vol, e_vecd_f, edot_vecd);

                     // derivatives with respect to lattice orientation changes
                     double dC_quat_dxi_T[ ecmech::nwvec * ecmech::qdim ];
                     double dDsm_dxi[ ecmech::ntvec * ecmech::nwvec ];
                     double dWsm_dxi[ ecmech::nwvec * ecmech::nwvec ];
                     eval_d_dxi_impl_quat(dC_quat_dxi_T, dDsm_dxi, dWsm_dxi,
                                          m_d_vecd_sm, m_w_veccp_sm,
                                          xi_f,
                                          m_lattice_rot_prob.m_xtal_ori_quat_n,
                                          C_matx, m_C_quat);

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
                     // d(B_xi)/d(e_vecds_f)
                     //
                     m_lattice_strain_prob.template get_deriv_omega_wrt_elast_strain<nDimSolve, m_i_sub_r>(Jacobian2, edot_vecd, ee_fac, dpl_deps_skew, A_e_M35);
                     // d(B_xi)/d(xi_f)
                     //
                     m_lattice_rot_prob.template get_deriv_omega_wrt_omega<nDimSolve>(Jacobian2, dWsm_dxi);

                     double cauchy_stress_lattice[ ecmech::nsvec ];
                     m_lattice_strain_prob.m_thermo_elast_n.getCauchy(cauchy_stress_lattice, T_vecds, m_lattice_strain_prob.m_inv_det_vol);
                     get_material_tangent_stiffness<ThermoElastN, nDimSolve, m_i_sub_r>(m_mtan_sI, Jacobian2,
                     dC_quat_dxi_T, qr5x5_ls,
                     m_C_quat, C_matx,
                     cauchy_stress_lattice,
                     m_lattice_strain_prob.m_inv_det_vol,
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
                                       const double* const elas_strain
                                      )
            {
               get_slip_contributions(pl_disipation_rate, effective_shear_rate, gdot,
                                      m_lattice_strain_prob.m_inv_det_vol, elas_strain, m_kin_vals,
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

            const double* const m_e_vecd_n;
            const double* const m_C_quat;
            const double* const m_d_vecd_sm; // d_vecds_sm would be fine too -- but do not use m_d_vecd_sm[iSvecS];
            const double* const m_w_veccp_sm;

            static const int m_nXnDim = nDimSys * nDimSys;
            static const int m_i_sub_e = 0; // ntvec
            static const int m_i_sub_r = ecmech::ntvec; // nwvec
            // for mtan (material tangent stiffnes)
            double* m_mtan_sI; // null if not wanting tangent evaluation
      }; // class EvptnNRUpdstProblem

#endif

   } // namespace evptn
} // namespace ecmech

#endif // ECMECH_EVPTN_H
