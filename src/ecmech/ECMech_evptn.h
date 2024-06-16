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
               _lattice_strain_prob(thermoElastN, dt, detV, eVref, p_EOS, tK, e_vecd_n),
               _lattice_rot_prob(dt, Cn_quat),
               _eVref(eVref),
               _p_EOS(p_EOS),
               _tK(tK),
               _h_state(h_state),
               _d_vecd_sm(d_vecd_sm), // vel_grad_sm%d_vecds
               _w_veccp_sm(w_veccp_sm), // vel_grad_sm%w_veccp
               _mtan_sI(nullptr)
            {
               _hdn_scale = _kinetics.getVals(_kin_vals, _p_EOS, _tK, _h_state);

               double adots_ref = _kinetics.getFixedRefRate(_kin_vals);
               double eff = vecNorm<ntvec>(_d_vecd_sm); // do not worry about factor of sqrt(twothird)
               if (eff < epsdot_scl_nzeff * adots_ref) {
                  _epsdot_scale_inv = one / adots_ref;
               }
               else {
                  _epsdot_scale_inv = fmin(one / eff, 1e6 * _lattice_strain_prob.m_dt);
               }
               //
               _rotincr_scale_inv = _lattice_strain_prob.m_inv_dt * _epsdot_scale_inv;
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
            double getDtRi() const { return _lattice_strain_prob.m_inv_dt; }

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
               _lattice_rot_prob.stateFromX(quat,  &(x[_i_sub_r]));
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
                  for (size_t ijJ = 0; ijJ<_nXnDim; ++ijJ) {
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
               vecsVapb<ntvec>(e_vecd_f, edot_vecd, _lattice_strain_prob.m_elast_dev_vec_n);
               vecsVsa<ntvec>(edot_vecd, _lattice_strain_prob.m_inv_dt); // _now_ edot_vecd has edot_vecd
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
               get_c_quat(C_quat, A_quat, _lattice_rot_prob.m_xtal_ori_quat_n);
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
               double pl_vecd[ecmech::ntvec] = { 0.0 };
               double pl_wvec[ecmech::nwvec] = { 0.0 }; // \pcDhat

               get_slip_rate_terms(dgdot_dtau, pl_vecd, pl_wvec, T_vecds, _kin_vals, _slipGeom, _kinetics);

               // Higher-order terms related to the elasticity stuff that's used in the residuals and jacobian calculation
               double A_e_M35[ecmech::nwvec * ecmech::ntvec];
               double ee_wvec[ecmech::nwvec];
               double ee_fac;
               elasticity_higher_order_terms(A_e_M35, ee_wvec, ee_fac, _lattice_strain_prob.m_inv_a_vol, e_vecd_f, edot_vecd);

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
                  get_slip_rate_deriv_terms(dpl_deps_symm, dpl_deps_skew, dgdot_dtau, _lattice_strain_prob.m_inv_a_vol, _slipGeom, _lattice_strain_prob.m_thermo_elast_n);

                  //
                  //
                  // derivatives with respect to lattice orientation changes
                  double dC_quat_dxi_T[ ecmech::nwvec * ecmech::qdim ];
                  double dDsm_dxi[ ecmech::ntvec * ecmech::nwvec ];
                  double dWsm_dxi[ ecmech::nwvec * ecmech::nwvec ];
                  eval_d_dxi_impl_quat(dC_quat_dxi_T, dDsm_dxi, dWsm_dxi,
                                       _d_vecd_sm, _w_veccp_sm,
                                       xi_f, 
                                       _lattice_rot_prob.m_xtal_ori_quat_n,
                                       C_matx, C_quat);

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
                     _lattice_strain_prob.m_thermo_elast_n.getCauchy(cauchy_stress_lattice, T_vecds, _lattice_strain_prob.m_inv_det_vol);
                     get_material_tangent_stiffness<ThermoElastN, nDimSys, _i_sub_r>(_mtan_sI, Jacobian,
                                                                                     dC_quat_dxi_T, qr5x5_ls,
                                                                                     C_quat, C_matx,
                                                                                     cauchy_stress_lattice,
                                                                                     _lattice_strain_prob.m_inv_det_vol,
                                                                                     _lattice_strain_prob.m_inv_a_vol,
                                                                                     _lattice_strain_prob.m_thermo_elast_n);
                  }

                  // SCALING
                  {
                     double scaleFactorJ;
                     for (size_t iJ = 0; iJ<_i_sub_r; ++iJ) {
                        // Jacobian(i_sub_e:i_sup_e,i_sub_e:i_sup_e) = jacob_ee * epsdot_scale_inv  * e_scale ! resid, x
                        scaleFactorJ = _epsdot_scale_inv * ecmech::e_scale;
                        for (size_t jJ = 0; jJ<_i_sub_r; ++jJ) { // <=_i_sup_e
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
                        for (size_t jJ = 0; jJ<_i_sub_r; ++jJ) { // <=_i_sup_e
                           int ijJ = ECMECH_NN_INDX(iJ, jJ, nDimSys);
                           Jacobian[ ijJ ] *= scaleFactorJ;
                        }

                        // Jacobian(i_sub_r:i_sup_r,i_sub_r:i_sup_r) = jacob_rr * rotincr_scale_inv * r_scale
                        scaleFactorJ = _rotincr_scale_inv * ecmech::r_scale;
                        for (size_t jJ = _i_sub_r; jJ<nDimSys; ++jJ) { // <_i_sup_r
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
                                      _lattice_strain_prob.m_inv_det_vol, elas_strain, _kin_vals,
                                      _slipGeom, _kinetics, _lattice_strain_prob);
            }
                              

         private:

            const SlipGeom &_slipGeom;
            const Kinetics &_kinetics;
            const EvptnLatticeStrainProblem<ThermoElastN> _lattice_strain_prob;
            const EvptnLatticeRotationProblem<ecmech::ntvec> _lattice_rot_prob;

            double _eVref, _p_EOS, _tK, _a_V;

            double _hdn_scale;
            double _epsdot_scale_inv, _rotincr_scale_inv;

            double _kin_vals[Kinetics::nVals];

            const double* const _h_state;
            const double* const _d_vecd_sm; // d_vecds_sm would be fine too -- but do not use _d_vecd_sm[iSvecS];
            const double* const _w_veccp_sm;

            static constexpr size_t _nXnDim = nDimSys * nDimSys;
            static constexpr size_t _i_sub_e = 0; // ntvec
            static constexpr size_t _i_sub_r = ecmech::ntvec; // nwvec

            // for mtan (material tangent stiffnes)
            double* _mtan_sI; // null if not wanting tangent evaluation
      }; // class EvptnUpdstProblem
   } // namespace evptn
} // namespace ecmech

#endif // ECMECH_EVPTN_H
