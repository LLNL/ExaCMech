// -*-c++-*-

#ifndef ECMECH_EVPTN_FI_H
#define ECMECH_EVPTN_FI_H

#include <cassert>

#include "ECMech_core.h"
#include "ECMech_util.h"
#include "ECMech_eosSimple.h"

#include "SNLS_lup_solve.h"
#include "SNLS_TrDLDenseG.h"

#include "RAJA/RAJA.hpp"

namespace ecmech {
   namespace evptn {
      template<class SlipGeom, class Kinetics, class ThermoElastN>
      class EvptnUpdstFIProblem
      {
         public:

            static const int nDimSys = ecmech::ntvec + ecmech::nwvec + Kinetics::nH;

            // constructor
            __ecmech_hdev__
            EvptnUpdstFIProblem(const SlipGeom& slipGeom,
                              const Kinetics& kinetics,
                              const ThermoElastN& thermoElastN,
                              double dt,
                              double detV, double eVref, double p_EOS, double tK,
                              const double* const h_state_n,
                              const double* const e_vecd_n,
                              const double* const Cn_quat,
                              const double* const d_vecd_sm, // okay to pass d_vecds_sm, but d_vecd_sm[iSvecS] is not used
                              const double* const w_veccp_sm,
                              const double hdn,
                              const double adots_ref
                              )
               : _slipGeom(slipGeom),
               _kinetics(kinetics),
               _thermoElastN(thermoElastN),
               _dt(dt),
               _detV(detV),
               _eVref(eVref),
               _p_EOS(p_EOS),
               _tK(tK),
               _h_state_n(h_state_n),
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
                // This is based on the initial state values
                // with this being a fully coupled solve of the hardness,
                // elastic strains, and lattice orientations these will change
                // every iteration.
               _hdn_scale = hdn;
               double eff = vecNorm<ntvec>(_d_vecd_sm); // do not worry about factor of sqrt(twothird)
               if (eff < epsdot_scl_nzeff * adots_ref) {
                  _epsdot_scale_inv = one / adots_ref;
               }
               else {
                  _epsdot_scale_inv = fmin(one / eff, 1e6 * _dt);
               }
               //
               _rotincr_scale_inv = _dt_ri * _epsdot_scale_inv;
               for(int iH = 0; iH < _kinetics.nH; iH++)
               {
                    _x_scale[iH] = fmax(_h_state_n[iH], 1.0); // TO_DO -- generalize this to not max with 1
                    // NOTE : see comment below about changing Jacobian calculation if _res_scale != one / s_scale
                    _res_scale[iH] = one / _x_scale[iH];
               }
            }

            // deconstructor
            __ecmech_hdev__
            ~EvptnUpdstFIProblem() {}

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
            double getShrateEff() const { return _shrate_eff_contrib; }

            __ecmech_hdev__
            inline
            double getDisRate() const { return _dp_dis_rate_contrib; }

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
                            double* const hard,
                            const double* const x) {
               double e_vecd_delta[ecmech::ntvec];
               vecsVxa<ntvec>(e_vecd_delta, ecmech::e_scale, &(x[_i_sub_e]) );
               vecsVapb<ntvec>(e_vecd, e_vecd_delta, _e_vecd_n);

               double xi_f[nwvec];
               vecsVxa<nwvec>(xi_f, ecmech::r_scale, &(x[_i_sub_r]) );
               //
               double A_quat[ecmech::qdim];
               emap_to_quat(A_quat, xi_f);
               //
               // double C_quat[ecmech::qdim] ;
               // get_c_quat(C_quat, A_quat, _Cn_quat) ;
               get_c_quat(quat, A_quat, _Cn_quat);
                // Obtain the final H update so if h was transformed within kinetics
                // we want the transformed back form of it.
                _kinetics.getHUpdate(_h_state_n, &x[_i_sub_h], _x_scale, hard, true);

            }

            __ecmech_hdev__
            inline
            void elastNEtoT(double* const T_vecds, // nsvec
                            const double* const e_vecd_f // ntvec
                            ) {
               //// do not need to use elaw_T_BT here as T and BT are the same
               //
               // specialize to cem%l_lin_lnsd
               // CALL elawn_T(s_meas, e_vecd_f, crys%elas, tK, .TRUE., a_V, &
               // & p_EOS, eVref, crys%i_eos_model, crys%eos_const &
               // &)
               double Ee_vecds[ecmech::nsvec];
               vecsVxa<ntvec>(Ee_vecds, _a_V_ri, e_vecd_f);
               //// tr_Ee = three * DLOG(a_V%r)
               //// CALL trace_to_vecds_s(s_meas%Ee_vecds(SVEC), tr_Ee)
               Ee_vecds[iSvecS] = sqr3 * log(_a_V); // could go into constructor
               //
               //// Kirchhoff stress from Ee_vecds
               // CALL elawn_lin_op(s_meas%T_vecds, s_meas%Ee_vecds, cem, tK, &
               // & p_EOS, eVref, i_eos_model, eos_const)
               _thermoElastN.eval(T_vecds, Ee_vecds, _tK, _p_EOS, _eVref);
            }

            __ecmech_hdev__
            inline
            void elastNEtoC(double* const C_vecds, // nsvec
                            const double* const e_vecd_f // ntvec
                            ) {
               double T_vecds[ecmech::nsvec];
               this->elastNEtoT(T_vecds, e_vecd_f);
               _thermoElastN.getCauchy(C_vecds, T_vecds, _detV_ri);
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
               double hard[Kinetics::nH];
               double vals_extra[Kinetics::nValsDerivs];
               //
               // Obtain the H update for the solver if h was transformed within kinetics
               // we don't want the transformed back form.
               _kinetics.getHUpdate(_h_state_n, &x[_i_sub_h], _x_scale, hard, false);
               //
               // Need to update kinetic values each iteration
               _hdn_scale = _kinetics.getVals(_kin_vals, _p_EOS, _tK, hard, vals_extra);
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
               vecsVMTa<ntvec>(d_vecd_lat, qr5x5_ls, _d_vecd_sm);
               // d_vecds_lat(SVEC) = vel_grad_sm%d_vecds(SVEC)
               //
               //// CALL rot_mat_vecd(A, qr5x5_A)
               //
               // CALL rot_mat_wveccp(C_matx, qr3x3_ls) // amounts to qr3x3_ls = C_matx
               // CALL matt_x_vec_3(qr3x3_ls, vel_grad_sm%w_veccp, w_vec_lat)
               double w_vec_lat[ecmech::nwvec]; // assumes nwvec = ndim
               vecsVMTa<ndim>(w_vec_lat, C_matx, _w_veccp_sm);

               //////////////////////////////
               // CALCULATIONS

               double T_vecds[ecmech::nsvec];
               this->elastNEtoT(T_vecds, e_vecd_f);
               //
               double taua[SlipGeom::nslip] = { 0.0 }; // crys%tmp4_slp
               double dgdot_dtau[SlipGeom::nslip] = { 0.0 }; // crys%tmp2_slp
               // Would have dimensions of nslip x nH
               double dgdot_dh[SlipGeom::nslip * Kinetics::nH] = { 0.0 }; // crys%tmp3_slp
               // Would have dimensions of nH x nslip
               double dh_dgdot[SlipGeom::nslip * Kinetics::nH] = { 0.0 };
               double pl_vecd[ecmech::ntvec] = { 0.0 };
               double pl_wvec[ecmech::nwvec] = { 0.0 }; // \pcDhat
               double dhdot_dh[Kinetics::nH * Kinetics::nH] = { 0.0 };
               double hdot[Kinetics::nH] = { 0.0 };
               if (SlipGeom::nslip > 0) {
                  // resolve stress onto slip systems
                  // CALL resolve_tau_a_n(crys%tmp4_slp, s_meas%T_vecds, crys)
                  vecsVaTM<ntvec, SlipGeom::nslip>(taua, T_vecds, _slipGeom.getP() );
                  //
                  // CALL plaw_eval(pl_vecd, pl_wvec, gss, crys, tK, ierr)
                  _kinetics.evalGdots(_gdot, dgdot_dtau, dgdot_dh, taua, _kin_vals, true, vals_extra);
                  //
                  // CALL sum_slip_def(pl_vecd, pl_wvec, crys%tmp1_slp, crys) ;
                  vecsVMa<ntvec, SlipGeom::nslip>(pl_vecd, _slipGeom.getP(), _gdot);
                  vecsVMa<nwvec, SlipGeom::nslip>(pl_wvec, _slipGeom.getQ(), _gdot);
                  // dgdot_dh may or may not be scaled by the below set of code to account for
                  // differences from the evaluation within _kinetics.evalGdots
                  _kinetics.getExtDerivs(hdot, dhdot_dh, dh_dgdot, dgdot_dh, hard, _gdot);
               }
               //
               //// shrate_l%gdot => crys%tmp1_slp

               // from e edot product term in spin (formerly neglected)
               //
               double A_e_M35[ecmech::nwvec * ecmech::ntvec];
               M35_d_AAoB_dA(A_e_M35, e_vecd_f);
               //
               double ee_wvec[ecmech::nwvec];
               vecsVMa<nwvec, ntvec>(ee_wvec, A_e_M35, edot_vecd);
               //
               double ee_fac = onehalf * _a_V_ri * _a_V_ri;

               // RESIDUAL B_S
               //
               for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                  resid[_i_sub_e + iTvec] = _epsdot_scale_inv * ( // SCALING
                     _a_V_ri * edot_vecd[iTvec] + pl_vecd[iTvec] - d_vecd_lat[iTvec]);
               }

               // RESIDUAL B_xi
               //
               for (int iWvec = 0; iWvec < ecmech::nwvec; ++iWvec) {
                  resid[_i_sub_r + iWvec] = _rotincr_scale_inv * // SCALING
                                            (xi_f[iWvec] - _dt * (w_vec_lat[iWvec] - pl_wvec[iWvec] + ee_fac * ee_wvec[iWvec]) );
               }

                for (int iH = 0; iH < Kinetics::nH; iH++) {
                    resid[_i_sub_h + iH] = (x[_i_sub_h + iH] * _x_scale[iH] - hdot[iH] * _dt) * _res_scale[iH];
                }

               _dp_dis_rate_contrib = zero;
               if (SlipGeom::nslip > 0) {
                  // CALL calc_pl_dis(dp_dis_rate_contrib, crys%tmp4_slp, crys%tmp1_slp, detV%ri)
                  _dp_dis_rate_contrib = _detV_ri * vecsyadotb<SlipGeom::nslip>(taua, _gdot);
               }

               //// need shrate%eff instead
               //// CALL calc_pl_eff(dp_def_rate_contrib, pl_vecd, detV%ri)
               // CALL setup_ss_shrate_vals(shrate_l, crys%tmp1_slp, zero, .TRUE.)
#if defined(ECMECH_USE_DPEFF)
               _shrate_eff_contrib = vecd_Deff(pl_vecd);
#else
               _shrate_eff_contrib = vecsssumabs<SlipGeom::nslip>(_gdot);
#endif
               //////////////////////////////////////////////////////////////////////
               // JACOBIAN, temperature
               //
               //   The Jacobian can be broken down into the following
               //   PDE blocks. For complex hardening models with a large number
               //   of hardening variables this can make the Jacobian quite large.
               //   dR_E/dX_E | dR_E/dX_W  | dR_E/dX_h 
               //   ----------|------------|----------
               //   dR_W/dX_E | dR_W/dX_W  | dR_W/dX_h
               //   ----------|------------|----------
               //   dR_H/dX_E | dR_H/dX_W  | dR_H/dX_h
               //
               //   For the current models within ExaCMech, dR_H / dX_W equals 0.

               if (doComputeJ) {
                  // use RAJA::View machinery to simplify indexing for blocks in the Jacobian matrix ;
                  // can always swap this out later if it ends up being too heavyweight ;
                  // RAJA defaults to "row-major" -- final dimension indexing the fastest
                  //
                  const int JDIM = 2;
                  //
                  // preliminaries
                  //
                  double dpl_deps_symm[ ecmech::ntvec * ecmech::ntvec ] = { 0.0 };
                  double dpl_deps_skew[ ecmech::nwvec * ecmech::ntvec ] = { 0.0 };
                  // We've lifted this out of the loop, since it's used in the hardening
                  // set of equations as well. 
                  double dgdot_deps[ ecmech::ntvec * SlipGeom::nslip ] = { 0.0 };
                  if (SlipGeom::nslip > 0) {
                     // CALL elawn_T_dif(s_meas, e_vecd_f, crys%elas, tK, a_V, &
                     // & p_EOS, eVref, crys%i_eos_model, crys%eos_const, &
                     // & dpEOS_dtK, .FALSE., .FALSE., .FALSE.)
                     // CALL eval_dtaua_deps_n(dtaua_deps, s_meas%dT_deps, crys)
                     //
                     double dtaua_deps[ ecmech::ntvec * SlipGeom::nslip ];
                     _thermoElastN.multDTDepsT(dtaua_deps, _slipGeom.getP(), _a_V_ri, SlipGeom::nslip);

                     // CALL plaw_eval_dif_sn(TVEC, &
                     // & dpl_deps_symm, dpl_deps_skew, dgdot_deps, &
                     // & dtaua_deps, gss, crys, s_meas, .FALSE.)
                     // This can be re-arranged to be the transpose operation
                     // Since all the matrix operations down below are nslip x ntvec
                     // rather than ntvec x nslip
                     for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                        for (int iSlip = 0; iSlip < SlipGeom::nslip; ++iSlip) {
                           int ijThis = ECMECH_NM_INDX(iTvec, iSlip, ecmech::ntvec, SlipGeom::nslip);
                           dgdot_deps[ijThis] = dgdot_dtau[iSlip] * dtaua_deps[ijThis];
                        }
                     }

                     // DO islip = 1, crys%nslip
                     // DO i_TVEC = 1, TVEC
                     // dpl_deps_symm(:,i_TVEC) +=
                     // & crys%P_ref_vec(:,islip) * dgdot_deps(i_TVEC,islip)
                     // dpl_deps_skew(:,i_TVEC) +=
                     // & crys%Q_ref_vec(:,islip) * dgdot_deps(i_TVEC,islip)
                     // END DO
                     // END DO
                     vecsMABT<ntvec, SlipGeom::nslip>(dpl_deps_symm, _slipGeom.getP(), dgdot_deps);
                     vecsMABT<nwvec, ntvec, SlipGeom::nslip>(dpl_deps_skew, _slipGeom.getQ(), dgdot_deps);
                  }
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
                  {
                     RAJA::View<double, RAJA::Layout<JDIM> > jacob_ee(Jacobian, nDimSys, nDimSys);

                     // dislocation plasticity;
                     // first contribution; overwrite
                     // fixme:
                     // Swap the ordering of these
                     for (int jTvec = 0; jTvec < ecmech::ntvec; ++jTvec) {
                        for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                           jacob_ee(iTvec, jTvec) = dpl_deps_symm[ECMECH_NN_INDX(iTvec, jTvec, ecmech::ntvec)];
                        }
                     }

                     // elastic rate
                     //
                     {
                        double adti = _a_V_ri * _dt_ri;
                        for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                           jacob_ee(iTvec, iTvec) += adti;
                        }
                     }
                  }

                  // d(B_S)/d(xi_f)
                  //
                  // jacob_er = -dDsm_dxi(:,:)
                  {
                     RAJA::View<double, RAJA::Layout<JDIM> > jacob_er(Jacobian, nDimSys, nDimSys);
                     //fixme:
                     //Swap ordering of these
                     for (int jWvec = 0; jWvec<ecmech::nwvec; ++jWvec) {
                        for (int iTvec = 0; iTvec<ecmech::ntvec; ++iTvec) {
                           // could also make dDsm_dxi into a RAJA view, but not really needed
                           jacob_er(iTvec,
                                    jWvec + _i_sub_r) = -dDsm_dxi[ ECMECH_NM_INDX(iTvec, jWvec, ecmech::ntvec, ecmech::nwvec) ];
                        }
                     }
                  }

                  //d(B_S) / d(h)
                  //
                  // jacob_eh = (dR_e/dD^p)(dD^p / dgdot)(dgdot / dh)
                  {
                      //dR_e / dD^p  => I_5x5       (5 x 5 Identity matrix) 
                      //dD^p / dgdot => [P_sym ...] (5 x nslip)
                      //dgdot / dh => [ ... ]       (nslip x nh)
                      //Just a bunch of matrix products to get down to the
                      //5 x nh matrix 
                      double dR_e_dH[ecmech::ntvec * Kinetics::nH];
                      vecsMAB<ecmech::ntvec, Kinetics::nH, SlipGeom::nslip>(dR_e_dH, _slipGeom.getP(), dgdot_dh);
                      RAJA::View<double, RAJA::Layout<JDIM> > jacob_eh(Jacobian, nDimSys, nDimSys);
                      for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                        for (int jH = 0; jH < Kinetics::nH; ++jH) {
                           // could also make dDsm_dxi into a RAJA view, but not really needed
                           jacob_eh(iTvec, jH + _i_sub_h) = dR_e_dH[ ECMECH_NM_INDX(iTvec, jH, ecmech::ntvec, Kinetics::nH) ];
                        }
                     }
                  } 

                  // d(B_xi)/d(e_vecds_f)
                  //
                  {
                     RAJA::View<double, RAJA::Layout<JDIM> > jacob_re(Jacobian, nDimSys, nDimSys);

                     double A_edot_M35[ecmech::nwvec * ecmech::ntvec];
                     M35_d_AAoB_dA(A_edot_M35, edot_vecd);

                     double dt_ee_fac = _dt * ee_fac;

                     for (int iWvec = 0; iWvec < ecmech::nwvec; ++iWvec) {
                        for (int jTvec = 0; jTvec < ecmech::ntvec; ++jTvec) {
                           int ijWT = ECMECH_NM_INDX(iWvec, jTvec, ecmech::nwvec, ecmech::ntvec);
                           jacob_re(iWvec + _i_sub_r, jTvec) =
                              _dt * dpl_deps_skew[ijWT] - dt_ee_fac * (A_e_M35[ijWT] * _dt_ri - A_edot_M35[ijWT]);
                        }
                     }
                  }

                  // d(B_xi)/d(xi_f)
                  //
                  {
                     RAJA::View<double, RAJA::Layout<JDIM> > jacob_rr(Jacobian, nDimSys, nDimSys);

                     for (int iWvec = 0; iWvec < ecmech::nwvec; ++iWvec) {
                        for (int jWvec = 0; jWvec < ecmech::nwvec; ++jWvec) {
                           int ijWW = ECMECH_NN_INDX(iWvec, jWvec, ecmech::nwvec);
                           jacob_rr(iWvec + _i_sub_r, jWvec + _i_sub_r) = -_dt * dWsm_dxi[ijWW];
                        }

                        jacob_rr(iWvec + _i_sub_r, iWvec + _i_sub_r) += one;
                     }
                  }

                  // d(B_xi)/d(H)
                  // jacob_rh = (dR_xi/dW^p)(dW^p / dgdot)(dgdot / dh) * dt
                  {
                      //dR_xi / dW^p  => I_3x3       (3 x 3 Identity matrix) 
                      //dW^p / dgdot => [W_sym ...] (3 x nslip)
                      //dgdot / dh => [ ... ]       (nslip x nh)
                      //Just a bunch of matrix products to get down to the
                      //3 x nh matrix 
                      double dR_xi_dX_h[ecmech::nwvec * Kinetics::nH];
                      vecsMAB<ecmech::nwvec, Kinetics::nH, SlipGeom::nslip>(dR_xi_dX_h, _slipGeom.getQ(), dgdot_dh);
                      RAJA::View<double, RAJA::Layout<JDIM> > jacob_rh(Jacobian, nDimSys, nDimSys);
                      for (int iWvec = 0; iWvec < ecmech::nwvec; ++iWvec) {
                        for (int jH = 0; jH < Kinetics::nH; ++jH) {
                           // could also make dDsm_dxi into a RAJA view, but not really needed
                           jacob_rh(iWvec + _i_sub_r, jH + _i_sub_h) = _dt * dR_xi_dX_h[ ECMECH_NM_INDX(iWvec, jH, ecmech::nwvec, Kinetics::nH) ];
                        }
                     }
                  }

                  // d(B_h) / d(e)
                  // jacob_he = dt * (dhdot/dgdot)(dgdot/de)
                  {
                     // nh x ntvec matrix
                     double dR_h_dX_e[Kinetics::nH * ecmech::ntvec];
                     vecsMABT<Kinetics::nH, ecmech::ntvec, SlipGeom::nslip>(dR_h_dX_e, dh_dgdot, dgdot_deps);
                     RAJA::View<double, RAJA::Layout<JDIM> > jacob_he(Jacobian, nDimSys, nDimSys);
                     for (int iH = 0; iH < Kinetics::nH; ++iH) {
                        for (int jE = 0; jE < ecmech::ntvec; ++jE) {
                           // could also make dDsm_dxi into a RAJA view, but not really needed
                           jacob_he(iH + _i_sub_h, jE) = -_dt * dR_h_dX_e[ ECMECH_NM_INDX(iH, jE, Kinetics::nH, ecmech::ntvec) ];
                        }
                     }
                  }

                  // d(B_h) / d(h)
                  // jacob_hh = I + (dhdot/dh)
                  {
                     // We can go ahead and scale things here, since we won't be including it in a potential
                     // tangent modulus calculation down below
                     // Multiply dsdot_ds terms by the negative outer product of x_scale and res_scale and dt
                     for (int i = _i_sub_h; i < nDimSys; i++) {
                        for (int j = _i_sub_h; j < nDimSys; j++) {
                           Jacobian[ECMECH_NN_INDX(i, j, nDimSys)] = dhdot_dh[ECMECH_NN_INDX(i - _i_sub_h, j - _i_sub_h, Kinetics::nH)] 
                                                                     * (-_x_scale[j - _i_sub_h] * _res_scale[i - _i_sub_h] * _dt);
                        }
                     }

                     // Now add in the identity term
                     // The below is based on the assumption that _res_scale = 1/_x_scale
                     // if this were to change in the future version than this would need to become
                     // Jacobian[ECMECH_NN_INDX(i, i, nDimSys)] += ecmech::one * _x_scale[i] * _r_scale[i]
                     for (int i = _i_sub_h; i < nDimSys; i++) {
                        Jacobian[ECMECH_NN_INDX(i, i, nDimSys)] += ecmech::one;
                     }
                  } 

                  if (_mtan_sI) { // l_eval_derivs
                     //
                     // material tangent, do before scaling of Jacobian
                     // eval_mtan()

                     // must solve a set of systems to get needed partial derivatives

                     // eval_mtan_pfrac_r(de_dI, dxi_dI)
                     // compared to Fortran coding, dI has reduced back down to being only the deviatoric part of the deformation rate ; UB_I = ntvec
                     const int nRHS = ecmech::ntvec;
                     static const int nsys2 = ecmech::ntvec + ecmech::nwvec; 
                     double pfrac_rhs_T[ nRHS * nsys2 ] = { 0.0 }; // transpose for use in SNLS_LUP_SolveX !
                     // derivatives end up in pfrac_rhs_T
                     // de_dI  is pfrac_rhs_T[ :, _i_sub_e:i_sup_e ] // ecmech::ntvec * ecmech::ntvec
                     // dxi_dI is pfrac_rhs_T[ :, _i_sub_r:i_sup_r ] // ecmech::nwvec * ecmech::ntvec
                     {
                        // RHS
                        //
                        {
                           // negatives cancel
                           // pfrac_rhs(i_sub_e:i_sup_e,I_subD:I_supD) = TRANSPOSE(qr5x5_ls)
                           for (int jE = 0; jE<nRHS; ++jE) {
                              for (int iE = 0; iE<ntvec; ++iE) { // ntvec, _not_ nDimSys // iE is same as index into nDimSys, give how d(resid)/d(_d_vecd_sm) works out
                                 // pfrac_rhs[ECMECH_NM_INDX(iE,jE,nDimSys,nRHS)] = qr5x5_ls[ECMECH_NN_INDX(jE,iE,ntvec)] ;
                                 pfrac_rhs_T [ECMECH_NM_INDX(jE, iE, nRHS, nsys2)] = qr5x5_ls[ECMECH_NN_INDX(jE, iE, ntvec)];
                              }
                           }
                        }

                        // SYSTEM
                        //
                        // Only want a fraction of the system to be consistent with
                        // the other ways of things.
                        // So, we should only take the portion that elastic strain
                        // and exponential mapping parts of the Jacobian down here.
                        // Anything else will contain the hardening contributions
                        // which we don't necessarily want for our current implementation
                        double pfrac_sys[ nsys2 * nsys2 ];

                        for (int i_jac = 0; i_jac < nsys2; i_jac++) {
                           for (int j_jac = 0; j_jac < nsys2; j_jac++) {
                              pfrac_sys[ECMECH_NN_INDX(i_jac, j_jac, nsys2)] = Jacobian[ECMECH_NN_INDX(i_jac, j_jac, nDimSys)];
                           }
                        }

                        int err = SNLS_LUP_SolveX<(ecmech::ntvec + ecmech::nwvec)>(pfrac_sys, pfrac_rhs_T, nRHS);
                        if (err != 0) {
                           ECMECH_FAIL(__func__, "error from SNLS_LUP_SolveX");
                        }
                     } // eval_mtan_pfrac_r

                     // in Fortran code dCn_quat_dI was used to store dC_quat_dI, but here just call it dC_quat_dI
                     double dC_quat_dI[ ecmech::qdim * nRHS ];
                     for (int ii_I = 0; ii_I<nRHS; ++ii_I) {
                        for (int ii_Q = 0; ii_Q<ecmech::qdim; ++ii_Q) {
                           int iiQI = ECMECH_NM_INDX(ii_Q, ii_I, ecmech::qdim, nRHS);
                           dC_quat_dI[iiQI] = 0.0;
                           for (int ii_W = 0; ii_W<ecmech::nwvec; ++ii_W) {
                              // dC_quat_dI[iiQI] += dC_quat_dxi_T(ii_W,ii_Q) * dxi_dI(ii_W,ii_I)
                              dC_quat_dI[iiQI] +=
                                 dC_quat_dxi_T[ECMECH_NM_INDX(ii_W, ii_Q, ecmech::nwvec, ecmech::qdim)] *
                                 pfrac_rhs_T[ECMECH_NM_INDX(ii_I, _i_sub_r + ii_W, nRHS, nsys2)];
                           }
                        }
                     }


                     // contribution through e
                     //
                     // double temp_M6I[ ecmech::nsvec*nRHS ] ; // (SVEC,UB_I)
                     double temp_M6[ ecmech::nsvec2 ]; // (SVEC,UB_I)
                     // dsigClat_def(:,:) = detVi * s_meas%dT_deps(:,:)
                     // temp_M6I = MATMUL(dsigClat_def(:,1:TVEC), de_dI(:,:))
                     {
                        // TODO : get rid of the need for this memory copy (with transpose) into de_dI
                        double de_dI[ ecmech::ntvec * nRHS ]; // nRHS=ecmech::ntvec, but put nRHS in here for clarity, and thus use ECMECH_NM_INDX instead of ECMECH_NN_INDX
                        for (int iTvec = 0; iTvec<ecmech::ntvec; ++iTvec) {
                           for (int jTvec = 0; jTvec<nRHS; ++jTvec) {
                              de_dI[ECMECH_NM_INDX(iTvec, jTvec, ecmech::ntvec,
                                                   nRHS)] = pfrac_rhs_T[ECMECH_NM_INDX(jTvec, iTvec, nRHS, nsys2)];
                           }
                        }

                        _thermoElastN.multCauchyDif(temp_M6, de_dI, _detV_ri, _a_V_ri);
                     }
                     //
                     // CALL qr6x6_pre_mul(mtan_sI, temp_M6I, qr5x5_ls, UB_I, .FALSE.)
                     // UB_I=nRHS ; but here do nsvec instead of nRHS so that there is less monkeying with memory later
                     qr6x6_pre_mul<nsvec, false>(_mtan_sI, temp_M6, qr5x5_ls);

                     //
                     // dxi_dI has already been folded into dC_quat_dI;
                     // the use of the following here would be incomplete:
                     //
                     // ! get dsigClat_dxi
                     // CALL eval_d_dxi_Slat(dsigClat_dxi, sigC_vecds_lat, dC_matx_dxi, C_matx)
                     // mtan_sI(1:TVEC,:) = mtan_sI(1:TVEC,:) + &
                     // & MATMUL(dsigClat_dxi(1:TVEC,:), dxi_dI(:,:))
                     //
                     // contribution: dSlat_dCmatx(i,p,q) . dCmatx_dCquat(p,q,r) . dCquat_dI(r,j)
                     //
                     double dsigClat_dCquat[ ecmech::ntvec * ecmech::qdim ];
                     {
                        double dCmatx_dCquat[ ecmech::ndim * ecmech::ndim * ecmech::qdim ];
                        d_quat_to_tensor(dCmatx_dCquat, C_quat);

                        double sigC_vecds_lat[ ecmech::nsvec ];
                        _thermoElastN.getCauchy(sigC_vecds_lat, T_vecds, _detV_ri);

                        double dsigClat_dCmatx[ ecmech::ntvec * ecmech::ndim * ecmech::ndim ];
                        d_rot_mat_vecd_smop(dsigClat_dCmatx, C_matx, sigC_vecds_lat);
                        //
                        vecsMAB<ntvec, qdim, ndim*ndim>(dsigClat_dCquat, dsigClat_dCmatx, dCmatx_dCquat);
                     }
                     //
                     double dsigClat_dI[ ecmech::ntvec * nRHS ];
                     vecsMAB<ntvec, nRHS, qdim>(dsigClat_dI, dsigClat_dCquat, dC_quat_dI);
                     //
                     for (int ii_T = 0; ii_T < ecmech::ntvec; ++ii_T) {
                        for (int ii_I = 0; ii_I<nRHS; ++ii_I) {
                           // NOTE : only looping over ntvec, but mtan_sI is nsvec in the first dimension
                           _mtan_sI[ECMECH_NN_INDX(ii_T, ii_I,
                                                   ecmech::nsvec)] += dsigClat_dI[ECMECH_NM_INDX(ii_T, ii_I, ecmech::ntvec, nRHS)];
                        }
                     }
                  } // l_eval_derivs

                  // SCALING
                  {
                     double scaleFactorJ;
                     for (int iJ = 0; iJ < _i_sub_r; ++iJ) {
                        // Jacobian(i_sub_e:i_sup_e,i_sub_e:i_sup_e) = jacob_ee * epsdot_scale_inv  * e_scale ! resid, x
                        scaleFactorJ = _epsdot_scale_inv * ecmech::e_scale;
                        for (int jJ = 0; jJ < _i_sub_r; ++jJ) { // <=_i_sup_e
                           int ijJ = ECMECH_NN_INDX(iJ, jJ, nDimSys);
                           Jacobian[ ijJ ] *= scaleFactorJ;
                        }

                        // Jacobian(i_sub_e:i_sup_e,i_sub_r:i_sup_r) = jacob_er * epsdot_scale_inv  * r_scale
                        scaleFactorJ = _epsdot_scale_inv * ecmech::r_scale * ecmech::zero;
                        for (int jJ = _i_sub_r; jJ< _i_sub_h; ++jJ) { // <_i_sup_r
                           int ijJ = ECMECH_NN_INDX(iJ, jJ, nDimSys);
                           Jacobian[ ijJ ] *= scaleFactorJ;
                        }

                        // Jacobian(i_sub_e:i_sup_e,i_sub_h:i_sup_h) = jacob_eh * epsdot_scale_inv  * h_scale[j]
                        scaleFactorJ = _epsdot_scale_inv;
                        for (int jJ = _i_sub_h; jJ < nDimSys; ++jJ) { // <_i_sup_r
                           int ijJ = ECMECH_NN_INDX(iJ, jJ, nDimSys);
                           Jacobian[ ijJ ] *= scaleFactorJ * _x_scale[jJ - _i_sub_h];
                        }
                     }

                     for (int iJ = _i_sub_r; iJ < _i_sub_h; ++iJ) {
                        // Jacobian(i_sub_r:i_sup_r,i_sub_e:i_sup_e) = jacob_re * rotincr_scale_inv * e_scale
                        scaleFactorJ = _rotincr_scale_inv * ecmech::e_scale * ecmech::zero;
                        for (int jJ = 0; jJ<_i_sub_r; ++jJ) { // <=_i_sup_e
                           int ijJ = ECMECH_NN_INDX(iJ, jJ, nDimSys);
                           Jacobian[ ijJ ] *= scaleFactorJ;
                        }

                        // Jacobian(i_sub_r:i_sup_r,i_sub_r:i_sup_r) = jacob_rr * rotincr_scale_inv * r_scale
                        scaleFactorJ = _rotincr_scale_inv * ecmech::r_scale;
                        for (int jJ = _i_sub_r; jJ < _i_sub_h; ++jJ) { // <_i_sup_r
                           int ijJ = ECMECH_NN_INDX(iJ, jJ, nDimSys);
                           Jacobian[ ijJ ] *= scaleFactorJ;
                        }

                        // Jacobian(i_sub_r:i_sup_r,i_sub_h:i_sup_h) = jacob_rh * rotincr_scale_inv * h_scale[j]
                        scaleFactorJ = _rotincr_scale_inv;
                        for (int jJ = _i_sub_h; jJ < nDimSys; ++jJ) { // <_i_sup_r
                           int ijJ = ECMECH_NN_INDX(iJ, jJ, nDimSys);
                           Jacobian[ ijJ ] *= scaleFactorJ * _x_scale[jJ - _i_sub_h];
                        }
                     }

                     for (int iJ = _i_sub_h; iJ < nDimSys; ++iJ) {
                        // Jacobian(i_sub_h:i_sup_h,i_sub_e:i_sup_e) = jacob_he * h_res_scale[i] * e_scale
                        scaleFactorJ = _res_scale[iJ - _i_sub_h] * ecmech::e_scale;
                        for (int jJ = 0; jJ < _i_sub_r; ++jJ) { // <=_i_sup_e
                           int ijJ = ECMECH_NN_INDX(iJ, jJ, nDimSys);
                           Jacobian[ ijJ ] *= scaleFactorJ;
                        }
                        // Rotation term is 0
                        // Hardening term already scaled
                     }
                  } // SCALING
               } // doComputeJ

               return true;
            } // computeRJ

            __ecmech_hdev__
            inline
            const double* getGdot() const { return _gdot; }

         private:

            const SlipGeom &_slipGeom;
            const Kinetics &_kinetics;
            const ThermoElastN &_thermoElastN;

            double _dt, _detV, _eVref, _p_EOS, _tK, _a_V;
            double _dt_ri, _a_V_ri, _detV_ri;

            double _hdn_scale;
            double _epsdot_scale_inv, _rotincr_scale_inv;

            double _gdot[SlipGeom::nslip]; // crys%tmp1_slp

            double _kin_vals[Kinetics::nVals];
            double _x_scale[Kinetics::nH];
            double _res_scale[Kinetics::nH];

            const double* const _h_state_n;
            const double* const _e_vecd_n;
            const double* const _Cn_quat;
            const double* const _d_vecd_sm; // d_vecds_sm would be fine too -- but do not use _d_vecd_sm[iSvecS];
            const double* const _w_veccp_sm;

            static const int _nXnDim = nDimSys * nDimSys;
            static const int _i_sub_e = 0; // ntvec
            static const int _i_sub_r = ecmech::ntvec; // nwvec
            static const int _i_sub_h = _i_sub_r + ecmech::emapdim; // kinetics.nH

            double _dp_dis_rate_contrib, _shrate_eff_contrib;

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
      bool getResponseFISngl(const SlipGeom& slipGeom,
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

         double Cstr_vecds_lat[ecmech::nsvec];
         //
         double* e_vecd_u = &(hist[iHistLbE]);
         double* quat_u = &(hist[iHistLbQ]);
         double vNew = volRatio[1];
         {
            // get hardness in terms of how the solver wants it
            //
            double h_state_n[Kinetics::nH];
            for (int iH = 0; iH < Kinetics::nH; iH++) {
               h_state_n[iH] = h_state[iH];
            }

            double kin_vals[Kinetics::nVals];
            const double hdn = kinetics.getVals(kin_vals, pEOS, tkelv, h_state_n);
            const double adots_ref = kinetics.getFixedRefRate(kin_vals);

            kinetics.setH0Ext(h_state_n);
            EvptnUpdstFIProblem<SlipGeom, Kinetics, ThermoElastN> prob(slipGeom, kinetics, elastN,
                                                                       dt,
                                                                       vNew, eNew, pEOS, tkelv,
                                                                       h_state_n, e_vecd_n, quat_n,
                                                                       d_vecd_sm, w_veccp_sm,
                                                                       hdn, adots_ref);

            snls::SNLSTrDlDenseG<EvptnUpdstFIProblem<SlipGeom, Kinetics, ThermoElastN> > solver(prob);

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
            // solver.setOutputlevel(10);
            snls::SNLSStatus_t status = solver.solve( );
            //
            if (status != snls::converged) {
#ifdef __cuda_host_only__
               ECMECH_WARN(__func__, "Back-up fully implicit solver failed to converge -- will rerun to get output for debugging");

               // rerun to get more output for debugging
               //
               // get more output
               solver.setOutputlevel(10);
               //
               // reset initial guess
               for (int iX = 0; iX < prob.nDimSys; ++iX) {
                  solver._x[iX] = 0e0;
               }

               //
               // redo solve
               solver.solve( );
#endif
               // False is for the CUDA run so we could catch this and fail if need be
               // after the fact
               return false;
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
            prob.stateFromX(e_vecd_u, quat_u, h_state, solver._x);

            //
            {
               const double* gdot_u = prob.getGdot();
               for (int i_gdot = 0; i_gdot < SlipGeom::nslip; i_gdot++) {
                  gdot[i_gdot] = gdot_u[i_gdot];
               }
            }
            //
            hist[iHistA_shrateEff] = prob.getShrateEff();
            hist[iHistA_shrEff] += hist[iHistA_shrateEff] * dt;
            //
            {
               double dEff = vecd_Deff(d_vecd_sm);
               double flow_strength = prob.getHdnScale();
               if (dEff > idp_tiny_sqrt) {
                  flow_strength = prob.getDisRate() / dEff;
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
      } // getResponseFISngl
   } // namespace evptn
} // namespace ecmech

#endif // ECMECH_EVPTN_FI_H
