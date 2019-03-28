// -*-c++-*-

#ifndef ECMECH_EVPTN_H
#define ECMECH_EVPTN_H

#include <cassert>

#include "ECMech_core.h"
#include "ECMech_util.h"

#include "SNLS_lup_solve.h"

#include "RAJA/RAJA.hpp"

namespace ecmech {

/**
 * for cubic cyrstal symmetry
 *
 * in Fortran mdef coding, corresponds to cem%l_lin_lnsd
 *
 */
class ThermoElastNCubic
{
public:
   static const int nParams = 3 ;
   
   // constructor and destructor
   __ecmech_hdev__
   inline ThermoElastNCubic() {};
   inline ~ThermoElastNCubic() {};
   
   __ecmech_hdev__
   inline void setParams( const real8* const params ) {
      
      int iParam = 0 ;
      real8 c11 = params[iParam++] ;
      real8 c12 = params[iParam++] ;
      real8 c44 = params[iParam++] ;
      assert( iParam == nParams ) ;
      
      _K_diag[0] = c11-c12 ;
      _K_diag[1] = c11-c12 ;
      _K_diag[2] = two*c44 ;
      _K_diag[3] = two*c44 ;
      _K_diag[4] = two*c44 ;
      // _K_vecds_s = c11+two*c12 ; // not used
      
   }
   
   __ecmech_hdev__
   inline void eval( real8* const T_vecds,
                     const real8* const Ee_vecds,
                     real8 , // tK
                     real8 p_EOS,
                     real8 // eVref
                     ) const {

      real8 ln_J    = sqr3 * Ee_vecds[iSvecS] ; // vecds_s_to_trace
      real8 J       = exp(ln_J) ;
      real8 Ts_bulk = -sqr3 * J * p_EOS ;

      vecsVAdiagB<ntvec>( T_vecds, _K_diag, Ee_vecds ) ;
      T_vecds[iSvecS] = Ts_bulk ; // _K_vecds_s * Ee_vecds(SVEC)
   }

   /**
    * dT_deps[0:ntvec-1,:]^T * A, for non-square A[ntvec,q] (with q likely being nSlip)
    * so that even if T_vecds[iSvecS] depends on Ee_vecds, that is not in the result
    *
    * combines calls to elawn_T_dif and eval_dtaua_deps_n
    *
    * for cubic, dT_deps is diag(K_diag * a_V%ri) (symmetric) ; dT_deps[iSvecS,:] = 0
    */
   __ecmech_hdev__
   inline void multDTDepsT( real8* const P, // ntvec*q
                            const real8* const A, // ntvec*q
                            real8 a_V_ri,
                            int q) const {
      for ( int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec ) {
         real8 dTdepsThis = _K_diag[iTvec] * a_V_ri ;
         for ( int iQ=0; iQ < q; ++iQ ) {
            P[ECMECH_NM_INDX(iTvec,iQ,ecmech::ntvec,q)] = dTdepsThis * A[ECMECH_NM_INDX(iTvec,iQ,ecmech::ntvec,q)] ;
         }
      }
   }

   __ecmech_hdev__
   inline void getCauchy( real8* const sigC_vecds_lat,
                          const real8* const T_vecds,
                          real8 detVi ) const
   {
      for ( int iSvec = 0; iSvec < ecmech::nsvec; ++iSvec ) {      
         sigC_vecds_lat[iSvec] = detVi * T_vecds[iSvec] ;
      }
   }
   
   /**
    * like dsigC_de * A, with disgC_de[nsvec,ntvec] having come from elawn_Cauchy_dif
    * for A[ntvec,ntvec]
    *
    * NOTE : dsigC_de[nsvec,ntvec] with nsvec in the first dimension
    * because in general distorational deformation can produce
    * pressure -- for example in materials with hexagonal symmetry,
    * even if it does not happen in cubic symmetry
    */
   __ecmech_hdev__
   inline void multCauchyDif( real8* const M65,
                              const real8* const A,
                              real8 detVi,
                              real8 a_V_ri
                              ) const {
      // CALL vecds_s_to_trace(tr_ln_V, s_meas%Ee_vecds(SVEC))
      // detV = DEXP(tr_ln_V)
      // detVi = one / detV
      
      // dsigC_de(:,:) = detVi * s_meas%dT_deps(:,:)
      // for cubic, dT_deps is diag(K_diag * a_V%ri) (symmetric) ; dT_deps[iSvecS,:] = 0
      // M65_ij = dd_ii A_ij
      for ( int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec ) {
         real8 vFact = detVi * a_V_ri * _K_diag[iTvec] ;
         for ( int jTvec = 0; jTvec < ecmech::ntvec; ++jTvec ) {
            M65[ECMECH_NM_INDX(iTvec,jTvec,ecmech::nsvec,ecmech::ntvec)] = vFact * A[ECMECH_NN_INDX(iTvec,jTvec,ecmech::ntvec)] ;
         }
      }
      for ( int jTvec = 0; jTvec < ecmech::ntvec; ++jTvec ) {
         M65[ECMECH_NM_INDX(iSvecS,jTvec,ecmech::nsvec,ecmech::ntvec)] = 0.0 ;
      }
      
   }
   
private :
   real8 _K_diag[ecmech::ntvec] ;
   // real8 _K_vecds_s ; // not used 
};


template< class SlipGeom, class Kinetics, class ThermoElastN >
class EvptnUpdstProblem
{
public:
   
static const int nDimSys = ecmech::ntvec + ecmech::nwvec ;

// constructor
__ecmech_hdev__
EvptnUpdstProblem(const SlipGeom& slipGeom,
                  const Kinetics& kinetics,
                  const ThermoElastN& thermoElastN,
                  real8 dt, 
                  real8 detV, real8 eVref, real8 p_EOS, real8 tK,
                  const real8* const h_state,
                  const real8* const e_vecd_n,
                  const real8* const Cn_quat,
                  const real8* const d_vecd_sm, // okay to pass d_vecds_sm, but d_vecd_sm[iSvecS] is not used
                  const real8* const w_veccp_sm
                  ) 
   : _slipGeom(slipGeom),
     _kinetics(kinetics),
     _thermoElastN(thermoElastN),
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
   _dt_ri = 1.0 / _dt ;
   _detV_ri = 1.0 / _detV ;
   _a_V = pow(detV, onethird) ;
   _a_V_ri = 1.0 / _a_V ;

   _kinetics.getVals(_kin_vals, _p_EOS, _tK, _h_state) ;
   
   real8 adots_ref = _kinetics.getFixedRefRate(_kin_vals) ;
   real8 eff = vecNorm< ntvec >( _d_vecd_sm ) ; // do not worry about factor of sqrt(twothird)
   if (eff < epsdot_scl_nzeff*adots_ref ) {
      _epsdot_scale_inv = one / adots_ref ;
   } else {
      _epsdot_scale_inv = fmin( one / eff, 1e6 * _dt ) ;
   }
   //
   _rotincr_scale_inv = _dt_ri * _epsdot_scale_inv ;
   
}

__ecmech_hdev__
inline
bool computeRJ( real8* const resid,
                real8* const Jacobian,
                const real8* const x ) {
   bool doComputeJ = (Jacobian != nullptr) ;
      
   if ( doComputeJ ) {

      // zero the Jacobian so that do not need to worry about zero
      // entries in the midst of other things later
      //
      for ( int ijJ=0; ijJ<_nXnDim; ++ijJ ) {
         Jacobian[ijJ] = 0.0;
      }
         
   }
   //
   for ( int iR=0; iR<nDimSys; ++iR ) {
      resid[iR] = 0.0 ;
   }

   //////////////////////////////
   //  PULL VALUES out of x, with scalings
   //
   real8 edot_vecd[ecmech::ntvec] ;
   vecsVxa<ntvec>( edot_vecd, ecmech::e_scale, &(x[_i_sub_e]) ) ; // edot_vecd is now the delta, _not_ yet edot_vecd
   // e_vecd_f is end-of-step
   real8 e_vecd_f[ntvec] ;
   vecsVapb<ntvec>( e_vecd_f, edot_vecd, _e_vecd_n ) ;
   vecsVsa<ntvec>( edot_vecd, _dt_ri ) ; // _now_ edot_vecd has edot_vecd
   //
   real8 xi_f[nwvec] ;
   vecsVxa<nwvec>( xi_f, ecmech::r_scale, &(x[_i_sub_r]) ) ;
   //
   // not done in EvpC :
   // 	CALL exp_map_cpvec(A, xi_f)
   // 	CALL get_c(c, A, C_n)
   //
   real8 A_quat[ecmech::qdim] ;
   emap_to_quat(A_quat, xi_f) ;
   //
   real8 C_quat[ecmech::qdim] ;
   get_c_quat(C_quat, A_quat, _Cn_quat) ;
   //
   real8 C_matx[ecmech::ndim * ecmech::ndim] ;
   quat_to_tensor(C_matx, C_quat) ;
   //
   real8 qr5x5_ls[ecmech::ntvec * ecmech::ntvec] ;
   get_rot_mat_vecd(qr5x5_ls, C_matx) ;
   //
   // CALL matt_x_vec_5(qr5x5_ls, vel_grad_sm%d_vecds(1:TVEC), d_vecd_lat)
   real8 d_vecd_lat[ecmech::ntvec] ;
   vecsVMTa<ntvec>(d_vecd_lat, qr5x5_ls, _d_vecd_sm) ;
   // d_vecds_lat(SVEC) = vel_grad_sm%d_vecds(SVEC)
   //
   // // CALL rot_mat_vecd(A, qr5x5_A)
   // 
   // CALL rot_mat_wveccp(C_matx, qr3x3_ls) // amounts to qr3x3_ls = C_matx
   // CALL matt_x_vec_3(qr3x3_ls, vel_grad_sm%w_veccp, w_vec_lat)
   real8 w_vec_lat[ecmech::nwvec] ; // assumes nwvec = ndim
   vecsVMTa<ndim>(w_vec_lat, C_matx, _w_veccp_sm) ;

   //////////////////////////////
   // CALCULATIONS

   // // do not need to use elaw_T_BT here as T and BT are the same
   //
   // specialize to cem%l_lin_lnsd
   // CALL elawn_T(s_meas, e_vecd_f, crys%elas, tK, .TRUE., a_V, &
   //      & p_EOS, eVref, crys%i_eos_model, crys%eos_const &
   //      &)
   real8 Ee_vecds[ecmech::nsvec];
   vecsVxa<ntvec>(Ee_vecds, _a_V_ri, e_vecd_f) ;
   // // tr_Ee = three * DLOG(a_V%r)
   // // CALL trace_to_vecds_s(s_meas%Ee_vecds(SVEC), tr_Ee)
   Ee_vecds[iSvecS] = sqr3 * log(_a_V) ; // could go into constructor
   //
   // // Kirchhoff stress from Ee_vecds
   //    CALL elawn_lin_op(s_meas%T_vecds, s_meas%Ee_vecds, cem, tK, &
   //         & p_EOS, eVref, i_eos_model, eos_const)
   real8 T_vecds[ecmech::nsvec];
   _thermoElastN.eval(T_vecds, Ee_vecds, _tK, _p_EOS, _eVref) ;
   //
   real8 taua[SlipGeom::nslip] = {0.0} ; // crys%tmp4_slp
   real8 dgdot_dtau[SlipGeom::nslip] = {0.0} ; // crys%tmp2_slp
   real8 dgdot_dg[SlipGeom::nslip] = {0.0} ; // crys%tmp3_slp
   real8 pl_vecd[ecmech::ntvec] = {0.0} ;
   real8 pl_wvec[ecmech::nwvec] = {0.0} ; // \pcDhat
   if ( SlipGeom::nslip > 0 ) {
      
      //  resolve stress onto slip systems
      // CALL resolve_tau_a_n(crys%tmp4_slp, s_meas%T_vecds, crys)
      vecsVaTM< ntvec, SlipGeom::nslip >( taua, T_vecds, _slipGeom.getP() ) ;
      //
      // CALL plaw_eval(pl_vecd, pl_wvec, gss, crys, tK, ierr)
      _kinetics.evalGdots( _gdot, dgdot_dtau, dgdot_dg, taua, _kin_vals ) ;
      //
      // CALL sum_slip_def(pl_vecd, pl_wvec, crys%tmp1_slp, crys) ;
      vecsVMa< ntvec, SlipGeom::nslip >( pl_vecd, _slipGeom.getP(), _gdot ) ;
      vecsVMa< nwvec, SlipGeom::nslip >( pl_wvec, _slipGeom.getQ(), _gdot ) ;
   }
   //
   // // shrate_l%gdot => crys%tmp1_slp

   // from e edot product term in spin (formerly neglected)
   //
   real8 A_e_M35[ecmech::nwvec * ecmech::ntvec] ;
   M35_d_AAoB_dA( A_e_M35, e_vecd_f ) ;
   //
   real8 ee_wvec[ecmech::nwvec] ;
   vecsVMa< nwvec, ntvec >( ee_wvec, A_e_M35, edot_vecd ) ;
   //
   real8 ee_fac = onehalf * _a_V_ri * _a_V_ri ;

   // RESIDUAL B_S
   //
   for ( int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec ) {
      resid[_i_sub_e+iTvec] = _epsdot_scale_inv * ( // SCALING
         _a_V_ri * edot_vecd[iTvec] + pl_vecd[iTvec] - d_vecd_lat[iTvec] ) ;
   }

   // RESIDUAL B_xi
   //
   for ( int iWvec = 0; iWvec < ecmech::nwvec; ++iWvec ) {
      resid[_i_sub_r+iWvec] = _rotincr_scale_inv * // SCALING
         ( xi_f[iWvec] - _dt * (w_vec_lat[iWvec] - pl_wvec[iWvec] + ee_fac * ee_wvec[iWvec]) ) ;
   }

   _dp_dis_rate_contrib = zero ;
   if ( SlipGeom::nslip > 0 ) {
      // CALL calc_pl_dis(dp_dis_rate_contrib, crys%tmp4_slp, crys%tmp1_slp, detV%ri)
      _dp_dis_rate_contrib = _detV_ri * vecsyadotb< SlipGeom::nslip >( taua, _gdot ) ;
   }

   // // need shrate%eff instead
   // // CALL calc_pl_eff(dp_def_rate_contrib, pl_vecd, detV%ri)
   // CALL setup_ss_shrate_vals(shrate_l, crys%tmp1_slp, zero, .TRUE.)
   _shrate_eff_contrib = vecsssumabs<SlipGeom::nslip>(_gdot) ;
    
   //////////////////////////////////////////////////////////////////////
   // JACOBIAN, fixed hardness and temperature
   //
   if ( doComputeJ ) {

      // use RAJA::View machinery to simplify indexing for blocks in the Jacobian matrix ;
      // can always swap this out later if it ends up being too heavyweight ;
      // RAJA defaults to "row-major" -- final dimension indexing the fastest
      // 
      const int JDIM = 2 ;

      // 
      // preliminaries
      // 
      real8 dpl_deps_symm[ ecmech::ntvec * ecmech::ntvec ] = {0.0} ;
      real8 dpl_deps_skew[ ecmech::nwvec * ecmech::ntvec ] = {0.0} ;
      if ( SlipGeom::nslip > 0 ) {
         
         // CALL elawn_T_dif(s_meas, e_vecd_f, crys%elas, tK, a_V, &
         //                  & p_EOS, eVref, crys%i_eos_model, crys%eos_const, &
         //                  & dpEOS_dtK, .FALSE., .FALSE., .FALSE.)
         // CALL eval_dtaua_deps_n(dtaua_deps, s_meas%dT_deps, crys)
         //
         real8 dtaua_deps[ ecmech::ntvec * SlipGeom::nslip ] ;
         _thermoElastN.multDTDepsT( dtaua_deps, _slipGeom.getP(), _a_V_ri, SlipGeom::nslip ) ;
         
         // CALL plaw_eval_dif_sn(TVEC, &
         //                       & dpl_deps_symm, dpl_deps_skew, dgdot_deps, &
         //                       & dtaua_deps, gss, crys, s_meas, .FALSE.)
         real8 dgdot_deps[ ecmech::ntvec * SlipGeom::nslip ] = {0.0} ;
         for ( int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec ) {
            for ( int iSlip=0; iSlip < SlipGeom::nslip; ++iSlip ) {
               int ijThis = ECMECH_NM_INDX(iTvec,iSlip,ecmech::ntvec,SlipGeom::nslip) ;
               dgdot_deps[ijThis] = dgdot_dtau[iSlip] * dtaua_deps[ijThis] ;
            }
         }
         // DO islip = 1, crys%nslip
         //   DO i_TVEC = 1, TVEC
         //     dpl_deps_symm(:,i_TVEC) += 
         //          & crys%P_ref_vec(:,islip) * dgdot_deps(i_TVEC,islip)
         //     dpl_deps_skew(:,i_TVEC) += 
         //          & crys%Q_ref_vec(:,islip) * dgdot_deps(i_TVEC,islip)
         //   END DO
         // END DO
         vecsMABT<        ntvec, SlipGeom::nslip >(dpl_deps_symm, _slipGeom.getP(), dgdot_deps) ;
         vecsMABT< nwvec, ntvec, SlipGeom::nslip >(dpl_deps_skew, _slipGeom.getQ(), dgdot_deps) ;
         
      }
      //
      //
      // derivatives with respect to lattice orientation changes
      real8 dC_quat_dxi_T[ ecmech::nwvec*ecmech::qdim ] ;
      real8 dDsm_dxi[ ecmech::ntvec*ecmech::nwvec ] ;
      real8 dWsm_dxi[ ecmech::nwvec*ecmech::nwvec ] ;
      eval_d_dxi_impl_quat( dC_quat_dxi_T, dDsm_dxi, dWsm_dxi,
                            _d_vecd_sm, _w_veccp_sm, 
                            xi_f, _Cn_quat, C_matx, C_quat ) ;

      // d(B_S)/d(e_vecd_f)
      //
      {
         RAJA::View< real8, RAJA::Layout<JDIM> > jacob_ee(Jacobian, nDimSys, nDimSys) ;
         
         // dislocation plasticity;
         // first contribution; overwrite
         for ( int jTvec = 0; jTvec < ecmech::ntvec; ++jTvec ) {
            for ( int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec ) {
               jacob_ee(iTvec,jTvec) = dpl_deps_symm[ECMECH_NN_INDX(iTvec,jTvec,ecmech::ntvec)] ;
            }
         }

         // elastic rate
         //
         {
            real8 adti = _a_V_ri * _dt_ri ;
            for ( int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec ) {
            jacob_ee(iTvec,iTvec) += adti ;
            }
         }
      }
   
      // d(B_S)/d(xi_f)
      //
      // jacob_er = -dDsm_dxi(:,:)
      {
         RAJA::OffsetLayout<JDIM> layout = RAJA::make_offset_layout<JDIM>({{ 0 , -_i_sub_r }}, {{nDimSys-1, -_i_sub_r+nDimSys-1}});
         RAJA::View<real8, RAJA::OffsetLayout<JDIM> > jacob_er(Jacobian, layout);
         for ( int jWvec=0; jWvec<ecmech::nwvec; ++jWvec) {
            for ( int iTvec=0; iTvec<ecmech::ntvec; ++iTvec) {
               // could also make dDsm_dxi into a RAJA view, but not really needed
               jacob_er(iTvec,jWvec) = -dDsm_dxi[ ECMECH_NM_INDX(iTvec,jWvec,ecmech::ntvec,ecmech::nwvec) ] ;
            }
         }
      }

      // d(B_xi)/d(e_vecds_f)
      //
      {
         RAJA::OffsetLayout<JDIM> layout = RAJA::make_offset_layout<JDIM>({{ -_i_sub_r , 0 }}, {{ -_i_sub_r+nDimSys-1, nDimSys-1 }});
         RAJA::View<real8, RAJA::OffsetLayout<JDIM> > jacob_re(Jacobian, layout);

         real8 A_edot_M35[ecmech::nwvec * ecmech::ntvec] ;
         M35_d_AAoB_dA( A_edot_M35, edot_vecd ) ;

         real8 dt_ee_fac = _dt * ee_fac ;
         
         for ( int iWvec = 0; iWvec < ecmech::nwvec; ++iWvec ) {
            for ( int jTvec = 0; jTvec < ecmech::ntvec; ++jTvec ) {

               int ijWT = ECMECH_NM_INDX(iWvec,jTvec,ecmech::nwvec,ecmech::ntvec) ;
               jacob_re(iWvec,jTvec) =
                  _dt * dpl_deps_skew[ijWT] - dt_ee_fac * ( A_e_M35[ijWT] * _dt_ri - A_edot_M35[ijWT] ) ;
            }
         }
      }
      
      // d(B_xi)/d(xi_f)
      // 
      {
         RAJA::OffsetLayout<JDIM> layout = RAJA::make_offset_layout<JDIM>({{ -_i_sub_r , -_i_sub_r }}, {{ -_i_sub_r+nDimSys-1, -_i_sub_r+nDimSys-1 }});
         RAJA::View<real8, RAJA::OffsetLayout<JDIM> > jacob_rr(Jacobian, layout);

         for ( int iWvec = 0; iWvec < ecmech::nwvec; ++iWvec ) {
            for ( int jWvec = 0; jWvec < ecmech::nwvec; ++jWvec ) {
               int ijWW = ECMECH_NN_INDX(iWvec,jWvec,ecmech::nwvec) ;               
               jacob_rr(iWvec,jWvec) = - _dt * dWsm_dxi[ijWW] ;
            }
            jacob_rr(iWvec,iWvec) += one ;
         }
      }

      if ( _mtan_sI ) { // l_eval_derivs
         // material tangent, do before scaling of Jacobian
         // eval_mtan()
         for (int iiST=0; iiST<ecmech::nsvec*ecmech::ntvec; ++iiST) {
            _mtan_sI[iiST] = 0.0 ;
         }

         // must solve a set of systems to get needed partial derivatives

         // eval_mtan_pfrac_r(de_dI, dxi_dI)
         // compared to Fortran coding, dI has reduced back down to being only the deviatoric part of the deformation rate ; UB_I = ntvec
         const int nRHS = ecmech::ntvec ;
         real8 pfrac_rhs_T[ nRHS*nDimSys ] = {0.0} ; // transpose for use in SNLS_LUP_SolveX !
         // derivatives end up in pfrac_rhs_T
         // de_dI  is pfrac_rhs_T[ :, _i_sub_e:i_sup_e ] // ecmech::ntvec * ecmech::ntvec
         // dxi_dI is pfrac_rhs_T[ :, _i_sub_r:i_sup_r ] // ecmech::nwvec * ecmech::ntvec
         {
          
            // RHS
            //
            {
            
               // negatives cancel
               // pfrac_rhs(i_sub_e:i_sup_e,I_subD:I_supD) = TRANSPOSE(qr5x5_ls)
               for (int jE=0; jE<nRHS; ++jE) {
                  for (int iE=0; iE<nDimSys; ++iE) {
                     // pfrac_rhs[ECMECH_NM_INDX(iE,jE,nDimSys,nRHS)] = qr5x5_ls[ECMECH_NN_INDX(jE,iE,ntvec)] ;
                     pfrac_rhs_T [ECMECH_NM_INDX(jE,iE,nRHS,nDimSys)] = qr5x5_ls[ECMECH_NN_INDX(jE,iE,ntvec)] ;
                  }
               }
            }
            
            // SYSTEM
            //
            real8 pfrac_sys[ _nXnDim ] ;
            for ( int iNXN=0; iNXN<_nXnDim; ++iNXN ) {
               pfrac_sys[iNXN] = Jacobian[iNXN] ;
            }
            int err = SNLS_LUP_SolveX(pfrac_sys, pfrac_rhs_T, nDimSys, nRHS) ;
            if ( err != 0 ) {
               ECMECH_FAIL(__func__,"error from SNLS_LUP_SolveX") ;
            }

         } // eval_mtan_pfrac_r

         // in Fortran code dCn_quat_dI was used to store dC_quat_dI, but here just call it dC_quat_dI
         real8 dC_quat_dI[ ecmech::qdim * nRHS ] ;
         for ( int ii_I=0; ii_I<nRHS; ++ii_I ) {
            for ( int ii_Q=0; ii_Q<ecmech::qdim; ++ii_Q) {
               int iiQI = ECMECH_NM_INDX(ii_Q,ii_I,ecmech::qdim,nRHS) ;
               dC_quat_dI[iiQI] = 0.0 ;
               for ( int ii_W=0; ii_W<ecmech::nwvec; ++ii_W ) {
                  // dC_quat_dI[iiQI] += dC_quat_dxi_T(ii_W,ii_Q) * dxi_dI(ii_W,ii_I)
                  dC_quat_dI[iiQI] +=
                     dC_quat_dxi_T[ECMECH_NM_INDX(ii_W,ii_Q,ecmech::nwvec,ecmech::qdim)] *
                     pfrac_rhs_T[ECMECH_NM_INDX(ii_I,_i_sub_r+ii_W,nRHS,nDimSys)] ;
               }
            }
         }

         
         // contribution through e
         //
         real8 temp_M6I[ ecmech::nsvec*nRHS ] ; // (SVEC,UB_I)
         // dsigClat_def(:,:) = detVi * s_meas%dT_deps(:,:)
         // temp_M6I = MATMUL(dsigClat_def(:,1:TVEC), de_dI(:,:))
         {
            // TODO : get rid of the need for this memory copy (with transpose) into de_dI
            real8 de_dI[ ecmech::ntvec*nRHS ]; // nRHS=ecmech::ntvec, but put nRHS in here for clarity, and thus use ECMECH_NM_INDX instead of ECMECH_NN_INDX
            for ( int iTvec=0; iTvec<ecmech::ntvec; ++iTvec) {
               for ( int jTvec=0; jTvec<nRHS; ++jTvec) {
                  de_dI[ECMECH_NM_INDX(iTvec,jTvec,ecmech::ntvec,nRHS)] = pfrac_rhs_T[ECMECH_NM_INDX(jTvec,iTvec,nRHS,nDimSys)] ;
               }
            }
            _thermoElastN.multCauchyDif(temp_M6I, de_dI, _detV_ri, _a_V_ri )  ;
         }
         //
         // CALL qr6x6_pre_mul(mtan_sI, temp_M6I, qr5x5_ls, UB_I, .FALSE.)
         qr6x6_pre_mul<nRHS,false>(_mtan_sI, temp_M6I, qr5x5_ls) ; // UB_I=nRHS
         // SUBROUTINE qr6x6_pre_mul(M_out, M_in, qr5x5, n, l_T)
         
         //
         // dxi_dI has already been folded into dC_quat_dI;
         // the use of the following here would be incomplete:
         //
         //	! get dsigClat_dxi
         //	CALL eval_d_dxi_Slat(dsigClat_dxi, sigC_vecds_lat, dC_matx_dxi, C_matx)
         //	mtan_sI(1:TVEC,:) = mtan_sI(1:TVEC,:) + &
         //	     & MATMUL(dsigClat_dxi(1:TVEC,:), dxi_dI(:,:))
         //
         // contribution: dSlat_dCmatx(i,p,q) . dCmatx_dCquat(p,q,r) . dCquat_dI(r,j)
         //
         real8 dsigClat_dCquat[ ecmech::ntvec * ecmech::qdim ] ;
         {
            real8 dCmatx_dCquat[ ecmech::ndim * ecmech::ndim * ecmech::qdim ] ;
            d_quat_to_tensor(dCmatx_dCquat, C_quat) ;

            real8 sigC_vecds_lat[ ecmech::nsvec ] ;
            _thermoElastN.getCauchy( sigC_vecds_lat, T_vecds, _detV_ri ) ;
         
            real8 dsigClat_dCmatx[ ecmech::ntvec * ecmech::ndim * ecmech::ndim ] ;
            d_rot_mat_vecd_smop(dsigClat_dCmatx, C_matx, sigC_vecds_lat) ;
            //
            vecsMAB< ntvec, qdim, ndim*ndim >( dsigClat_dCquat, dsigClat_dCmatx, dCmatx_dCquat ) ;
         }
         //
         real8 dsigClat_dI[ ecmech::ntvec * nRHS ] ;
         vecsMAB< ntvec, nRHS, qdim >( dsigClat_dI, dsigClat_dCquat, dC_quat_dI ) ;
         //
         for ( int ii_T = 0; ii_T < ecmech::ntvec; ++ii_T ) {            
            for ( int ii_I=0; ii_I<nRHS; ++ii_I ) {
               // NOTE : only looping over ntvec, but mtan_sI is nsvec in the first dimension
               _mtan_sI[ECMECH_NM_INDX(ii_T,ii_I,ecmech::nsvec,nRHS)] += dsigClat_dI[ECMECH_NM_INDX(ii_T,ii_I,ecmech::ntvec,nRHS)] ;
            }
         }

      } // l_eval_derivs

      // SCALING
      {
         real8 scaleFactorJ ; 
         for ( int iJ=0; iJ<_i_sub_r; ++iJ) {
            
            // Jacobian(i_sub_e:i_sup_e,i_sub_e:i_sup_e) = jacob_ee * epsdot_scale_inv  * e_scale ! resid, x
            scaleFactorJ = _epsdot_scale_inv  * ecmech::e_scale ;
            for ( int jJ=0; jJ<_i_sub_r; ++jJ) { // <=_i_sup_e
               int ijJ = ECMECH_NN_INDX(iJ,jJ,nDimSys) ;
               Jacobian[ ijJ ] *= scaleFactorJ ;
            }

            // Jacobian(i_sub_e:i_sup_e,i_sub_r:i_sup_r) = jacob_er * epsdot_scale_inv  * r_scale           
            scaleFactorJ = _epsdot_scale_inv  * ecmech::r_scale ;
            for ( int jJ=_i_sub_r; jJ<nDimSys; ++jJ) { // <_i_sup_r
               int ijJ = ECMECH_NN_INDX(iJ,jJ,nDimSys) ;
               Jacobian[ ijJ ] *= scaleFactorJ ;
            }
            
         }
         for ( int iJ=_i_sub_r; iJ<nDimSys; ++iJ) {

            // Jacobian(i_sub_r:i_sup_r,i_sub_e:i_sup_e) = jacob_re * rotincr_scale_inv * e_scale           
            scaleFactorJ = _rotincr_scale_inv * ecmech::e_scale ;
            for ( int jJ=0; jJ<_i_sub_r; ++jJ) { // <=_i_sup_e
               int ijJ = ECMECH_NN_INDX(iJ,jJ,nDimSys) ;
               Jacobian[ ijJ ] *= scaleFactorJ ;
            }

            // Jacobian(i_sub_r:i_sup_r,i_sub_r:i_sup_r) = jacob_rr * rotincr_scale_inv * r_scale           
            scaleFactorJ = _rotincr_scale_inv * ecmech::r_scale ;
            for ( int jJ=_i_sub_r; jJ<nDimSys; ++jJ) { // <_i_sup_r
               int ijJ = ECMECH_NN_INDX(iJ,jJ,nDimSys) ;
               Jacobian[ ijJ ] *= scaleFactorJ ;
            }
            
         }
         
      } // SCALING
  
   } // doComputeJ

   return true ;
   
} // computeRJ

__ecmech_hdev__
inline
const real8* getGdot() const { return _gdot; };
   
private:
   
   const SlipGeom &_slipGeom ;
   const Kinetics &_kinetics ;
   const ThermoElastN &_thermoElastN ;

   real8 _dt, _detV, _eVref, _p_EOS, _tK, _a_V ;
   real8 _dt_ri, _a_V_ri, _detV_ri ;

   real8 _epsdot_scale_inv, _rotincr_scale_inv ;

   real8 _gdot[SlipGeom::nslip] ; // crys%tmp1_slp

   real8 _kin_vals[Kinetics::nVals] ;

   const real8* const _h_state ;
   const real8* const _e_vecd_n ;
   const real8* const _Cn_quat ;
   const real8* const _d_vecd_sm ; // d_vecds_sm would be fine too -- but do not use _d_vecd_sm[iSvecS];
   const real8* const _w_veccp_sm ;
   
   static const int _nXnDim = nDimSys * nDimSys ;
   static const int _i_sub_e = 0; // ntvec
   static const int _i_sub_r = ecmech::ntvec ; // nwvec

   real8 _dp_dis_rate_contrib, _shrate_eff_contrib ;

   // for mtan (material tangent stiffnes)
   real8* _mtan_sI ; // null if not wanting tangent evaluation

}; // class EvptnUpdstProblem

} // namespace ecmech

#endif  // ECMECH_EVPTN_H
