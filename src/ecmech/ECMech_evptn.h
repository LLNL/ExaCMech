// -*-c++-*-

#ifndef ECMECH_EVPTN_H
#define ECMECH_EVPTN_H

#include <cassert>

#include "ECMech_core.h"
#include "ECMech_util.h"
#include "ECMech_eosSimple.h"

#include "SNLS_lup_solve.h"
#include "SNLS_TrDLDenseG.h"

#include "RAJA/RAJA.hpp"

namespace ecmech {

namespace evptn {
   
const int numHistAux = 3 ; // effective shearing rate, accumulated shear, nFEval
//
const int iHistLbA = 0 ;
const int iHistA_shrateEff = iHistLbA+0 ;
const int iHistA_shrEff    = iHistLbA+1 ;
const int iHistA_nFEval    = iHistLbA+2 ;
const int iHistLbE = numHistAux ;
const int iHistLbQ = numHistAux + ecmech::ntvec ;
const int iHistLbH = numHistAux + ecmech::ntvec + ecmech::qdim ; 
   
/*
 * just a container for a traits
 */
template< class SlipGeom, class Kinetics, class ThermoElastN, class EosModel >
class NumHist
{
public:
   // see n_rsv_matmod in F90 code
   static const int iHistLbGdot = iHistLbH + Kinetics::nH ;   
   static const int numHist = iHistLbH + Kinetics::nH + SlipGeom::nslip ;
}; // NumHist

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
   inline ThermoElastNCubic() : _K_bulkMod(-1.0), _K_gmod(-1.0) {};
   inline ~ThermoElastNCubic() {};
   
   __ecmech_hdev__
   inline void setParams( const std::vector<real8> & params // const real8* const params
                          ) {
      
      std::vector<real8>::const_iterator parsIt = params.begin();
      
      _c11 = *parsIt; ++parsIt;
      _c12 = *parsIt; ++parsIt;
      _c44 = *parsIt; ++parsIt;
      //
      int iParam = parsIt - params.begin();
      assert( iParam == nParams );
      
      _K_diag[0] = _c11 - _c12 ;
      _K_diag[1] = _c11 - _c12 ;
      _K_diag[2] = two*_c44 ;
      _K_diag[3] = two*_c44 ;
      _K_diag[4] = two*_c44 ;
      real8 K_vecds_s = _c11 + two*_c12 ;
      _K_bulkMod = onethird * K_vecds_s ;
      _K_gmod = (two*_c11-two*_c12+six*_c44)*0.2 ; // average of _K_diag entries
      
   }
   
   __ecmech_hdev__
   inline void getParams( std::vector<real8> & params
                          ) const {
      
      // do not clear params in case adding to an existing set
      int paramsStart = params.size() ;
      
      params.push_back(_c11) ;
      params.push_back(_c12) ;
      params.push_back(_c44) ;

      int iParam = params.size() - paramsStart;
      assert( iParam == nParams );
      
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
    * dT_deps[0:ntvec-1,:]^T * A, for non-square A[ntvec,p] (with p likely being nSlip)
    * so that even if T_vecds[iSvecS] depends on Ee_vecds, that is not in the result
    *
    * combines calls to elawn_T_dif and eval_dtaua_deps_n
    *
    * for cubic, dT_deps is diag(K_diag * a_V%ri) (symmetric) ; dT_deps[iSvecS,:] = 0
    */
   __ecmech_hdev__
   inline void multDTDepsT( real8* const P, // ntvec*p
                            const real8* const A, // ntvec*p
                            real8 a_V_ri,
                            int p) const {
      for ( int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec ) {
         real8 dTdepsThis = _K_diag[iTvec] * a_V_ri ;
         for ( int iP=0; iP < p; ++iP ) {
            P[ECMECH_NM_INDX(iTvec,iP,ecmech::ntvec,p)] = dTdepsThis * A[ECMECH_NM_INDX(iTvec,iP,ecmech::ntvec,p)] ;
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
    *
    * NOTE : M6[nsvec,nsvec] with nsvec in the second dimension
    * (instead of ntvec) to make things easier elsewhere
    */
   __ecmech_hdev__
   inline void multCauchyDif( real8* const M6,
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
            M6[ECMECH_NN_INDX(iTvec,jTvec,ecmech::nsvec)] = vFact * A[ECMECH_NN_INDX(iTvec,jTvec,ecmech::ntvec)] ;
         }
      }
      for ( int jTvec = 0; jTvec < ecmech::ntvec; ++jTvec ) {
         M6[ECMECH_NN_INDX(iSvecS,jTvec,ecmech::nsvec)] = 0.0 ;
      }
      for ( int iSvec = 0; iSvec < ecmech::nsvec; ++iSvec ) {
         M6[ECMECH_NN_INDX(iSvec,iSvecS,ecmech::nsvec)] = 0.0 ;
      }
      
   }
   
   __ecmech_hdev__
   inline real8 getBulkMod( ) const {
      if ( _K_bulkMod <= 0.0 ) {
         ECMECH_FAIL(__func__,"bulk modulus negative -- not initialized?") ;
      }
      return _K_bulkMod ;
   }
  
   __ecmech_hdev__
   inline real8 getGmod( real8 , // tK
                         real8 , // p_EOS
                         real8   // eVref
                         ) const {
      if ( _K_gmod <= 0.0 ) {
         ECMECH_FAIL(__func__,"effective shear modulus negative -- not initialized?") ;
      }
      return _K_gmod ;
   }
  
private :
   real8 _c11, _c12, _c44 ;
   real8 _K_diag[ecmech::ntvec] ;
   real8 _K_bulkMod, _K_gmod ;
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
void provideMTan(real8* mtan_sI ) { _mtan_sI = mtan_sI ; }
   
__ecmech_hdev__
inline
void clearMTan( ) { _mtan_sI = nullptr ; }
   
__ecmech_hdev__
inline
real8 getDtRi() const { return _dt_ri ; }
   
__ecmech_hdev__
inline
real8 getShrateEff() const { return _shrate_eff_contrib ; }
   
__ecmech_hdev__
inline
real8 getDisRate() const { return _dp_dis_rate_contrib ; }
   
/*
 * NOTES :
 *	() should be equivalent to what happens in computeRJ
 *	() not necessarily safe if e_vecd is the same memory as _e_vecd_n or quat is the same as _Cn_quat
 */
__ecmech_hdev__
inline
void stateFromX( real8* const e_vecd,
                 real8* const quat,
                 const real8* const x ) {

   real8 e_vecd_delta[ecmech::ntvec] ;
   vecsVxa<ntvec>( e_vecd_delta, ecmech::e_scale, &(x[_i_sub_e]) ) ;
   vecsVapb<ntvec>( e_vecd, e_vecd_delta, _e_vecd_n ) ;

   real8 xi_f[nwvec] ;
   vecsVxa<nwvec>( xi_f, ecmech::r_scale, &(x[_i_sub_r]) ) ;
   //
   real8 A_quat[ecmech::qdim] ;
   emap_to_quat(A_quat, xi_f) ;
   //
   // real8 C_quat[ecmech::qdim] ;
   // get_c_quat(C_quat, A_quat, _Cn_quat) ;
   // std::copy(C_quat, C_quat+ecmech:qdim, quat) ;
   get_c_quat(quat, A_quat, _Cn_quat) ;   
}

__ecmech_hdev__
inline
void elastNEtoT( real8* const T_vecds, // nsvec
                 const real8* const e_vecd_f // ntvec
                 ) {

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
   _thermoElastN.eval(T_vecds, Ee_vecds, _tK, _p_EOS, _eVref) ;

}

__ecmech_hdev__
inline
void elastNEtoC( real8* const C_vecds, // nsvec
                 const real8* const e_vecd_f // ntvec
                 ) {
   real8 T_vecds[ecmech::nsvec];
   this->elastNEtoT( T_vecds, e_vecd_f ) ;
   _thermoElastN.getCauchy( C_vecds, T_vecds, _detV_ri ) ;

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

   real8 T_vecds[ecmech::nsvec];
   this->elastNEtoT( T_vecds, e_vecd_f ) ;
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
         //
         // material tangent, do before scaling of Jacobian
         // eval_mtan()

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
                  for (int iE=0; iE<ntvec; ++iE) { // ntvec, _not_ nDimSys // iE is same as index into nDimSys, give how d(resid)/d(_d_vecd_sm) works out
                     // pfrac_rhs[ECMECH_NM_INDX(iE,jE,nDimSys,nRHS)] = qr5x5_ls[ECMECH_NN_INDX(jE,iE,ntvec)] ;
                     pfrac_rhs_T [ECMECH_NM_INDX(jE,iE,nRHS,nDimSys)] = qr5x5_ls[ECMECH_NN_INDX(jE,iE,ntvec)] ;
                  }
               }
            }
            
            // SYSTEM
            //
            real8 pfrac_sys[ _nXnDim ] ;
            std::copy(Jacobian, Jacobian+_nXnDim, pfrac_sys) ;
            int err = SNLS_LUP_SolveX<nDimSys>(pfrac_sys, pfrac_rhs_T, nRHS) ;
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
         // real8 temp_M6I[ ecmech::nsvec*nRHS ] ; // (SVEC,UB_I)
         real8 temp_M6[ ecmech::nsvec2 ] ; // (SVEC,UB_I)
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
            _thermoElastN.multCauchyDif(temp_M6, de_dI, _detV_ri, _a_V_ri )  ;
         }
         //
         // CALL qr6x6_pre_mul(mtan_sI, temp_M6I, qr5x5_ls, UB_I, .FALSE.)
         // UB_I=nRHS ; but here do nsvec instead of nRHS so that there is less monkeying with memory later
         qr6x6_pre_mul<nsvec,false>(_mtan_sI, temp_M6, qr5x5_ls) ; 
         
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
               _mtan_sI[ECMECH_NN_INDX(ii_T,ii_I,ecmech::nsvec)] += dsigClat_dI[ECMECH_NM_INDX(ii_T,ii_I,ecmech::ntvec,nRHS)] ;
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

/*
 * for steady-flow capability, might want to check out Dlsmm_getEnabled() stuff in EvpC.c
 *
 * convention for spin coming in should be consistent with w_veccp_sm convention
 */
template< class SlipGeom, class Kinetics, class ThermoElastN, class EosModel >
__ecmech_hdev__
inline
void getResponseSngl(const SlipGeom& slipGeom,
                     const Kinetics& kinetics,
                     const ThermoElastN& elastN,
                     const EosModel& eos,
                           real8    dt,
                           real8    tolerance, 
                     const real8  * d_svec_kk_sm,  // defRate,
                     const real8  * w_veccp_sm, // spin
                     const real8  * volRatio,
                           real8  * eInt,
                           real8  * stressSvecP,
                           real8  * hist,
                           real8  & tkelv,
                           real8  * sdd,
                           real8  * mtanSD,
                     
                     int outputLevel = 0 ) 
{

   static const int iHistLbGdot = NumHist<SlipGeom,Kinetics,ThermoElastN,EosModel>::iHistLbGdot ;

   // NOTE : mtanSD can be nullptr
   //
   bool haveMtan = ( mtanSD != nullptr ) ;

   // convert deformation rate convention
   //
   real8 d_vecd_sm[ecmech::ntvec] ;
   svecToVecd( d_vecd_sm, d_svec_kk_sm ) ;
#ifdef NEED_SCALAR_FLOW_STRENGTH
   real8 dEff = vecd_Deff( d_vecd_sm ) ;
#endif
   
   // pointers to state
   //
   real8* h_state = &(hist[iHistLbH]) ;
   real8* gdot    = &(hist[iHistLbGdot]) ;
   //
   // copies, to keep beginning-of-step state safe
   //
   real8 e_vecd_n[ecmech::ntvec] ;
   std::copy(hist+iHistLbE, hist+iHistLbE+ecmech::ntvec, e_vecd_n) ;
   real8 quat_n[ecmech::qdim] ;
   std::copy(hist+iHistLbQ, hist+iHistLbQ+ecmech::qdim,  quat_n) ;
   //
   // normalize quat just in case
   vecsVNormalize<qdim>(quat_n) ;

   // total increment in the deviatoric part of the strain energy using
   // trapezoidal rule integration
   //
   // just beginning-of-step stress part so far
   //
   real8 halfVMidDt = oneqrtr * (volRatio[0]+volRatio[1]) * dt ;
   real8 eDevTot = halfVMidDt * vecsInnerSvecDev( stressSvecP, d_svec_kk_sm ) ;
   
   // EOS
   //
   real8 eOld = eInt[0] ;
   real8 pOld = stressSvecP[6] ;
   real8 pEOS, eNew, bulkNew ;
   //
   // get tkelv from beginning-of-step to avoid tangent stiffness contributions
   {
      real8 pBOS ;
      real8 vOld = volRatio[0] ;
      eos.evalPT(pBOS, tkelv, vOld, eOld);
   }
   //
   {
      real8 tkelvNew ;
      updateSimple<EosModel>( eos, pEOS, tkelvNew, eNew, bulkNew,
                              volRatio, eOld, pOld ) ;
   }

   // update hardness state to the end of the step
   // gdot is still at beginning-of-step
   //
   real8 h_state_u[Kinetics::nH] ;
   kinetics.updateH( h_state_u, h_state, dt, gdot ) ;

   real8 Cstr_vecds_lat[ecmech::nsvec] ;
   //
   real8* e_vecd_u  = &(hist[iHistLbE]) ;
   real8* quat_u    = &(hist[iHistLbQ]) ;
   {
      real8 vNew = volRatio[1] ;
      EvptnUpdstProblem< SlipGeom, Kinetics, ThermoElastN > prob(slipGeom, kinetics, elastN,
                                                                 dt,
                                                                 vNew, eNew, pEOS, tkelv, 
                                                                 h_state_u, e_vecd_n, quat_n,
                                                                 d_vecd_sm, w_veccp_sm ) ;
   
      snls::SNLSTrDlDenseG< EvptnUpdstProblem<SlipGeom, Kinetics, ThermoElastN> > solver(prob) ;

      snls::TrDeltaControl deltaControl ;
      deltaControl._deltaInit = 1e0 ;
      {
         static const int maxIter = 100 ;
         solver.setupSolver(maxIter, tolerance, &deltaControl, outputLevel) ;
      }

      // set initial guess
      //
      real8* x = solver.getXPntr() ;
      for (int iX = 0; iX < prob.nDimSys; ++iX) {
         x[iX] = 0e0 ;
      }
   
      snls::SNLSStatus_t status = solver.solve( ) ;
      //
      if ( status != snls::converged ){
         
#ifdef __cuda_host_only__
         ECMECH_WARN(__func__,"Solver failed to converge -- will rerun to get output for debuggin");
         
         // rerun to get more output for debugging
         //
         // get more output
         solver.setOutputlevel( 10 ) ;
         //
         // reset initial guess
         for (int iX = 0; iX < prob.nDimSys; ++iX) { x[iX] = 0e0 ; }
         //
         // redo solve
         solver.solve( ) ;
#endif
         
         ECMECH_FAIL(__func__,"Solver failed to converge!");
      }
      // std::cout << "Function evaluations: " << solver.getNFEvals() << std::endl ;

      if ( haveMtan ) {

         real8 mtanSD_vecds[ ecmech::nsvec2 ] ;
         prob.provideMTan( mtanSD_vecds ) ; 
         solver.computeRJ() ;
         prob.clearMTan() ;

         // currently have derivative with-respsect-to deformation rate ;
         // to get derivative with-respsect-to strain increment,
         // multiply by 1/dt 
         //
         real8 dt_ri = prob.getDtRi() ;
         for ( int i=0; i<ecmech::nsvec2; ++i ) {
            mtanSD_vecds[i] = mtanSD_vecds[i] * dt_ri ;
         }
         
         // contribution to stiffness from EOS 
         // this is a bit crude, but should do the trick for now;
         // neglects effect of pEOS and vNew on workings of evptn
         // 
         mtanSD_vecds[ECMECH_NN_INDX(iSvecS,iSvecS,ecmech::nsvec)] = three * bulkNew ;

         // convert from vecds notation to svec notation
         //
         mtan_conv_sd_svec<true>(mtanSD, mtanSD_vecds) ;
         
      }

      // store updated state
      //
      prob.stateFromX(e_vecd_u, quat_u, x) ;
      std::copy(h_state_u, h_state_u+Kinetics::nH, h_state) ;
      {
         const real8* gdot_u = prob.getGdot() ;
         std::copy(gdot_u, gdot_u+SlipGeom::nslip, gdot) ;
      }
      hist[iHistA_shrateEff]  =  prob.getShrateEff() ;
      hist[iHistA_shrEff]    +=  hist[iHistLbA+0] * dt ;
      hist[iHistA_nFEval]     =  solver.getNFEvals() ; // does _not_ include updateH iterations

      //  get Cauchy stress
      //
      prob.elastNEtoC( Cstr_vecds_lat, e_vecd_u ) ;
      
#ifdef NEED_SCALAR_FLOW_STRENGTH
      real8 flow_strength = 0.0 ; // TO_DO : consider setting this instead to a reference value (dependent on the current state?)?
      if ( dEff > idp_tiny_sqrt ) {
         flow_strength = prob.getDisRate() / dEff ;
      }
#endif
      
   }

   real8 C_matx[ecmech::ndim * ecmech::ndim] ;
   quat_to_tensor( C_matx, quat_u ) ;   
   //
   real8 qr5x5_ls[ecmech::ntvec * ecmech::ntvec] ;
   get_rot_mat_vecd(qr5x5_ls, C_matx) ;
   //
   real8 Cstr_vecds_sm[ecmech::nsvec] ;
   vecsVMa<ntvec>(Cstr_vecds_sm, qr5x5_ls, Cstr_vecds_lat) ;
   Cstr_vecds_sm[iSvecS] = Cstr_vecds_lat[iSvecS] ;
   //
   // put end-of-step stress in stressSvecP
   vecdsToSvecP( stressSvecP, Cstr_vecds_sm) ;
   //
   // and now the second half of the trapezoidal integration
   //
   eDevTot += halfVMidDt * vecsInnerSvecDev( stressSvecP, d_svec_kk_sm ) ;

   // adjust sign on quat so that as close as possible to quat_o;
   // more likely to keep orientations clustered this way
   //
   if ( vecsyadotb<qdim>(quat_u, quat_n) < zero ) {
      for ( int iQ=0; iQ<ecmech::qdim; ++iQ ) { quat_u[iQ] = -quat_u[iQ] ; }
   }

   {
      real8 gmod = elastN.getGmod( tkelv, pEOS, eNew ) ;
      sdd[i_sdd_bulk] = bulkNew ;
      sdd[i_sdd_gmod] = gmod ;
   }
#ifdef DEBUG
   assert(ecmech::nsdd == 2) ;
#endif

   eNew = eNew + eDevTot; 
   //
   // could update pressure and temperature again, but do not bother

   eInt[ecmech::i_ne_total] = eNew ;
#ifdef DEBUG
   assert(ecmech::ne == 1) ;
#endif

} // getResponseSngl
   
} // namespace evptn

} // namespace ecmech

#endif  // ECMECH_EVPTN_H
