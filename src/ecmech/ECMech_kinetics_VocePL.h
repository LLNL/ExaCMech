// -*-c++-*-

#ifndef ECMECH_KINETICS_VOCEPL_H
#define ECMECH_KINETICS_VOCEPL_H

#include <cassert>
#include <cmath>

#include "SNLS_TrDLDenseG.h"

#include "ECMech_core.h"
#include "ECMech_util.h"

namespace ecmech {

/**
 * slip and hardening kinetics
 *
 * power-law slip kinetics with Voce hardening law -- meant to be about as simple as it gets
 *
 * relevant stuff from F90 coding :

  !    params: h0, tausi, taus0, xms, gamss0
  !    vals: sv_sat

!!!!!!!!!! hl_values_set

      ! saturation sv depends on shrate
      IF (shrate%eff > idp_tiny_sqrt) THEN
        hl_params%vals(1) = hl_params%params(3) * &
             & ((shrate%eff / hl_params%params(5))**hl_params%params(4))
      ELSE
        hl_params%vals(1) = hl_params%params(3)
      ENDIF
      hl_params%nl_control%res_scale = sv(1)

!!!!!!!!!! hard_law_ms
! ! func is sdot such that :
! sv(:) = svo(:) + dtime * sdot(:)

       temp2 = hl_params%vals(1) - hl_params%params(2)

       IF (PRESENT(dfdtK)) THEN
         dfdtK(1) = zero
       END IF

       IF (temp2 .LE. zero) THEN
          ! most likely, shrate is small enough that hs_sat has become
          ! small, but for small shrate have small rate of hs change
          IF (PRESENT(func)) func(1) = zero
          IF (PRESENT(dfunc)) dfunc(1,1) = zero
          IF (PRESENT(dfdshr)) dfdshr(1,1) = zero
          RETURN
       ENDIF

       IF (PRESENT(func) .OR. PRESENT(dfdshr)) THEN
         temp1 = hl_params%params(1) * ((hl_params%vals(1) - hs(1)) / &
            &        temp2)
         IF (PRESENT(func)) func(1) = temp1 * shrate%eff
         IF (PRESENT(dfdshr)) THEN
           dfdshr(1,1) = temp1 + &
                & hl_params%params(1) * ( &
                &   (hs(1) - hl_params%params(2)) / &
                &      (temp2*temp2)) * &
                & hl_params%params(4) * hl_params%vals(1)
         END IF
       END IF

       IF (PRESENT(dfunc)) &
            & dfunc(1,1) = - hl_params%params(1) / temp2 * shrate%eff

!!!!!!!!!! hs_to_hdn_scalar

      hdn_scalar = h_state(1)

!!!!!!!!!! hs_to_gss

      temp   = h_state(1)
      temp_b = zero

!!!!!!!!!! read_hard_law

       hl_params%nparams = 5
       nparams_get = 5
       hl_params%nvals = 1
       hl_params%nshr_eqv_nhdn = .true.
       hl_params%l_fast_kinetics = .FALSE.

       hl_params%scale_shr_hdn_rate = hl_params%params(1)
       hl_params%scale_strength = hl_params%hdn_init(1)
       hl_params%scale_hs(:) = hl_params%params(3) ! hl_params%hdn_init(:)
       hl_params%scale_hs_inv(:) = one / hl_params%scale_hs(:)
#ifdef USE_UPDSV
       hl_params%nl_control%res_scale = hl_params%hdn_init(1)
       CALL set_step_ctrl_nl(hl_params%nl_control, hl_params%hdn_init(1))

!!!!!!!!!! fix_h

      IF (hs(1) <= hs_min_p) THEN
        IF (l_m_small) THEN
          hs(1) = hl_params%hdn_init(1)
        ELSE
          CALL consider_ierr(1,location,CIERR_GEN_p,IERR_FATAL_p)
        END IF
      END IF

 */

class KineticsVocePL
{
public:
   static const int nH = 1 ;
   static const int nParams = 3 + 5 + nH ;
   static const int nVals = 1 ;
   
   // constructor
   __ecmech_hdev__
   KineticsVocePL(int nslip) : _nslip(nslip) {};

   __ecmech_host__
   inline void setParams( const real8* const params ) {

      int iParam = 0 ;

      //////////////////////////////
      // power-law stuff
      
      _mu     = params[iParam++] ;
      _xm     = params[iParam++] ;
      _gam_w  = params[iParam++] ;

      //    CALL fill_power_law(pl)
      // xmm  = xm - one ;
      _xnn       = one / _xm ;
      _xn        = _xnn - one ;
      // xMp1 = xnn + one
      //
      //    CALL set_t_min_max(pl)
      _t_min = pow(ecmech::gam_ratio_min, _xm) ;
      _t_max = pow(ecmech::gam_ratio_ovf, _xm) ;
       
      //////////////////////////////
      // Voce hardening stuff

      _h0     = params[iParam++] ;
      _tausi  = params[iParam++] ;
      _taus0  = params[iParam++] ;
      _xms    = params[iParam++] ;
      _gamss0 = params[iParam++] ;

      //////////////////////////////
      // nH

      _hdn_init  = params[iParam++] ;

      //////////////////////////////
      
      assert( iParam == nParams );
      
   };
   
   __ecmech_host__
   void getHistInfo(std::vector<std::string> & names,
                    std::vector<real8>       & init,
                    std::vector<bool>        & plot,
                    std::vector<bool>        & state) {
      names.push_back("h") ;
      init.push_back(_hdn_init) ;
      plot.push_back(true) ;
      state.push_back(true) ;
   }
   
private:

   const int _nslip ; // could template on this if there were call to do so

   // static const _nXnDim = nH*nH ; // do not bother

   //////////////////////////////
   // power-law stuff
   
   // parameters
   real8 _mu ; // may evetually set for current conditions
   real8 _xm ;
   real8 _gam_w ; // pl%adots, adots0

   // derived from parameters
   real8 _t_max, _t_min, _xn, _xnn ;

   //////////////////////////////
   // Voce hardening stuff

   real8 _h0, _tausi, _taus0, _xms, _gamss0 ;

   //////////////////////////////
   
   real8 _hdn_init ;
   
public:

__ecmech_hdev__
inline real8 getFixedRefRate(const real8* const // vals, not used
                             ) const
{
   return _gam_w ;
}

__ecmech_hdev__
inline
void
getVals( real8* const vals,
         real8 , // p, not currently used
         real8 , // tK, not currently used
         const real8* const h_state
         ) const
{
   vals[0] = h_state[0] ; // _gAll
   assert(vals[0] > zero);
}

__ecmech_hdev__
inline
void
evalGdots( real8* const gdot,
           real8* const dgdot_dtau,
           real8* const dgdot_dg,
           const real8* const tau,
           const real8* const vals
           ) const
{
   real8 gAll = vals[0] ; // gss%h(islip) // _gAll
   for ( int iSlip=0; iSlip<this->_nslip; ++iSlip ) {
      bool l_act ;
      this->evalGdot( gdot[iSlip], l_act, dgdot_dtau[iSlip], dgdot_dg[iSlip],
                      gAll,
                      tau[iSlip],
                      _mu // gss%ctrl%mu(islip)
                      ) ;
   }
}

/**
 * see kinetics_pl_d
 */
__ecmech_hdev__
inline
void
evalGdot(
                          real8 & gdot,
                          bool  & l_act,
                          real8 & dgdot_dtau,  // wrt resolved shear stress
                          real8 & dgdot_dg,    // wrt slip system strength
#if MORE_DERIVS
                          real8 & dgdot_dmu,   // wrt shear modulus, not through g
                          real8 & dgdot_dgamo, // wrt reference rate for thermal part
                          real8 & dgdot_dgamr, // wrt reference rate for drag limited part
                          real8 & dgdot_dtK,   // wrt temperature, with other arguments fixed
#endif
                          real8   gIn,
                          real8   tau,
                          real8   // mu, not currently used
#if MORE_DERIVS
                          ,
                          real8   tK
#endif
                         ) const
{
   //  zero things so that can more easily just return in inactive
   // // gdot_w = zero; gdot_r = zero; ! not used by l_linear or l_pl
   gdot = zero ;
   //
   dgdot_dtau  = zero ;
   dgdot_dg    = zero ;
#if MORE_DERIVS
   dgdot_dmu   = zero ;
   dgdot_dgamo = zero ;
   dgdot_dgamr = zero ;
   dgdot_dtK   = zero ;
#endif
   l_act = false ;

   real8 g_i = one / gIn ; // assume have checked gIn>0 elsewhere
   real8 t_frac = tau * g_i ; // has sign of tau
   real8 at = fabs(t_frac) ;

   if ( at > _t_min) {
      //
      l_act = true ;
      
      if ( at > _t_max) {
         // ierr = IERR_OVF_p 
         // set gdot big, evpp may need this for recovery
         gdot = ecmech::gam_ratio_ovffx * _gam_w ; 
         gdot = copysign(gdot, tau) ;
         // do not set any of deriviatives (they are, in truth, zero)
      }
      else {

         real8 abslog = log(at) ;
         real8 blog = _xn * abslog ;
         real8 temp = _gam_w * exp(blog) ; 

         gdot = temp * t_frac ;

         dgdot_dtau  =  temp * _xnn * g_i ; // note: always positive, = xnn * gdot/t
         dgdot_dg    = -dgdot_dtau * t_frac ; // = - gdot * xnn * g_i
#if MORE_DERIVS
         // dgdot_dmu   =  zero ; // already done
         dgdot_dgamo =  gdot / _gam_w ;
         // dgdot_dgamr = zero ; // already done
         // dgdot_dtK   = zero ; // already done
         
        // IF (pl%tdp%T_dep) THEN
        //   dgdot_dtK   =  gdot * pl%tdp%qoverr / (tK * tK)
        // END IF
        // IF (pl%tdp%T_dep_m) THEN !  .AND. pl%xm < XM_UB_p
        //   dxn_dtK = -(pl%xnn*pl%xnn) * pl%tdp%dxm_dtK_current
        //   dgdot_dtK = dgdot_dtK + &
        //        & gdot * abslog * dxn_dtK
        // END IF
#endif
      }

   }
   
} // evalGdot

__ecmech_hdev__
inline
void
updateH( real8* const hs_n,
         real8* const hs_o,
         real8 dt,
         const real8* const gdot ) ; 

__ecmech_hdev__
inline
void
getSdot( real8 &sdot,
         real8 &dsdot_ds,
         real8 h,
         real8 sv_sat,
         real8 shrate_eff) const {
   
   real8 temp2 = sv_sat - _tausi ;
   
   // IF (PRESENT(dfdtK)) THEN
   //   dfdtK(1) = zero
   // END IF

   sdot = 0.0 ;
   dsdot_ds = 0.0 ;
   //
   if ( temp2 <= zero ) {
      // most likely, shrate is small enough that hs_sat has become
      // small, but for small shrate have small rate of h change ;
      //
      // sdot and dsdot_ds already set to zero
   }
   else {
      real8 temp1 = _h0 * ((sv_sat - h) / temp2) ;
      sdot = temp1 * shrate_eff ;
      // real8 dfdshr = temp1 + _h0 * ( (h - _tausi) / (temp2*temp2)) * _xms * sv_sat ;
      dsdot_ds = - _h0 / temp2 * shrate_eff ;
   }
   
}

}; // class KineticsVocePL

class KineticsVocePL_HProblem
{
public:
   static const int nDimSys = KineticsVocePL::nH ;
   
   // constructor
   __ecmech_hdev__
   KineticsVocePL_HProblem(const KineticsVocePL* kinetics,
                           real8 h_o,
                           real8 sv_sat,
                           real8 dt,
                           real8 shrate_eff) :
      _kinetics(kinetics), _h_o(h_o), _sv_sat(sv_sat), _dt(dt), _shrate_eff(shrate_eff)
   {
      _x_scale   = _h_o ;
      _res_scale = one / _h_o ;
   }

   __ecmech_hdev__
   inline
   real8 getHn( const real8* const x ) const {
      return _h_o + x[0] * _x_scale ;
   }
   
   __ecmech_hdev__
   inline
   bool computeRJ( real8* const resid,
                   real8* const Jacobian,
                   const real8* const x ) {
      bool doComputeJ = (Jacobian != nullptr) ;
   
      real8 h_delta = x[0] * _x_scale ;
      real8 h = _h_o + h_delta ;

      real8 sdot, dsdot_ds ;
      _kinetics->getSdot(sdot, dsdot_ds, h, _sv_sat, _shrate_eff) ;
   
      resid[0] = (h_delta - sdot * _dt) * _res_scale ;

      if ( doComputeJ ) {
         Jacobian[0] = (one - dsdot_ds * _dt) * _res_scale * _x_scale ;
      }

      return true ;
      
   } // computeRJ

private :
   const KineticsVocePL* _kinetics ;
   const real8 _h_o, _sv_sat, _dt, _shrate_eff ;
   real8 _x_scale, _res_scale ; 
   
}; // class KineticsVocePL_HProblem

__ecmech_hdev__
inline
void
KineticsVocePL::updateH( real8* const hs_n,
                         real8* const hs_o,
                         real8 dt,
                         const real8* const gdot ) 
{

   // recompute effective shear rate here versus using a stored value
   real8 shrate_eff = vecsssumabs_n(gdot, _nslip) ; // could switch to template if template class on _nslip

   real8 sv_sat = _taus0 ;
   if ( shrate_eff > ecmech::idp_tiny_sqrt ) {
      sv_sat = _taus0 * pow((shrate_eff / _gamss0 ), _xms) ;
   }
   
   KineticsVocePL_HProblem prob(this, hs_o[0], sv_sat, dt, shrate_eff) ;
   snls::SNLSTrDlDenseG<KineticsVocePL_HProblem> solver(prob) ;

   snls::TrDeltaControl deltaControl ;
   deltaControl._deltaInit = 1e0 ;
   {
      int maxIter = 100 ;
      real8 tolerance = 1e-10 ;
      solver.setupSolver(maxIter, tolerance, &deltaControl) ;
   }
   
   real8* x = solver.getXPntr() ;
   for (int iX = 0; iX < prob.nDimSys; ++iX) {
      x[iX] = 0e0 ;
   }
   
   snls::SNLSStatus_t status = solver.solve( ) ;
   if ( status != snls::converged ){
      ECMECH_FAIL(__func__,"Solver failed to converge!");
   }

   hs_n[0] = prob.getHn(x) ;
         
}

} // namespace ecmech

#endif  // ECMECH_KINETICS_VOCEPL_H
