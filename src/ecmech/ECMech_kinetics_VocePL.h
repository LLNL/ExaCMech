// -*-c++-*-

#ifndef ECMECH_KINETICS_VOCEPL_H
#define ECMECH_KINETICS_VOCEPL_H

#include <cassert>

#include "ECMech_core.h"

namespace ecmech {

/**
 * slip and hardening kinetics
 *
 * power-law slip kinetics with Voce hardening law -- meant to be about as simple as it gets
 */
class KineticsVocePL
{
public:
   static const int nH = 1 ;
   static const int nParams = 3 ;

   // constructor
   __ecmech_hdev__
   KineticsVocePL(int nslip) : _nslip(nslip) {};

   __ecmech_hdev__
   inline void setParams( const real8* const params ) {

      int iParam = 0 ;
      _mu     = params[iParam++] ;
      _xm     = params[iParam++] ;
      _gam_w  = params[iParam++] ;
      assert( iParam == nParams );

      //    CALL fill_power_law(pl)
      // xmm  = xm - one ;
      _xnn       = one / _xm ;
      _xn        = _xnn - one ;
      // xMp1 = xnn + one
      //
      //    CALL set_t_min_max(pl)
      _t_min = pow(ecmech::gam_ratio_min, _xm) ;
      _t_max = pow(ecmech::gam_ratio_ovf, _xm) ;
       
   };
   
private:

   const int _nslip ; // could template on this if there were call to do so

   // parameters
   real8 _mu ; // may evetually set for current conditions
   real8 _xm ;
   real8 _gam_w ; // pl%adots, adots0

   // derived from parameters
   real8 _t_max, _t_min, _xn, _xnn ;

   // current values
   real8 _g ;
   
public:

__ecmech_hdev__
inline
void
setVals( real8 p,
         real8 tK,
         const real8* const h_state
         )
{
   _g = h_state[0] ;
   assert(_g > zero);
}

__ecmech_hdev__
inline
void
evalGdots( real8* const gdot,
           real8* const dgdot_dtau,
           real8* const dgdot_dg,
           const real8* const tau
           ) const
{
   for ( int iSlip=0; iSlip<this->_nslip; ++iSlip ) {
      bool l_act ;
      this->evalGdot( gdot[iSlip], l_act, dgdot_dtau[iSlip], dgdot_dg[iSlip],
                      _g, // gss%h(islip)
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
                          real8   mu
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

   real8 g_i = one / _g ; // assume have checked _g>0 elsewhere
   real8 t_frac = tau * g_i ; // has sign of tau
   real8 at = fabs(t_frac) ;

   if ( at > _t_min) {
      //
      l_act = true ;
      
      if ( at > _t_max) {
         // ierr = IERR_OVF_p 
         // set gdot big, evpp may need this for recovery
         gdot = ecmech::gam_ratio_ovffx * _gam_w ; 
         gdot = std::copysign(gdot, tau) ;
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

}; // class KineticsVocePL

} // namespace ecmech

#endif  // ECMECH_KINETICS_VOCEPL_H
