// -*-c++-*-

#ifndef ECMECH_KINETICS_KMBALD_H
#define ECMECH_KINETICS_KMBALD_H

#include "ECMech_const.h"

namespace ecmech {

/**
 * slip and hardening kinetics
 * based on a single Kocks-Mecking dislocation density
 * balanced thermally activated MTS-like slip kinetics with phonon drag effects
 *
 * if withGAthermal then 
 *	 see subroutine kinetics_mtspwr_d in mdef : (l_mts, l_mtsp, l_plwr)
 *	   ! like kinetics_mtswr_d, but with pl%tau_a (possible associated
 *	   ! with the Peierls barrier) being the thermally activated part and
 *	   ! g being athermal
 * 	   ! use balanced and pegged MTS model;
 * 	   ! add to it a low rate sensitivity power law model to take over for high stresses;
 * 	   ! and combine with drag limited kinetics
 * else then see subroutine kinetics_mtswr_d in mdef
 *
 *   ! note: gdot_w, gdot_r are always positive by definition
 *   !
 *   ! tK should only be used for derivative calculations
 *
 * templated on p and q being 1 or not;
 * might eventually template on number of slip systems, but to not do so just yet
 
 TODO : ...*** state update for hardness vector?
 
 */
template< bool withGAthermal, bool pOne, bool qOne > // l_p_1, l_q_1
class KineticsKMBalD
{
public:
   static const int nH = 1 ;
   static const int nParams = 11 ;

   // constructor
   __ecmech_hdev__
   KineticsKMBalD(int nslip) : _nslip(nslip) {};

   __ecmech_hdev__
   void setParams( const real8* const params ) {

      int iParam = 0 ;
      _mu     = params[iParam++] ;
      _tK_ref = params[iParam++] ;
      _c_1    = params[iParam++] ;
      _tau_a  = params[iParam++] ;
      _p      = params[iParam++] ;
      _q      = params[iParam++] ;
      _gam_wo = params[iParam++] ;
      _gam_ro = params[iParam++] ;
      _wrD    = params[iParam++] ;
      _go     = params[iParam++] ;
      _s      = params[iParam++] ;
      assert( iParam == nParams );

      if ( pOne ) {
         assert(_p == one);
      }
      if ( qOne ) {
         assert(_q == one);
      }
      
      // plaw_from_elawRef
      //
      //    pl%xm = getMtsxmEffective(pl, mu_ref, T_ref)
      real8 xm = one / (two * ((_c_1 / tK_ref) * mu_ref * _p * _q)) ;
      //
      //    CALL fill_power_law(pl)
      // xmm  = xm - one ;
      _xnn       = one / xm ;
      _xn        = _xnn - one ;
      // xMp1 = xnn + one
      //
      //    CALL set_t_min_max(pl)
      _t_min = pow(ecmech::gam_ratio_min, xm) ;
      _t_max = pow(ecmech::gam_ratio_ovf, xm) ;
       
   };
   
private:

   const int _nslip ; // could template on this if there were call to do so

   // parameters
   real8 _mu ; // may evetually set for current conditions
   real8 _tK_ref ;
   real8 _tau_a ; // if withGAthermal the is Peierls barrier
   real8 _p ; // only used if pOne is false
   real8 _q ; // only used if qOne is false
   real8 _gam_ro ;
   real8 _gam_wo ; // adots0
   real8 _c_1 ;
   real8 _wrD ;
   real8 _go, _s ;

   // derived from parameters
   real8 _t_max, _t_min, _xn, _xnn ;

   // current values
   real8 _g ;
   real8 _gam_r, gam_w ;
   real8 _c_t ;

public:

/**
 * @brief Akin to hs_to_gss, power_law_tdep_vals, and plaw_from_hs
 *
 * Could eventually bring in additional pressure and temperature dependence through the dependence of _mu on such ;
 * see use of mu_factors in Fortran code
 */
__ecmech_hdev__
inline
void
setVals( real8 p,
         real8 tK,
         const real8* const h_state
         )
{
   _c_t = _c_1 / tK ;
   real8 sqrtDDens = exp(onehalf * h_state[0]) ;
   _g = _go + _s * sqrtDDens ;
   _gam_w = _gam_wo / sqrtDDens ;
   _gam_r = _gam_ro * sqrtDDens * sqrtDDens ;
   if ( withGAthermal ) {
      assert(_tau_a > 0) ;
   } else {
      assert(_g > zero);
   }
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
 * like mts_dG, but with output args first
 */
__ecmech_hdev__
inline
void
get_mts_dG(real8 &exp_arg,
           real8 &mts_dfac,
           real8 c_e, real8 denom_i, real8 t_frac) {

   mts_dfac = c_e * denom_i ;

   real8 p_func ;
   if ( pOne ) {
      p_func = t_frac ;
   }
   else { 
      if ( fabs(t_frac) < idp_tiny_sqrt ) {
        // !! p_dfac is either zero or blows up
        // !IF (pl%p > one) THEN ! no longer allowed
        // !   mts_dfac = zero
        // !ELSE
        // ! blows up, but just set big
        mts_dfac = mts_dfac * 1e10 ;
        // !END IF
      }
      else {
         p_func = pow( fabs(t_frac), _p) ;
         p_func = std::copysign(p_func, t_frac) ;
         mts_dfac = mts_dfac * 
            _p * p_func / t_frac ; // always positive
      }
   }

   real8 q_arg = one - p_func ;
   real8 pq_fac ;
   if ( q_arg < idp_tiny_sqrt) {
      //  peg
      q_arg = zero ;
      mts_dfac = zero ;
      pq_fac = zero ;
   }
   else {
      if ( qOne ) {
         pq_fac = q_arg ;
      }
      else {
        real8 temp = pow( fabs(q_arg), _q ) ;
        mts_dfac = mts_dfac * 
           pl%q * temp / fabs(q_arg) ; // always positive
        pq_fac = std:copysign( temp, q_arg ) ;
      }
   }

    exp_arg = -c_e * pq_fac ;
    
}
   
/**
 * see subroutine kinetics_mtspwr_d in mdef 
 */
__ecmech_hdev__
inline
void
evalGdot(
                          real8 & gdot,
                          bool  & l_act
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
                          rael8   mu,
#if MORE_DERIVS
                          rael8   tK,
#endif
                         ) const
{
   static const real8 gdot_w_pl_scaling = 10.0 ;
   static const real8 one = 1.0, zero=0.0 ;

   // zero things so that can more easily just return if inactive
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

   real8 g_i ;
   real8 gAth ;
   if ( withGAthermal ) {
      gAth = gIn ;
      g_i = one / _tau_a ;
   }
   else {
      gAth = _tau_a ;
      if (tau == zero) {
         return ;
      }
      g_i = one / gIn ;
   }
   real8 at_0 = fmax(zero,fabs(tau) - gAth) * g_i ;
   
   // calculate drag limited kinetics
   //
   real8 gdot_r, dgdot_r ;
#if MORE_DERIVS
   real8 dgdotr_dtK ;
#endif
   {
      real8 exp_arg = (fabs(tau) - gAth)/_wrD ;
      real8 temp ;
      if ( exp_arg < gam_ratio_min) { // ! IF (gdot_r < gam_ratio_min) THEN
         //  note that this should catch tau <= g
         return ;
      }
      else if (exp_arg < idp_eps_sqrt) {
         // linear expansion is cheaper and more accurate
         gdot_r = _gam_r * exp_arg ;
         temp = one - exp_arg ; // still use temp below as approximation to exp(-fabs(tau)/_wrD)
      }
      else { 
         temp = exp(-exp_arg) ;
         gdot_r = _gam_r * (one - temp) ;
      }
      real8 dgdot_r = _gam_r * temp / _wrD ;
#if MORE_DERIVS
      real8 dgdotr_dtK ;
      if ( withGAthermal ) { 
         dgdotr_dtK = -_gam_r * temp * exp_arg * _wrDT / _wrD ;
      } else {
         dgdotr_dtK = -_gam_r * temp * fabs(tau) * _wrDT / (_wrD * _wrD) ;
      }
#endif
   }
   //
   if (at_0 > _t_max) {
      // have overflow of thermally activated kinetics, purely drag limited

      gdot        =  gdot_r ;

      dgdot_dtau  =  dgdot_r ;
      if ( withGAthermal ) {
         dgdot_dg =  zero ;
      } else {
         dgdot_dg = -std::copysign(dgdot_r, tau) ;
      }
#if MORE_DERIVS
      dgdot_dmu   =  zero ;
      dgdot_dgamo =  zero ;
      dgdot_dtK   =  std::copysign(dgdotr_dtK, tau) ;
      dgdot_dgamr =  std::copysign(gdot, tau) / _gam_r ;
#endif
      gdot   = std::copysign(gdot,tau) ;

      l_act = true ;
      return ;

   }

   real8 gdot_w, dgdot_w ;
   real8 dgdot_wg ; // only used if !withGAthermal
   //
   // calculate thermally activated kinetics
   {
      real8 c_e = _c_t * mu ;
      //
      real8 t_frac = (fabs(tau) - gAth) * g_i ;
      real8 exp_arg, mts_dfac ;
      get_mts_dG(exp_arg, mts_dfac, c_e, g_i, t_frac) ;
      // 
      if ( exp_arg < ln_gam_ratio_min ) {
         // effectively zero due to thermally activated kinetics
         l_act = false ;
         return ;
      }
      //
#if MORE_DERIVS
      dgdotw_dmu = zero ;
      dgdotw_dtK = zero ;
#endif         
      //
      // !IF (exp_arg > ln_gam_ratio_ovf) THEN
      // !END IF
      // ! do not need to check the above condition because have pegged the MTS part of the kinetics
      // 
      gdot_w = _gam_w * exp(exp_arg) ;
      dgdot_w = mts_dfac * gdot_w ;
      if ( !withGAthermal ) {
         dgdot_wg = dgdot_w * t_frac ;
      }
#if MORE_DERIVS
      dgdotw_dmu = gdot_w * (-exp_arg/mu) ;
      dgdotw_dtK = gdot_w * ( exp_arg/tK) ; // negatives cancel
#endif
      //
      real8 t_frac_m = (-fabs(tau) - gAth) * g_i ;
      real8 exp_arg_m, mts_dfac_m ;
      get_mts_dG(exp_arg_m, mts_dfac_m, c_e, g_i, t_frac_m) ;
      // 
      if ( exp_arg_m > ln_gam_ratio_min ) {
         // non-vanishing contribution from balancing MTS-like kinetics
         real8 gdot_w_m = _gam_w * exp(exp_arg_m) ;
         gdot_w = gdot_w - gdot_w_m ;
         real8 contrib = mts_dfac_m * gdot_w_m ;
         dgdot_w = dgdot_w - contrib ; // sign used to be the other way, but suspect that was a bug
         if ( !withGAthermal ) {
            dgdot_wg = dgdot_wg - contrib * t_frac_m ;
         }
         if ( fabs(gdot_w/_gam_w) < gam_ratio_min ) {
            // effectively zero from roundoff
            l_act = false ;
            return ;
         }
      }
   }

   if ( at_0 > _t_min) {
      // need power-law part

      real8 abslog = log(at_0) ;
      real8 blog = _xn * abslog ;
      real8 temp = (_gam_w * gdot_w_pl_scaling) * exp(blog) ;

      real8 gdot_w_pl = temp * at_0 ; // not signed ! std::copysign(at_0,tau) 
      gdot_w = gdot_w + gdot_w_pl ;

      real8 contrib = temp * _xnn * g_i ;
      dgdot_w = dgdot_w + contrib ;
      if ( !withGAthermal ) {
         dgdot_wg = dgdot_wg + contrib * at_0 ;
      }

   }

   l_act = true ;
   //
   {
      gdot   = one/(one/gdot_w + one/gdot_r) ;
      real8 gdrdiv2 = one/(gdot_r*gdot_r) ;
      real8 gdwdiv2 = one/(gdot_w*gdot_w) ;
      dgdot_dtau  = (gdot*gdot)*( dgdot_w*gdwdiv2 + dgdot_r*gdrdiv2 ) ;
      //
      real8 temp  =   gdot * std::copysign(gdot,tau) * gdwdiv2 ;
      // neglect difference in at_0 versus t_frac for dgdot_dg evaluation
      if ( withGAthermal ) {
         dgdot_dg = - temp * dgdot_w ;  // opposite sign as signed gdot
      } else {
         dgdot_dg = - temp * dgdot_wg ; // opposite sign as signed gdot
      }
#if MORE_DERIVS
      dgdot_dgamo =   temp * (gdot_w / _gam_w) ;
      dgdot_dmu   =   temp * dgdotw_dmu ;
      dgdot_dtK   =   temp * dgdotw_dtK ;
#endif
      //
      temp        =   gdot * std::copysign(gdot,tau) * gdrdiv2 ;
      if ( withGAthermal ) {
         dgdot_dg = - temp * dgdot_r + dgdot_dg ; // opposite sign as signed gdot
      }
#if MORE_DERIVS
      dgdot_dgamr =   temp * (gdot_r / pl%gam_r) ;
      dgdot_dtK   = dgdot_dtK + temp * dgdotr_dtK ;
#endif
   }

   gdot   = std::copysign(gdot,tau) ;
   
} // evalGdot

} // class KineticsKMBalD

} // namespace ecmech

#endif  // ECMECH_KINETICS_KMBALD_H
