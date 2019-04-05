// -*-c++-*-

#ifndef ECMECH_KINETICS_KMBALD_H
#define ECMECH_KINETICS_KMBALD_H

#include <cassert>
#include <cmath>

#include <string>
#include <vector>

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
   static const int nParams = 11 + 4 + nH ;
   static const int nVals = 4 ;
   static const int nEvolVals = 2 ;
   
   // constructor
   __ecmech_hdev__
   KineticsKMBalD(int nslip) : _nslip(nslip) {};

   __ecmech_hdev__
   void setParams( const std::vector<real8> & params // const real8* const params
                   ) {

      std::vector<real8>::const_iterator parsIt = params.begin();

      //////////////////////////////
      // power-law stuff
      
      _mu_ref = *parsIt; ++parsIt;
      _tK_ref = *parsIt; ++parsIt;
      _c_1    = *parsIt; ++parsIt;
      _tau_a  = *parsIt; ++parsIt;
      _p      = *parsIt; ++parsIt;
      _q      = *parsIt; ++parsIt;
      _gam_wo = *parsIt; ++parsIt;
      _gam_ro = *parsIt; ++parsIt;
      _wrD    = *parsIt; ++parsIt;
      _go     = *parsIt; ++parsIt;
      _s      = *parsIt; ++parsIt;

      if ( pOne ) {
         assert(_p == one);
      }
      if ( qOne ) {
         assert(_q == one);
      }
      
      // plaw_from_elawRef
      //
      //    pl%xm = getMtsxmEffective(pl, mu_ref, T_ref)
      real8 xm = one / (two * ((_c_1 / _tK_ref) * _mu_ref * _p * _q)) ;
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

      //////////////////////////////
      // Kocks-Mecking stuff

      _k1      = *parsIt; ++parsIt;
      _k2o     = *parsIt; ++parsIt;
      _ninv    = *parsIt; ++parsIt;
      _gamma_o = *parsIt; ++parsIt;

      //////////////////////////////
      // nH

      _hdn_init  = *parsIt; ++parsIt;

      _hdn_min   = 1e-4*_hdn_init ;

      //////////////////////////////
      
      int iParam = parsIt - params.begin();
      assert( iParam == nParams );
      
   };

   __ecmech_hdev__
   void getParams( std::vector<real8> & params
                   ) const {

      // do not clear params in case adding to an existing set
      int paramsStart = params.size() ;

      //////////////////////////////
      // power-law stuff
      
      params.push_back(_mu_ref ) ;
      params.push_back(_tK_ref ) ;
      params.push_back(_c_1    ) ;
      params.push_back(_tau_a  ) ;
      params.push_back(_p      ) ;
      params.push_back(_q      ) ;
      params.push_back(_gam_wo ) ;
      params.push_back(_gam_ro ) ;
      params.push_back(_wrD    ) ;
      params.push_back(_go     ) ;
      params.push_back(_s      ) ;

      //////////////////////////////
      // Kocks-Mecking stuff

      params.push_back(_k1      ) ;
      params.push_back(_k2o     ) ;
      params.push_back(_ninv    ) ;
      params.push_back(_gamma_o ) ;

      //////////////////////////////
      // nH

      params.push_back(_hdn_init ) ;

      //////////////////////////////
      
      int iParam = params.size() - paramsStart;
      assert( iParam == nParams );
      
   };
   
#ifdef __cuda_host_only__
   __ecmech_host__
   void getHistInfo(std::vector<std::string> & names,
                    std::vector<real8>       & init,
                    std::vector<bool>        & plot,
                    std::vector<bool>        & state) const {
      names.push_back("rho_dd") ;
      init.push_back(_hdn_init) ;
      plot.push_back(true) ;
      state.push_back(true) ;
   }
#endif
   
private:

   const int _nslip ; // could template on this if there were call to do so

   //////////////////////////////
   // MTS-like stuff
   
   // parameters
   real8 _mu_ref ; // may evetually set for current conditions
   real8 _tK_ref ;
   real8 _tau_a ; // if withGAthermal then is Peierls barrier
   real8 _p ; // only used if pOne is false
   real8 _q ; // only used if qOne is false
   real8 _gam_ro ;
   real8 _gam_wo ; // adots0
   real8 _c_1 ;
   real8 _wrD ;
   real8 _go, _s ;

   // derived from parameters
   real8 _t_max, _t_min, _xn, _xnn ;

   //////////////////////////////
   // Kocks-Mecking stuff

   real8 _k1, _k2o, _ninv, _gamma_o ;

   //////////////////////////////
   
   real8 _hdn_init, _hdn_min ;
   
public:

__ecmech_hdev__
inline
real8
getFixedRefRate(const real8* const vals) const
{
   
   return vals[1]+vals[2] ; // _gam_w + _gam_r ;
}

/**
 * @brief Akin to hs_to_gss, power_law_tdep_vals, and plaw_from_hs
 *
 * Could eventually bring in additional pressure and temperature dependence through the dependence of _mu on such ;
 * see use of mu_factors in Fortran code
 */
__ecmech_hdev__
inline
void
getVals( real8* const vals,
         real8 , // p, not used
         real8 tK,
         const real8* const h_state
         ) const
{
   vals[3] = _c_1 / tK ; // _c_t
   real8 sqrtDDens = sqrt(h_state[0]) ;
   // real8 sqrtDDens = exp(onehalf * h_state[0]) ; // this is for h_state[0] storing the log of the dislocation density
   vals[0] = _go + _s * sqrtDDens ; // _gAll
   vals[1] = _gam_wo / sqrtDDens ; // _gam_w
   vals[2] = _gam_ro * sqrtDDens * sqrtDDens ; // _gam_r
   if ( withGAthermal ) {
      assert(_tau_a > 0) ;
   } else {
      assert(vals[0] > zero);
   }
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
   for ( int iSlip=0; iSlip<this->_nslip; ++iSlip ) {
      bool l_act ;
      this->evalGdot( gdot[iSlip], l_act, dgdot_dtau[iSlip], dgdot_dg[iSlip],
                      vals,
                      tau[iSlip],
                      _mu_ref // gss%ctrl%mu(islip)
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
           real8 c_e, real8 denom_i, real8 t_frac) const {

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
         p_func = copysign(p_func, t_frac) ;
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
           _q * temp / fabs(q_arg) ; // always positive
        pq_fac = copysign( temp, q_arg ) ;
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
                          bool  & l_act,
                          real8 & dgdot_dtau,  // wrt resolved shear stress
                          real8 & dgdot_dg,    // wrt slip system strength
#if MORE_DERIVS
                          real8 & dgdot_dmu,   // wrt shear modulus, not through g
                          real8 & dgdot_dgamo, // wrt reference rate for thermal part
                          real8 & dgdot_dgamr, // wrt reference rate for drag limited part
                          real8 & dgdot_dtK,   // wrt temperature, with other arguments fixed
#endif
                          const real8* const vals,
                          real8   tau,
                          real8   mu
#if MORE_DERIVS
                          ,
                          real8   tK
#endif
                         ) const
{
   static const real8 gdot_w_pl_scaling = 10.0 ;
   static const real8 one = 1.0, zero=0.0 ;

   real8 gIn   = vals[0] ;
   real8 gam_w = vals[1] ;
   real8 gam_r = vals[2] ;
   real8 c_t   = vals[3] ;

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
         gdot_r = gam_r * exp_arg ;
         temp = one - exp_arg ; // still use temp below as approximation to exp(-fabs(tau)/_wrD)
      }
      else { 
         temp = exp(-exp_arg) ;
         gdot_r = gam_r * (one - temp) ;
      }
      dgdot_r = gam_r * temp / _wrD ;
#if MORE_DERIVS
      real8 dgdotr_dtK ;
      if ( withGAthermal ) { 
         dgdotr_dtK = -gam_r * temp * exp_arg * _wrDT / _wrD ;
      } else {
         dgdotr_dtK = -gam_r * temp * fabs(tau) * _wrDT / (_wrD * _wrD) ;
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
         dgdot_dg = -copysign(dgdot_r, tau) ;
      }
#if MORE_DERIVS
      dgdot_dmu   =  zero ;
      dgdot_dgamo =  zero ;
      dgdot_dtK   =  copysign(dgdotr_dtK, tau) ;
      dgdot_dgamr =  copysign(gdot, tau) / gam_r ;
#endif
      gdot   = copysign(gdot,tau) ;

      l_act = true ;
      return ;

   }

   real8 gdot_w, dgdot_w ;
   real8 dgdot_wg ; // only used if !withGAthermal
   //
   // calculate thermally activated kinetics
   {
      real8 c_e = c_t * mu ;
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
      gdot_w = gam_w * exp(exp_arg) ;
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
         real8 gdot_w_m = gam_w * exp(exp_arg_m) ;
         gdot_w = gdot_w - gdot_w_m ;
         real8 contrib = mts_dfac_m * gdot_w_m ;
         dgdot_w = dgdot_w - contrib ; // sign used to be the other way, but suspect that was a bug
         if ( !withGAthermal ) {
            dgdot_wg = dgdot_wg - contrib * t_frac_m ;
         }
         if ( fabs(gdot_w/gam_w) < gam_ratio_min ) {
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
      real8 temp = (gam_w * gdot_w_pl_scaling) * exp(blog) ;

      real8 gdot_w_pl = temp * at_0 ; // not signed ! copysign(at_0,tau) 
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
      real8 temp  =   gdot * copysign(gdot,tau) * gdwdiv2 ;
      // neglect difference in at_0 versus t_frac for dgdot_dg evaluation
      if ( withGAthermal ) {
         dgdot_dg = - temp * dgdot_w ;  // opposite sign as signed gdot
      } else {
         dgdot_dg = - temp * dgdot_wg ; // opposite sign as signed gdot
      }
#if MORE_DERIVS
      dgdot_dgamo =   temp * (gdot_w / gam_w) ;
      dgdot_dmu   =   temp * dgdotw_dmu ;
      dgdot_dtK   =   temp * dgdotw_dtK ;
#endif
      //
      temp        =   gdot * copysign(gdot,tau) * gdrdiv2 ;
      if ( withGAthermal ) {
         dgdot_dg = - temp * dgdot_r + dgdot_dg ; // opposite sign as signed gdot
      }
#if MORE_DERIVS
      dgdot_dgamr =   temp * (gdot_r / gam_r) ;
      dgdot_dtK   = dgdot_dtK + temp * dgdotr_dtK ;
#endif
   }

   gdot   = copysign(gdot,tau) ;
   
} // evalGdot

__ecmech_hdev__
inline
void
updateH( real8* const hs_u,
         const real8* const hs_o,
         real8 dt,
         const real8* const gdot,
         int outputLevel = 0 ) const
{

   // do not yet both with l_overdriven and setting-to-saturation machinery as in Fortran coding
   
   // update is done on log(h) -- h treated as a nomralized (unitless) dislocation density
   real8 log_hs_u ;
   real8 log_hs_o = log(fmax(hs_o[0],_hdn_min)) ;
   updateH1<KineticsKMBalD>(this,
                            log_hs_u, log_hs_o, dt, gdot,
                            outputLevel) ;
   hs_u[0] = exp(log_hs_u) ;
   
}

__ecmech_hdev__
inline
void
getEvolVals( real8* const evolVals,
             const real8* const gdot
             ) const
{
   // recompute effective shear rate here versus using a stored value
   real8 shrate_eff = vecsssumabs_n(gdot, _nslip) ; // could switch to template if template class on _nslip

   real8 k2 = _k2o ;
   if ( shrate_eff > ecmech::idp_tiny_sqrt ) {
      k2 = _k2o * pow((_gamma_o/shrate_eff), _ninv) ;
   }

   evolVals[0] = shrate_eff ;
   evolVals[1] = k2 ;
}

__ecmech_hdev__
inline
void
getSdot1( real8 &sdot,
          real8 &dsdot_ds,
          real8 h,
          const real8* const evolVals) const
{

   real8 shrate_eff = evolVals[0] ;
   real8 k2         = evolVals[1] ;

   // IF (PRESENT(dfdtK)) THEN
   //   dfdtK(1) = zero
   // END IF

   // sdot = 0.0 ;
   // dsdot_ds = 0.0 ;
   //
   // if ( shrate_eff <= zero ) {
   //    // do not get any evolution, and will get errors if proceed with calculations below
   // }
   // else {
      real8 temp_hs_a = exp(-onehalf * h) ;
      real8 temp1 = _k1 * temp_hs_a - k2 ;
      sdot = temp1 * shrate_eff ;
      // dfdshr = temp1 + _ninv * k2 ;
      dsdot_ds = (-_k1 * onehalf * temp_hs_a) * shrate_eff ;
   // }
   
}

}; // class KineticsKMBalD

} // namespace ecmech

#endif  // ECMECH_KINETICS_KMBALD_H
