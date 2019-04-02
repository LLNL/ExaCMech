// -*-c++-*-

#ifndef ECMECH_KINETICS_VOCEPL_H
#define ECMECH_KINETICS_VOCEPL_H

#include <cassert>
#include <cmath>

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
   static const int nParams = 3 + 5 + nH ;
   static const int nVals = 1 ;
   static const int nEvolVals = 2 ;
   
   // constructor
   __ecmech_hdev__
   KineticsVocePL(int nslip) : _nslip(nslip) {};

   __ecmech_host__
   inline void setParams( const std::vector<real8> & params // const real8* const params
                          ) {

      std::vector<real8>::const_iterator parsIt = params.begin();

      //////////////////////////////
      // power-law stuff
      
      _mu     = *parsIt; ++parsIt;
      _xm     = *parsIt; ++parsIt;
      _gam_w  = *parsIt; ++parsIt;

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

      _h0     = *parsIt; ++parsIt;
      _tausi  = *parsIt; ++parsIt;
      _taus0  = *parsIt; ++parsIt;
      _xms    = *parsIt; ++parsIt;
      _gamss0 = *parsIt; ++parsIt;

      //////////////////////////////
      // nH

      _hdn_init  = *parsIt; ++parsIt;

      //////////////////////////////

      int iParam = parsIt - params.begin();
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
updateH( real8* const hs_u,
         const real8* const hs_o,
         real8 dt,
         const real8* const gdot,
         int outputLevel = 0  )
{
   real8 hs_u_1 ;
   updateH1<KineticsVocePL>(this,
                            hs_u_1, hs_o[0], dt, gdot,
                            outputLevel) ;
   hs_u[0] = hs_u_1 ;
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

   real8 sv_sat = _taus0 ;
   if ( shrate_eff > ecmech::idp_tiny_sqrt ) {
      sv_sat = _taus0 * pow((shrate_eff / _gamss0 ), _xms) ;
   }

   evolVals[0] = shrate_eff ;
   evolVals[1] = sv_sat ;
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
   real8 sv_sat     = evolVals[1] ;
   
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

} // namespace ecmech

#endif  // ECMECH_KINETICS_VOCEPL_H
