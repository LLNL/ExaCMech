// -*-c++-*-

#ifndef ECMECH_KINETICS_H
#define ECMECH_KINETICS_H

#include "SNLS_TrDLDenseG.h"

#include "ECMech_core.h"
#include "ECMech_util.h"

/**

// classes need to have member functions:

// constructor
__ecmech_hdev__
Kinetics(int nslip) : _nslip(nslip) {};

__ecmech_hdev__
inline real8 getFixedRefRate(const real8* const vals) const

__ecmech_hdev__
inline void setParams( const std::vector<real8> & params) // const real8* const params

// Akin to hs_to_gss, power_law_tdep_vals, and plaw_from_hs
//
__ecmech_hdev__
inline
void
getVals( real8* const vals,
         real8 p,
         real8 tK,
         const real8* const h_state
         ) const

__ecmech_hdev__
inline
void
evalGdots( real8* const gdot,
           real8* const dgdot_dtau,
           real8* const dgdot_dg,
           const real8* const tau,
           const real8* const vals
           ) const

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
                          real8   mu,
#if MORE_DERIVS
                          real8   tK,
#endif
                         ) const
           
 */

namespace ecmech {

template< class Kinetics >
class Kinetics_H1Problem
{
public:
   static const int nDimSys = Kinetics::nH ;
   
   // constructor
   __ecmech_hdev__
   Kinetics_H1Problem(const Kinetics* const kinetics,
                      real8 h_o,
                      real8 dt,
                      const real8* const evolVals) :
      _kinetics(kinetics), _h_o(h_o), _dt(dt), _evolVals(evolVals)
   {
      _x_scale   = fmax(_h_o, 1.0) ;
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
      _kinetics->getSdot1(sdot, dsdot_ds, h, _evolVals) ;
   
      resid[0] = (h_delta - sdot * _dt) * _res_scale ;

      if ( doComputeJ ) {
         Jacobian[0] = (one - dsdot_ds * _dt) * _res_scale * _x_scale ;
      }

      return true ;
      
   } // computeRJ

private :
   const Kinetics* _kinetics ;
   const real8 _h_o, _dt ;
   const real8* const _evolVals ;
   real8 _x_scale, _res_scale ; 
   
}; // class Kinetics_H1Problem

template< class Kinetics >
__ecmech_hdev__
inline
void
updateH1( const Kinetics* const kinetics,
          real8 &hs_n,
          real8 hs_o,
          real8 dt,
          const real8* const gdot ) 
{

   real8 evolVals[Kinetics::nEvolVals] ;
   kinetics->getEvolVals(evolVals, gdot) ;
   
   Kinetics_H1Problem<Kinetics> prob(kinetics, hs_o, dt, evolVals) ;
   snls::SNLSTrDlDenseG<Kinetics_H1Problem<Kinetics>> solver(prob) ;

   snls::TrDeltaControl deltaControl ;
   deltaControl._deltaInit = 1e0 ;
   {
      int maxIter = 100 ;
      real8 tolerance = 1e-10 ;
      solver.setupSolver(maxIter, tolerance, &deltaControl) ;
   }
   
   real8* x = solver.getXPntr() ;
   // for (int iX = 0; iX < prob.nDimSys; ++iX) {
   //    x[iX] = 0e0 ;
   // }
   x[0] = 0.0 ;
   
   snls::SNLSStatus_t status = solver.solve( ) ;
   if ( status != snls::converged ){
      ECMECH_FAIL(__func__,"Solver failed to converge!");
   }

   hs_n = prob.getHn(x) ;
         
}

} // namespace ecmech

#include "ECMech_kinetics_KMBalD.h"
#include "ECMech_kinetics_VocePL.h"

#endif  // ECMECH_KINETICS_H
