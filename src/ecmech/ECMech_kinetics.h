// -*-c++-*-

#ifndef ECMECH_KINETICS_H
#define ECMECH_KINETICS_H

#include "SNLS_TrDLDenseG.h"

#include "ECMech_core.h"
#include "ECMech_util.h"

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
          const real8* const gdot,
          int outputLevel = 0) 
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
      solver.setupSolver(maxIter, tolerance, &deltaControl, outputLevel) ;
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
