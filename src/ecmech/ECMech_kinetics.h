// -*-c++-*-

#ifndef ECMECH_KINETICS_H
#define ECMECH_KINETICS_H

/**

// classes need to have member functions:

// constructor
__ecmech_hdev__
Kinetics(int nslip) : _nslip(nslip) {};

__ecmech_hdev__
void setParams( const real8* const params ) ;

// Akin to hs_to_gss, power_law_tdep_vals, and plaw_from_hs
//
__ecmech_hdev__
inline
void
setVals( real8 p,
         real8 tK,
         const real8* const h_state
         )

__ecmech_hdev__
inline
void
evalGdots( real8* const gdot,
           real8* const dgdot_dtau,
           real8* const dgdot_dg,
           const real8* const tau
           ) const

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
           
 */
#include "ECMech_kinetics_KMBalD.h"
#include "ECMech_kinetics_VocePL.h"

#endif  // ECMECH_KINETICS_H
