#include "SNLS_TrDLDenseG.h"
#include "ECMech_evptn.h"

namespace ecmech {

__ecmech_hdev__
bool
ECMechEvptnUpdstProblem::computeRJ( real8* const r,
                               real8* const J,
                               const real8* const x )
{

   bool doComputeJ = (J != NULL) ;
      
   if ( doComputeJ ) {

      // zero the Jacobian so that do not need to worry about zero
      // entries in the midst of other things later
      //
      for ( int ijJ=0; ijJ<_nXnDim; ++ijJ ) {
         J[ijJ] = 0.0;
      }
         
   }

   // ...*** ; do stuff

   // ~/mdef/build/matlib/matEvpc/evptn.F90 : eval_rJ
   //	... state_evptn_er_mod versus state_evptn_e_mod !!!
   // 	... different class instance for doing forward-Euler rotation? (or for keeping rotation fixed?)
   //		... using the same struct for the data?
   
}

// template snls::SNLSTrDlDenseG<ECMechEvptnUpdstProblem> {}; // not needed?

typedef snls::SNLSTrDlDenseG<ECMechEvptnUpdstProble> ECMechUpdstSolver ;

} // namespace ecmech
