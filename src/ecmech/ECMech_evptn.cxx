#include "SNLS_TrDLDenseG.h"

#include "ECMech_evptn.h"
#include "ECMech_kinetics.h"
#include "ECMech_slipgeom.h"

namespace ecmech {

   // some convenience stuff
   //
   typedef EvptnUpdstProblem< SlipGeomFCCA, KineticsVocePL, ThermoElastNCubic > EvptnUpsdtProblem_FCCAVocePl ; 
   typedef snls::SNLSTrDlDenseG<EvptnUpsdtProblem_FCCAVocePl> EvptnSolver_FCCAVocePl ;

} // namespace ecmech
