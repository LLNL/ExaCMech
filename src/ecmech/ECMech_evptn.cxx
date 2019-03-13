#include "SNLS_TrDLDenseG.h"
#include "ECMech_evptn.h"

namespace ecmech {

// some convenience stuff
//
// case for some specific kinetics, symmetry, whatever ...
typedef EvptnUpdstProblem<...***> EvptnAProblem ; 
typedef snls::SNLSTrDlDenseG<EvptnAProblem> EvptnASolver ;

} // namespace ecmech
