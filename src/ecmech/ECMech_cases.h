// -*-c++-*-

#include "SNLS_TrDLDenseG.h"

#include "ECMech_evptn.h"
#include "ECMech_evptnWrap.h"
#include "ECMech_kinetics.h"
#include "ECMech_slipgeom.h"

namespace ecmech {

   typedef KineticsKMBalD<true,false,false>  Kin_KMBalD_TFF ;
   typedef KineticsKMBalD<false,false,false> Kin_KMBalD_FFF ;

   typedef evptn::EvptnUpdstProblem< SlipGeomFCC, KineticsVocePL, evptn::ThermoElastNCubic > EvptnUpsdtProblem_FCC_A ;
   typedef snls::SNLSTrDlDenseG<EvptnUpsdtProblem_FCC_A> EvptnSolver_FCC_A ;
   typedef evptn::matModel< SlipGeomFCC, KineticsVocePL, evptn::ThermoElastNCubic, EosModelConst<false> > matModelEvptn_FCC_A ;

   typedef evptn::EvptnUpdstProblem< SlipGeomFCC, Kin_KMBalD_FFF, evptn::ThermoElastNCubic > EvptnUpsdtProblem_FCC_B ;
   typedef snls::SNLSTrDlDenseG<EvptnUpsdtProblem_FCC_B> EvptnSolver_FCC_B ;
   typedef evptn::matModel< SlipGeomFCC, Kin_KMBalD_FFF, evptn::ThermoElastNCubic, EosModelConst<false> > matModelEvptn_FCC_B ;
   
}
