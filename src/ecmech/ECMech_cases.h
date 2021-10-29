// -*-c++-*-

#include "SNLS_TrDLDenseG.h"

#include "ECMech_evptn.h"
#include "ECMech_evptnWrap.h"
#include "ECMech_kinetics.h"
#include "ECMech_slipgeom.h"

namespace ecmech {
   typedef KineticsKMBalD<true, false, false, false, 1> Kin_KMBalD_TFF;
   typedef KineticsKMBalD<false, false, false, false, 1> Kin_KMBalD_FFF;

   typedef KineticsOrowanD<false, false, false, true, false, 1, SlipGeomFCC> Kin_OroD_Iso_FCC;
   typedef KineticsOrowanD<true, false, false, true, false, 1, SlipGeomBCC<12>> Kin_OroD_Iso_BCC;
   typedef KineticsOrowanD<true, false, false, false, false, 1, SlipGeomBCC<12>> Kin_OroD_Aniso_BCC;
   typedef KineticsOrowanD<true, false, false, true, false, 1, SlipGeomBCC<24>> Kin_OroD_Iso_BCC_24;
   typedef KineticsOrowanD<true, false, false, false, false, 1, SlipGeomBCC<24>> Kin_OroD_Aniso_BCC_24;

   typedef KineticsVocePL<false> Kin_FCC_A;
   typedef evptn::EvptnUpdstProblem<SlipGeomFCC, Kin_FCC_A, evptn::ThermoElastNCubic> EvptnUpsdtProblem_FCC_A;
   typedef snls::SNLSTrDlDenseG<EvptnUpsdtProblem_FCC_A> EvptnSolver_FCC_A;
   typedef evptn::matModel<SlipGeomFCC, Kin_FCC_A, evptn::ThermoElastNCubic, EosModelConst<false> > matModelEvptn_FCC_A;

   typedef KineticsVocePL<true> Kin_FCC_AH;
   typedef evptn::EvptnUpdstProblem<SlipGeomFCC, Kin_FCC_AH, evptn::ThermoElastNCubic> EvptnUpsdtProblem_FCC_AH;
   typedef snls::SNLSTrDlDenseG<EvptnUpsdtProblem_FCC_AH> EvptnSolver_FCC_AH;
   typedef evptn::matModel<SlipGeomFCC, Kin_FCC_AH, evptn::ThermoElastNCubic, EosModelConst<false> > matModelEvptn_FCC_AH;

   typedef Kin_KMBalD_FFF Kin_FCC_B;
   typedef evptn::EvptnUpdstProblem<SlipGeomFCC, Kin_FCC_B, evptn::ThermoElastNCubic> EvptnUpsdtProblem_FCC_B;
   typedef snls::SNLSTrDlDenseG<EvptnUpsdtProblem_FCC_B> EvptnSolver_FCC_B;
   typedef evptn::matModel<SlipGeomFCC, Kin_FCC_B, evptn::ThermoElastNCubic, EosModelConst<false> > matModelEvptn_FCC_B;

   typedef Kin_OroD_Iso_FCC Kin_FCC_C;
   typedef evptn::EvptnUpdstProblem<SlipGeomFCC, Kin_FCC_C, evptn::ThermoElastNCubic> EvptnUpsdtProblem_FCC_C;
   typedef snls::SNLSTrDlDenseG<EvptnUpsdtProblem_FCC_B> EvptnSolver_FCC_C;
   typedef evptn::matModel<SlipGeomFCC, Kin_FCC_C, evptn::ThermoElastNCubic, EosModelConst<false> > matModelEvptn_FCC_C;

   typedef SlipGeomBCC<12> SlipGeom_BCC_A;
   typedef SlipGeomBCC<24> SlipGeom_BCC_B;
   typedef Kin_KMBalD_TFF Kin_BCC_A;
   typedef evptn::EvptnUpdstProblem<SlipGeom_BCC_A, Kin_BCC_A, evptn::ThermoElastNCubic> EvptnUpsdtProblem_BCC_A;
   typedef snls::SNLSTrDlDenseG<EvptnUpsdtProblem_BCC_A> EvptnSolver_BCC_A;
   typedef evptn::matModel<SlipGeom_BCC_A, Kin_BCC_A, evptn::ThermoElastNCubic, EosModelConst<false> > matModelEvptn_BCC_A;

   typedef Kin_OroD_Iso_BCC Kin_BCC_B;
   typedef evptn::EvptnUpdstProblem<SlipGeom_BCC_A, Kin_BCC_B, evptn::ThermoElastNCubic> EvptnUpsdtProblem_BCC_B;
   typedef snls::SNLSTrDlDenseG<EvptnUpsdtProblem_BCC_B> EvptnSolver_BCC_B;
   typedef evptn::matModel<SlipGeom_BCC_A, Kin_BCC_B, evptn::ThermoElastNCubic, EosModelConst<false> > matModelEvptn_BCC_B;
   typedef evptn::matModel<SlipGeom_BCC_A, Kin_OroD_Aniso_BCC, evptn::ThermoElastNCubic, EosModelConst<false> > matModelEvptn_BCC_C;
   typedef evptn::matModel<SlipGeom_BCC_B, Kin_OroD_Iso_BCC_24, evptn::ThermoElastNCubic, EosModelConst<false> > matModelEvptn_BCC_B_24;
   typedef evptn::matModel<SlipGeom_BCC_B, Kin_OroD_Aniso_BCC_24, evptn::ThermoElastNCubic, EosModelConst<false> > matModelEvptn_BCC_C_24;

   typedef SlipGeomHCPaBRYcaY1 SlipGeom_HCP_A;
   typedef KineticsKMBalD<true, true, true, true, SlipGeom_HCP_A::nslip> Kin_HCP_A;
   typedef evptn::EvptnUpdstProblem<SlipGeom_HCP_A,
                                    Kin_HCP_A,
                                    evptn::ThermoElastNHexag> EvptnUpsdtProblem_HCP_A;
   typedef snls::SNLSTrDlDenseG<EvptnUpsdtProblem_HCP_A> EvptnSolver_HCP_A;
   typedef evptn::matModel<SlipGeom_HCP_A,
                           Kin_HCP_A,
                           evptn::ThermoElastNHexag,
                           EosModelConst<false> > matModelEvptn_HCP_A;

   __ecmech_host__
   matModelBase* makeMatModel(const std::string &modelName);
}
