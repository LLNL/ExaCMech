#include "SNLS_TrDLDenseG.h"

#include "ECMech_evptn.h"
#include "ECMech_cases.h"
#include "ECMech_kinetics.h"
#include "ECMech_slipgeom.h"
#include "ECMech_eosSimple.h"
#include "ECMech_util.h"

#ifndef KIN_TYPE
#define KIN_TYPE 1
#endif

int main(int argc, char *argv[])
{
   int outputLevel = 1 ;
   if ( argc > 1 ) {
      outputLevel = atoi(argv[1]) ;
   }
   
#include "setup_base.h"
   
   // some convenience stuff
   using namespace ecmech ;
#if KIN_TYPE
   typedef Kin_KMBalD_FFF Kinetics  ;
   typedef EvptnUpsdtProblem_FCC_B Prob ;
   typedef EvptnSolver_FCC_B Solver ;
#else
   typedef KineticsVocePL Kinetics ;
   typedef EvptnUpsdtProblem_FCC_A Prob ;
   typedef EvptnSolver_FCC_A Solver ;
#endif   
   
   typedef ecmech::SlipGeomFCC SlipGeom ;
   SlipGeom slipGeom ;

#if KIN_TYPE
   Kinetics kinetics(slipGeom.nslip) ;
#include "setup_kin_KMBalD_FFF.h"
#else
   Kinetics kinetics(slipGeom.nslip) ;
#include "setup_kin_VocePL.h"
#endif   

   typedef evptn::ThermoElastNCubic ThermoElastN ;
   ThermoElastN elastN ;
#include "setup_elastn.h"

   //////////////////////////////
   
   real8 p = 0.0, tK = 300.0 ;
   std::vector<real8> h_state_vec ;
   real8* h_state ;
   {
      std::vector<std::string> names ;
      std::vector<bool>        plot ;
      std::vector<bool>        state ;
      kinetics.getHistInfo( names, h_state_vec, plot, state ) ;
      h_state = &(h_state_vec[0]) ;
   }

#include "setup_conditions.h"   
   real8 detV = 1.0 ;
   real8 eVref = 0.0 ;
   real8 e_vecd_n[ecmech::ntvec] = {0.0} ;
   real8 Cn_quat[ecmech::qdim] = {1.0, 0.0, 0.0, 0.0} ;
   
   Prob prob(slipGeom, kinetics, elastN,
             dt,
             detV, eVref, p, tK,
             h_state, e_vecd_n, Cn_quat,
             d_vecd_sm, w_veccp_sm ) ;
   
   Solver solver(prob) ;

   snls::TrDeltaControl deltaControl ;
   deltaControl._deltaInit = 1e0 ;
   {
      int maxIter = 100 ;
      solver.setupSolver(maxIter, tolerance, &deltaControl, outputLevel) ;
   }
   
   real8* x = solver.getXPntr() ;
   for (int iX = 0; iX < prob.nDimSys; ++iX) {
      x[iX] = 0e0 ;
   }
   
   snls::SNLSStatus_t status = solver.solve( ) ;
   if ( status != snls::converged ){
      ECMECH_FAIL(__func__,"Solver failed to converge!");
   }
   std::cout << "Function evaluations: " << solver.getNFEvals() << std::endl ;

   std::cout << "Slip system shearing rates : " ;
   printVec<slipGeom.nslip>(prob.getGdot(), std::cout) ;

   //////////////////////////////////////////////////////////////////////
   //
   // this gives the same as above? does updateH internally, but with
   // beginning-of-step shearing rates which are all set to zero

   typedef ecmech::EosModelConst<false> EosModel ;
   EosModel eos ;
#include "setup_eos.h"
   
   real8 eInt[ecmech::ne] = {0.0} ;
   real8 stressSvecP[ecmech::nsvp] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                      0.0};
   static const int iHistLbGdot = evptn::NumHist<SlipGeom,Kinetics,ThermoElastN,EosModel>::iHistLbGdot ;
   static const int numHist     = evptn::NumHist<SlipGeom, Kinetics, ThermoElastN, EosModel>::numHist ;
   real8 hist[numHist] = {0.0} ;
   std::copy(Cn_quat, Cn_quat+ecmech::qdim, hist+evptn::iHistLbQ) ;
   std::copy(h_state, h_state+kinetics.nH,  hist+evptn::iHistLbH) ;
   real8* gdot = &(hist[iHistLbGdot]) ; // already zerod
   // do not bother with other stuff (like e_vecd_n) that is all zero above
   //
   real8 tkelv ;
   real8 sdd[ecmech::nsdd] ;
   real8 mtanSD[ecmech::nsvec2] ;
   //
   evptn::getResponseSngl<SlipGeom, Kinetics, ThermoElastN, EosModel>
      ( slipGeom, kinetics, elastN, eos,
        dt, 
        tolerance,
        d_svec_kk_sm, w_veccp_sm, volRatio,
        eInt, stressSvecP, hist,
        tkelv, sdd, mtanSD ) ;
   std::cout << "Function evaluations: " << hist[evptn::iHistA_nFEval] << std::endl ;
   
   std::cout << "Slip system shearing rates : " ;
   printVec<slipGeom.nslip>(gdot, std::cout) ;
      
}

