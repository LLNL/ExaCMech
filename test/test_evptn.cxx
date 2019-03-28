#include "SNLS_TrDLDenseG.h"

#include "ECMech_evptn.h"
#include "ECMech_cases.h"
#include "ECMech_kinetics.h"
#include "ECMech_slipgeom.h"
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
   
   // some convenience stuff
   using namespace ecmech ;
#if KIN_TYPE
   typedef Kin_KMBalD_FFF KinType  ;
   typedef EvptnUpsdtProblem_FCC_B Prob ;
   typedef EvptnSolver_FCC_B Solver ;
#else
   typedef KineticsVocePL KinType ;
   typedef EvptnUpsdtProblem_FCC_A Prob ;
   typedef EvptnSolver_FCC_A Solver ;
#endif   
   
   ecmech::SlipGeomFCC slipGeom ;

#if KIN_TYPE
   KinType kinetics(slipGeom.nslip) ;
#include "setup_kin_KMBalD_FFF.h"
#else
   KinType kinetics(slipGeom.nslip) ;
#include "setup_kin_VocePL.h"
#endif   

   ecmech::ThermoElastNCubic elastN ;
   {
      real8 c11 = 300e-2, c12 = 100e-2, c44 = 100e-2 ;
      real8 params[elastN.nParams] = { c11, c12, c44 } ;
      elastN.setParams( params ) ;
   }

   //////////////////////////////
   
   real8 p = 0.0, tK = 300.0 ;
   real8 h_state[kinetics.nH] = { 100e-5 } ;

   real8 dt = 1e-1 ;
   real8 detV = 1.0 ;
   real8 eVref = 0.0 ;
   real8 e_vecd_n[ecmech::ntvec] = {0.0} ;
   real8 Cn_quat[ecmech::qdim] = {1.0, 0.0, 0.0, 0.0} ;
   real8 d_vecd_sm[ecmech::ntvec] = {0.0, 1.0, 0.0, 0.0, 0.0} ;
   real8 w_veccp_sm[ecmech::nwvec] = {0.0, 0.0, 0.5} ;
   
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
      real8 tolerance = 1e-10 ;
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
   
}
