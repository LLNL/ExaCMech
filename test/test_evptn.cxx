#include "SNLS_TrDLDenseG.h"

#include "ECMech_evptn.h"
#include "ECMech_evptn_cases.h"
#include "ECMech_kinetics.h"
#include "ECMech_slipgeom.h"
#include "ECMech_util.h"

#define KIN_KMBAL 1

int main(int , // argc,
         char ** // *argv[]
         )
{
   // some convenience stuff
   using namespace ecmech ;
#if KIN_KMBAL
   typedef Kin_KMBalD_FFF KinType  ;
   typedef EvptnUpsdtProblem_FCC_B Prob ;
   typedef EvptnSolver_FCC_B Solver ;
#else
   typedef KineticsVocePL KinType ;
   typedef EvptnUpsdtProblem_FCC_A Prob ;
   typedef EvptnSolver_FCC_A Solver ;
#endif   
   
   ecmech::SlipGeomFCC slipGeom ;

#if KIN_KMBAL
   KinType kinetics(slipGeom.nslip) ;
   {
      real8
         mu     = 1.0,
         tK_ref = 300.,
         c_1    = 20000.,
         tau_a  = 0.004,
         p      = 0.28,
         q      = 1.34,
         gam_wo = 20.,
         gam_ro = 1e3,
         wrD    = 0.02,
         go     = 10e-5,
         s      = 5e-5 ;
      real8 params[kinetics.nParams] = { mu, tK_ref, c_1, tau_a, p, q, gam_wo, gam_ro, wrD, go, s } ;
      kinetics.setParams( params ) ;
   }
#else
   KinType kinetics(slipGeom.nslip) ;
   {
      real8 mu = 1.0, xm = 0.01, gam_w = 1.0 ; 
      real8 params[kinetics.nParams] = { mu, xm, gam_w } ;
      kinetics.setParams( params ) ;
   }
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
      int outputLevel = 10 ;
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
