#include "SNLS_TrDLDenseG.h"

#include "ECMech_evptn.h"
#include "ECMech_kinetics.h"
#include "ECMech_slipgeom.h"

int main(int argc, char *argv[])
{
   // some convenience stuff
   using namespace ecmech ;
   typedef EvptnUpdstProblem< SlipGeomFCCA, KineticsVocePL, ThermoElastNCubic > EvptnUpsdtProblem_FCCAVocePl ; 
   
   ecmech::SlipGeomFCCA slipGeom ;

   ecmech::KineticsVocePL kinetics(slipGeom.nslip) ;
   {
      real8 mu = 1.0, xm = 0.1, gam_w = 1.0 ; 
      real8 params[kinetics.nParams] = { mu, xm, gam_w } ;
      kinetics.setParams( params ) ;
   }

   ecmech::ThermoElastNCubic elastN ;
   {
      real8 c11 = 300e-2, c12 = 100e-2, c44 = 100e-2 ;
      real8 params[elastN.nParams] = { c11, c12, c44 } ;
      kinetics.setParams( params ) ;
   }

   //////////////////////////////
   
   real8 p = 0.0, tK = 300.0 ;
   {
      real8 h_state[kinetics.nH] = { 100e-5 } ;
      kinetics.setVals( p, tK, h_state );
   }

   real8 dt = 1e-3 ;
   real8 detV = 1.0 ;
   real8 eVref = 0.0 ;
   real8 e_vecd_n[ecmech::ntvec] = {0.0} ;
   real8 Cn_quat[ecmech::qdim] = {1.0, 0.0, 0.0, 0.0} ;
   real8 d_vecd_sm[ecmech::ntvec] = {0.0, 1.0, 0.0, 0.0, 0.0} ;
   real8 w_veccp_sm[ecmech::nwvec] = {0.0, 0.0, 0.5} ;
   
   EvptnUpsdtProblem_FCCAVocePl prob(slipGeom, kinetics, elastN,
                                     dt,
                                     detV, eVref, p, tK,
                                     e_vecd_n, Cn_quat,
                                     d_vecd_sm, w_veccp_sm ) ;
   
   snls::SNLSTrDlDenseG<EvptnUpsdtProblem_FCCAVocePl> solver(prob) ;
   
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
   std::cout << "Function evaluations: " << solver.getNFEvals() << "\n";    
   
}