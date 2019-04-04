#include "SNLS_TrDLDenseG.h"

#include "ECMech_cases.h"
#include "ECMech_evptnWrap.h"

#define STACK_PARAMS

int main(int argc, char *argv[])
{
   int outputLevel = 1 ;
   if ( argc > 1 ) {
      outputLevel = atoi(argv[1]) ;
   }

   using namespace ecmech ;

   matModelEvptn_FCC_B* mmodel = new matModelEvptn_FCC_B() ;
   matModelBase* mmb = dynamic_cast<matModelBase*>(mmodel) ;

#include "setup_base.h"
   std::vector<int>           opts ; // none
   std::vector<std::string>   strs ; // none
   std::vector<real8>         params{ rho0, cvav, tolerance } ;
#include "setup_elastn.h"
#include "setup_kin_KMBalD_FFF.h"
#include "setup_eos.h"
   //
   mmb->initFromParams( opts, params, strs ) ;
   //
   mmb->complete() ;

   std::vector<real8>       hist_vec ;
   {
      std::vector<std::string> names ;
      std::vector<bool>        plot ;
      std::vector<bool>        state ;
      mmb->getHistInfo(names, hist_vec, plot, state ) ;
   }
   int numHist = hist_vec.size() ; // should equal mmodel->numHist
   real8* hist = &(hist_vec[0]);
   static const int iHistLbGdot = mmodel->iHistLbGdot ;
   real8* gdot = &(hist[iHistLbGdot]) ; 
   real8* h_state = &(hist[ecmech::evptn::iHistLbH]) ; 
   
   std::cout << "Initial hist : " ;
   printVec(hist, numHist, std::cout) ;
      
   const int nPassed = 1 ; // just do a single point here as a simple example
   //
   mmodel->setOutputLevel( outputLevel ) ; // would not normally do this in a production setting
   //
   {
#include "setup_conditions.h"      
      real8 eInt[ecmech::ne] = {0.0} ;
      real8 stressSvecP[ecmech::nsvp] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                         0.0};
      real8 tkelv[nPassed] ;
      real8 sdd[ecmech::nsdd*nPassed] ;
      mmb->getResponse( dt, d_svec_kk_sm, w_veccp_sm, volRatio,
                        eInt, stressSvecP, hist, tkelv, sdd, nullptr,
                        nPassed ) ;

      std::cout << "Function evaluations: " << hist[evptn::iHistA_nFEval] << std::endl ;
   
   }
      
   std::cout << "Updated hist : " ;
   printVec(hist, numHist, std::cout) ;
      
   std::cout << "Hardness state : " ;
   printVec<mmodel->nH>(h_state, std::cout) ;
      
   std::cout << "Slip system shearing rates : " ;
   printVec<mmodel->nslip>(gdot, std::cout) ;
   
   delete mmodel ;

   exit(0) ;
}

