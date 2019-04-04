#include "SNLS_TrDLDenseG.h"

#include "ECMech_cases.h"
#include "ECMech_evptnWrap.h"

#define STACK_PARAMS

#ifndef KIN_TYPE
#define KIN_TYPE 1
#endif

#ifndef AS_DRIVER
#define AS_DRIVER 0
#endif

int main(int argc, char *argv[])
{
#if AS_DRIVER
   int outputLevel = 0 ;
#else
   int outputLevel = 1 ;
#endif   
   if ( argc > 1 ) {
      outputLevel = atoi(argv[1]) ;
   }

   using namespace ecmech ;

#if KIN_TYPE
   matModelEvptn_FCC_B* mmodel = new matModelEvptn_FCC_B() ;
#else
   matModelEvptn_FCC_A* mmodel = new matModelEvptn_FCC_A() ;
#endif
   matModelBase* mmb = dynamic_cast<matModelBase*>(mmodel) ;
   
#include "setup_base.h"
   std::vector<int>           opts ; // none
   std::vector<std::string>   strs ; // none
   std::vector<real8>         params{ rho0, cvav, tolerance } ;
#if KIN_TYPE

#include "setup_elastn.h"
#include "setup_kin_KMBalD_FFF.h"
#include "setup_eos.h"

#else

#include "setup_elastn.h"
#include "setup_kin_VocePL.h"
#include "setup_eos.h"
   
#endif
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
   real8* hist = &(hist_vec[0]);
   
   real8* h_state = &(hist[ecmech::evptn::iHistLbH]) ; 

   const int nPassed = 1 ; // just do a single point here as a simple example

   mmodel->setOutputLevel( outputLevel ) ; // would not normally do this in a production setting
   
#if AS_DRIVER
   real8 relRate = 1e-6 ;
   real8 d_svec_kk_sm[ecmech::nsvp] = {-0.5*relRate, -0.5*relRate, 1.0*relRate,
                                       0.0, 0.0, 0.0,
                                       0.0} ;
   // vecsVsa<ecmech::nsvp>(d_svec_kk_sm, sqr2b3) ; // nope, choose not to do that here
   //
   real8 d_vecd_sm[ecmech::ntvec] ;
   svecToVecd(d_vecd_sm, d_svec_kk_sm) ;

   // dt value in setup_conditions.h is meant to stress the implementation --
   // here go with a smaller value to be able to make a nicer curve
   real8 dt = 0.002/relRate ;
   int nStep = 100 ;

   real8 w_veccp_sm[ecmech::nwvec] = {0.0, 0.0, 0.0} ;

   real8 volRatio[ecmech::nvr] = {1.0, 1.0, 0.0, 0.0} ;

   real8 eInt[ecmech::ne] = {0.0} ;
   real8 stressSvecP[ecmech::nsvp] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                      0.0};
   real8 tkelv[nPassed] ;
   real8 sdd[ecmech::nsdd*nPassed] ;
   
   std::cout << "# time, Axial deviatoric stress, h[0], p : " << std::endl ;
   real8 time = 0.0 ;
   //
   for ( int iStep=0; iStep<nStep; ++iStep ) {
      //
      time += dt ;

      // update current relative volume from the volumetric deformation rate
      //
      volRatio[0] = volRatio[1] ;
      volRatio[1] = volRatio[1] * exp( d_svec_kk_sm[ecmech::iSvecS] * dt ) ;
      volRatio[3] = volRatio[1] - volRatio[0] ;
      volRatio[2] = volRatio[3] / ( dt * 0.5 *(volRatio[0]+volRatio[1]) )  ;
      
      mmb->getResponse( dt, d_svec_kk_sm, w_veccp_sm, volRatio,
                        eInt, stressSvecP, hist, tkelv, sdd, nullptr,
                        nPassed ) ;

      std::cout << time << " "
                << std::setprecision(14) << stressSvecP[2] << " "
                << std::setprecision(14) << h_state[0] << " "
                << std::setprecision(14) << stressSvecP[iSvecP] << " "
                << std::endl ;
   }
      
#else

   int numHist = hist_vec.size() ; // should equal mmodel->numHist
   
   static const int iHistLbGdot = mmodel->iHistLbGdot ;
   real8* gdot = &(hist[iHistLbGdot]) ; 
   
   std::cout << "Initial hist : " ;
   printVec(hist, numHist, std::cout) ;

#include "setup_conditions.h"
   
   {
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

#endif
   
   delete mmodel ;

   exit(0) ;
}

