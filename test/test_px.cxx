#include "SNLS_TrDLDenseG.h"

#include "ECMech_cases.h"
#include "ECMech_evptnWrap.h"

#include <random>

#define STACK_PARAMS

#ifndef KIN_TYPE
#define KIN_TYPE 1
#endif

int main(int argc, char *argv[])
{
   
   // can adjust these to change the computational workload
   //
   const int nPassed = 16 ; 
   const int nStep = 100 ;
   //
   const real8 weight = 1.0 / (real8)(nPassed) ;

   int outputLevel = 0 ;
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

   mmodel->setOutputLevel( outputLevel ) ; // would not normally do this in a production setting
   
   std::vector<real8> histInit_vec ;
   {
      std::vector<std::string> names ;
      std::vector<bool>        plot ;
      std::vector<bool>        state ;
      mmb->getHistInfo(names, histInit_vec, plot, state ) ;
   }

   //
   // set up hist and other state information
   //
   const int numHist = mmodel->numHist ;
   real8 V_hist[numHist * nPassed] ;
   {
      std::default_random_engine gen ;
      std::normal_distribution<real8> distrib(0.0,1.0) ;
      for ( int iPassed = 0; iPassed < nPassed; iPassed++ ) {
         real8* hist = &(V_hist[numHist*iPassed]) ;
         for ( int iHist = 0; iHist < numHist; iHist++ ) {
            hist[iHist] = histInit_vec[iHist];
         }
         real8* q_state = &(hist[ecmech::evptn::iHistLbQ]) ;
         q_state[0] = distrib(gen) ;
         q_state[1] = distrib(gen) ;
         q_state[2] = distrib(gen) ;
         q_state[3] = distrib(gen) ;
         vecsVNormalize<ecmech::qdim>(q_state) ;
      }
   }

   real8 relRate = 1e-6 ;
   real8 dt = 0.002/relRate ;

   real8 d_svec_kk_sm[ecmech::nsvp] = {-0.5*relRate, -0.5*relRate, 1.0*relRate,
                                       0.0, 0.0, 0.0,
                                       0.0} ;
   real8 V_d_svec_kk_sm[ecmech::nsvp*nPassed] ;
   for ( int iPassed = 0; iPassed < nPassed; iPassed++ ) {
      int pOffsetSVP = ecmech::nsvp * iPassed;
      std::copy(d_svec_kk_sm, d_svec_kk_sm+ecmech::nsvp, &(V_d_svec_kk_sm[pOffsetSVP])) ;
   }

   std::vector<real8> w_veccp_sm_vec(ecmech::nwvec*nPassed, 0.0) ;
   real8* V_w_veccp_sm = &(w_veccp_sm_vec[0]);

   std::vector<real8> volRatio_vec(ecmech::nvr*nPassed, 1.0) ; // not really 1 for all entries, but this works given what happens below
   real8* V_volRatio = &(volRatio_vec[0]) ;

   std::vector<real8> eInt_vec(ecmech::ne*nPassed, 0.0) ;
   real8* V_eInt = &(eInt_vec[0]) ;
   
   std::vector<real8> stressSvecP_vec(ecmech::nsvp*nPassed, 0.0) ;
   real8* V_stressSvecP = &(stressSvecP_vec[0]) ;
   
   real8 V_tkelv[nPassed] ;
   real8 V_sdd[ecmech::nsdd*nPassed] ;
   //
   std::cout << "# time, Axial deviatoric stress : " << std::endl ;
   real8 time = 0.0 ;
   //
   for ( int iStep=0; iStep<nStep; ++iStep ) {
      //
      time += dt ;

      // update current relative volume from the volumetric deformation rate
      //
      for ( int iPassed = 0; iPassed < nPassed; iPassed++ ) {
         int pOffsetSVP = ecmech::nsvp * iPassed;
         int pOffsetVR  = ecmech::nvr  * iPassed;
         V_volRatio[0 + pOffsetVR] = V_volRatio[1 + pOffsetVR];
         V_volRatio[1 + pOffsetVR] = V_volRatio[0 + pOffsetVR] * exp( V_d_svec_kk_sm[ecmech::iSvecP + pOffsetSVP] * dt );
         V_volRatio[3 + pOffsetVR] = V_volRatio[1 + pOffsetVR] - V_volRatio[0 + pOffsetVR] ;
         V_volRatio[2 + pOffsetVR] = V_volRatio[3 + pOffsetVR] / ( dt * 0.5*(V_volRatio[0+pOffsetVR]+V_volRatio[1+pOffsetVR]) )  ;
      }
      
      mmb->getResponse( dt,
                        V_d_svec_kk_sm, V_w_veccp_sm, V_volRatio,
                        V_eInt, V_stressSvecP, V_hist, V_tkelv, V_sdd, nullptr,
                        nPassed ) ;

      real8 sAvg = 0.0 ;
      for ( int iPassed=0; iPassed<nPassed; ++iPassed ) {
         sAvg += weight * V_stressSvecP[iPassed*ecmech::nsvp+2] ;
      }
      std::cout << time << " "
                << std::setprecision(14) << sAvg
                << std::endl ;
      
   }
      
   delete mmodel ;

   exit(0) ;
}

