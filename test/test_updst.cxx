#include <gtest/gtest.h>

#include "SNLS_TrDLDenseG.h"

#include "ECMech_cases.h"
#include "ECMech_evptnWrap.h"

#define STACK_PARAMS

#ifndef DO_FD_CHECK_MTAN
#define DO_FD_CHECK_MTAN 0
#endif

#if DO_FD_CHECK_MTAN
// if doing finite-difference check of tangent stiffness, then change some other defaults
#define NON_I_QUAT 1
#define KIN_TYPE 0
#define XM_MUSHY 1
#endif

#ifndef KIN_TYPE
#define KIN_TYPE 1
#endif

#ifndef NON_I_QUAT
#define NON_I_QUAT 0
#endif

static int outputLevel = 0 ;

#include "test_expectedVals.h"

TEST(ecmech, updst_a)
{
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
#if NON_I_QUAT
   real8* q_state = &(hist[ecmech::evptn::iHistLbQ]) ;
   {
      real8 th = 0.2 ;
      real8 et = 0.7 ;
      q_state[0] = cos(0.5*th);
      q_state[1] = sin(0.5*th)*cos(et) ;
      q_state[2] = sin(0.5*th)*sin(et) ;
      q_state[3] = 0.0 ;
   }
#endif
   
   // int numHist = hist_vec.size() ; // should equal mmodel->numHist

   const int nPassed = 1 ; // just do a single point here as a simple example

   mmodel->setOutputLevel( outputLevel ) ; // would not normally do this in a production setting

   static const int iHistLbGdot = mmodel->iHistLbGdot ;
   real8* gdot = &(hist[iHistLbGdot]) ; 
   
   std::cout << "Initial hist : " ;
   printVec(hist, mmodel->numHist, std::cout) ;

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
   printVec(hist, mmodel->numHist, std::cout) ;
   
   std::cout << "Hardness state : " ;
   printVec<mmodel->nH>(&(hist[ecmech::evptn::iHistLbH]), std::cout) ;
   
   std::cout << "Slip system shearing rates : " ;
   printVec<mmodel->nslip>(gdot, std::cout) ;

   EXPECT_TRUE( hist[evptn::iHistA_nFEval] == expectedNFEvals ) << "Not the expected number of function evaluations" ;
   EXPECT_LT( fabs( hist[evptn::iHistLbE+1] - expectedE2 ) , 1e-10 ) <<
      "Did not get expected value for lattice strain component" ;
   EXPECT_LT( fabs( hist[evptn::iHistLbQ] - expectedQ1 ) , 1e-8 ) <<
      "Did not get expected value for quat_1" ;
   EXPECT_LT( fabs( gdot[1] - expectedGdotVal ) , 1e-8 ) <<
      "Did not get expected value for gdot[1]" ;
   
   delete mmodel ;
}

TEST(ecmech, driver_a)
{
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
#if NON_I_QUAT
   real8* q_state = &(hist[ecmech::evptn::iHistLbQ]) ;
   {
      real8 th = 0.2 ;
      real8 et = 0.7 ;
      q_state[0] = cos(0.5*th);
      q_state[1] = sin(0.5*th)*cos(et) ;
      q_state[2] = sin(0.5*th)*sin(et) ;
      q_state[3] = 0.0 ;
   }
#endif
   
   // int numHist = hist_vec.size() ; // should equal mmodel->numHist

   const int nPassed = 1 ; // just do a single point here as a simple example

   mmodel->setOutputLevel( outputLevel ) ; // would not normally do this in a production setting
   
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
   
#if !(DO_FD_CHECK_MTAN)
   std::cout << "# time, Axial deviatoric stress, h[0], p : " << std::endl ;
#endif
   real8 time = 0.0 ;
   //
   for ( int iStep=0; iStep<nStep; ++iStep ) {
      //
      time += dt ;

      // update current relative volume from the volumetric deformation rate
      //
      volRatio[0] = volRatio[1] ;
      volRatio[1] = volRatio[0] * exp( d_svec_kk_sm[ecmech::iSvecP] * dt ) ;
      volRatio[3] = volRatio[1] - volRatio[0] ;
      volRatio[2] = volRatio[3] / ( dt * 0.5 *(volRatio[0]+volRatio[1]) )  ;
      
      mmb->getResponse( dt, d_svec_kk_sm, w_veccp_sm, volRatio,
                        eInt, stressSvecP, hist, tkelv, sdd, nullptr,
                        nPassed ) ;

#if !(DO_FD_CHECK_MTAN)
      std::cout << time << " "
                << std::setprecision(14) << stressSvecP[2] << " "
                << std::setprecision(14) << hist[ecmech::evptn::iHistLbH+0] << " "
                << std::setprecision(14) << stressSvecP[iSvecP] << " "
                << std::endl ;

      // std::cout << "hist : " ;
      // printVec(hist, numHist, std::cout) ;
#endif
      
   }
#if KIN_TYPE && !(DO_FD_CHECK_MTAN)
   EXPECT_LT( fabs( stressSvecP[2] - 0.006664661118275 ) , 1e-10 ) <<
      "Did not get expected value for stress component" ;
   EXPECT_LT( fabs( hist[ecmech::evptn::iHistLbH+0] - 88.61845050083 ) , 1e-8 ) <<
      "Did not get expected value for history variable" ;
   EXPECT_LT( fabs( stressSvecP[iSvecP] -  0.00332602112947 ) , 1e-10 ) <<
      "Did not get expected value for stress component" ;
#endif

#if DO_FD_CHECK_MTAN
   {
      //
      // do another step, and do finite differencing to check mtanSD

      std::vector<real8> hist_ref(hist, hist+mmodel->numHist) ;
      std::vector<real8> eInt_ref(eInt, eInt+ecmech::ne) ;
      std::vector<real8> stressSvecP_ref(stressSvecP, stressSvecP+ecmech::nsvp) ;
      real8 v_ref = volRatio[1] ;

      volRatio[0] = v_ref ;
      volRatio[1] = volRatio[0] * exp( d_svec_kk_sm[ecmech::iSvecP] * dt ) ;
      volRatio[3] = volRatio[1] - volRatio[0] ;
      volRatio[2] = volRatio[3] / ( dt * 0.5 *(volRatio[0]+volRatio[1]) )  ;
      
      real8 mtanSD_an[ecmech::nsvec2] ;
      mmb->getResponse( dt, d_svec_kk_sm, w_veccp_sm, volRatio,
                        eInt, stressSvecP, hist, tkelv, sdd, mtanSD_an, 
                        nPassed ) ;


      real8 stressSvec[ecmech::nsvec] ;
      svecpToSvec( stressSvec, stressSvecP ) ;
      
      std::cout << "mtanSD_an : " << std::endl ;
      printMat<ecmech::nsvec>( mtanSD_an, std::cout ) ;
      
      real8 d_svec_kk_sm_pert[ecmech::nsvp] ;
      const real8 pertVal = 1e-8*relRate ;
      real8 mtanSD_fd[ecmech::nsvec2] ;
      //
      real8 eInt_pert[ecmech::ne] ;
      real8 stressSvecP_pert[ecmech::nsvp] ;
      //
      for ( int jSvec=0; jSvec<ecmech::nsvec; ++jSvec ) {
   
         std::copy(d_svec_kk_sm, d_svec_kk_sm+ecmech::nsvp, d_svec_kk_sm_pert) ;
         if ( jSvec < 3 ) {
            d_svec_kk_sm_pert[jSvec] += pertVal ;
            real8 d_kk = d_svec_kk_sm_pert[0] + d_svec_kk_sm_pert[1] + d_svec_kk_sm_pert[2] ;
            d_svec_kk_sm_pert[ecmech::iSvecP] += d_kk ;
            d_svec_kk_sm_pert[0] += (- ecmech::onethird * d_kk) ;
            d_svec_kk_sm_pert[1] += (- ecmech::onethird * d_kk) ;
            d_svec_kk_sm_pert[2] += (- ecmech::onethird * d_kk) ;
         }
         else {
            // factor of 2 to go with l_ddsdde_gamma being true in call to mtan_conv_sd_svec ;
            d_svec_kk_sm_pert[jSvec] += 0.5 * pertVal ;
         }
         //
         volRatio[0] = v_ref ;
         volRatio[1] = volRatio[0] * exp( d_svec_kk_sm_pert[ecmech::iSvecP] * dt ) ;
         volRatio[3] = volRatio[1] - volRatio[0] ;
         volRatio[2] = volRatio[3] / ( dt * 0.5 *(volRatio[0]+volRatio[1]) )  ;

         std::copy(eInt_ref.begin(), eInt_ref.end(), eInt_pert) ;
         std::copy(stressSvecP_ref.begin(), stressSvecP_ref.end(), stressSvecP_pert) ;
         
         real8 tkelv_pert[nPassed] ;
         real8 sdd_pert[ecmech::nsdd*nPassed] ;
   
         std::copy(hist_ref.begin(), hist_ref.end(), hist) ; // make hist equal to hist_ref again
         
         mmb->getResponse( dt, d_svec_kk_sm_pert, w_veccp_sm, volRatio,
                           eInt_pert, stressSvecP_pert, hist, tkelv_pert, sdd_pert, nullptr,
                           nPassed ) ;
   
         real8 stressSvec_pert[ecmech::nsvec] ;
         svecpToSvec( stressSvec_pert, stressSvecP_pert ) ;
         //
         for ( int iSvec=0; iSvec<ecmech::nsvec; ++iSvec ) {
            // divide by dt because tangent gets converted to a per-strain-increment type quantity
            mtanSD_fd[ECMECH_NN_INDX(iSvec,jSvec,ecmech::nsvec)] = ( stressSvec_pert[iSvec] - stressSvec[iSvec] ) / pertVal / dt ;
         }
         
      }
      std::cout << "mtanSD_fd : " << std::endl ;
      printMat<ecmech::nsvec>( mtanSD_fd, std::cout ) ;

      // do not bother restoring things to evaluation at non-perturbed condition

      for ( int iiMtan=0; iiMtan<ecmech::nsvec2; ++iiMtan ) {
         EXPECT_LT( fabs( mtanSD_fd[iiMtan] - mtanSD_an[iiMtan] ) , 1e-3 ) <<
            "Analytic and finite-differenced mtan differ by more than expected" ;
      }
      
   }
#endif // DO_FD_CHECK_MTAN

   delete mmodel ;
}

int main(int argc, char *argv[])
{
   ::testing::InitGoogleTest(&argc, argv);
   if ( argc > 1 ) {
      outputLevel = atoi(argv[1]) ;
   }
   std::cout << "got outputLevel : " << outputLevel << std::endl ;

  return RUN_ALL_TESTS();
}
