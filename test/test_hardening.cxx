#include <gtest/gtest.h>

#include "SNLS_TrDLDenseG.h"

#include "ECMech_kinetics.h"
#include "ECMech_cases.h"

static int outputLevel = 1 ;

TEST(ecmech,hard_voce_a)
{
   using namespace ecmech ;

   const real8 hUpdtVal = 0.001016620868315 ;
   const real8 hUpdtTol = 1e-11 ;
   
   const int nslip = 12 ;
   real8 dt = 1e-1 ;
   real8 gdot[nslip] = { 1.0/nslip } ;
   
   {
      KineticsVocePL kinetics(nslip) ;
#include "setup_kin_VocePL.h"
   
      std::vector<real8>       init ;
      {
         std::vector<std::string> names ;
         std::vector<bool>        plot ;
         std::vector<bool>        state ;
         kinetics.getHistInfo( names, init, plot, state ) ;
      }
      real8 hs_u[kinetics.nH] ;
      int nFEvals = kinetics.updateH( hs_u, &(init[0]), dt, gdot, outputLevel ) ;
      std::cout << "Converged with nFEvals : " << nFEvals << std::endl ;
      EXPECT_TRUE( nFEvals == 2 ) << "Not the expected number of function evaluations" ;
   
      std::cout << "Updated hardness state : " ;
      printVec<kinetics.nH>(hs_u, std::cout) ;
      EXPECT_LT( fabs(hs_u[0]-hUpdtVal) , hUpdtTol ) << "Did not get expected value" ;
   }
}

TEST(ecmech,hard_kmbaldfff_a)
{
   using namespace ecmech ;

   const real8 hUpdtVal = 0.6633659171982 ;
   const real8 hUpdtTol = 1e-8 ;

   const int nslip = 12 ;
   real8 dt = 1e-1 ;
   real8 gdot[nslip] = { 1.0/nslip } ;
   
   {
      Kin_KMBalD_FFF kinetics(nslip) ;
#include "setup_kin_KMBalD_FFF.h"
   
      std::vector<real8>       init ;
      {
         std::vector<std::string> names ;
         std::vector<bool>        plot ;
         std::vector<bool>        state ;
         kinetics.getHistInfo( names, init, plot, state ) ;
      }
      real8 hs_u[kinetics.nH] ;
      int nFEvals = kinetics.updateH( hs_u, &(init[0]), dt, gdot, outputLevel ) ;
      std::cout << "Converged with nFEvals : " << nFEvals << std::endl ;
      EXPECT_TRUE( nFEvals == 5 ) << "Not the expected number of function evaluations" ;
      
      std::cout << "Updated hardness state : " ;
      printVec<kinetics.nH>(hs_u, std::cout) ;
      EXPECT_LT( fabs(hs_u[0]-hUpdtVal) , hUpdtTol ) << "Did not get expected value" ;
   }
   
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
