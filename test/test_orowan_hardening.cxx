#include <gtest/gtest.h>

#include "SNLS_TrDLDenseG.h"

#include "ECMech_kinetics.h"
#include "ECMech_cases.h"

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>

static int outputLevel = 1;

TEST(ecmech, hard_orowan_fcc)
{
   using namespace ecmech;
   const double hUpdtVal1 = 1.6211613943304e+16;
   const double hUpdtVal2 = 4.1497322905107e+17;

   const double hUpdtTol = 1e-11;

   const int nslip = 12;
   double dt = 1e-2;
   double gdot[nslip] = { 0.1 };

   {

   Kin_OroD_Iso_FCC kinetics(nslip);
#include "setup_kin_OroD_Iso_FCC.h"

      std::vector<double>       init;
      {
         std::vector<std::string> names;
         std::vector<bool>        plot;
         std::vector<bool>        state;
         kinetics.getHistInfo(names, init, plot, state);
      }
      double hs_u[kinetics.nH];
      int nFEvals = kinetics.updateH(hs_u, &(init[0]), dt, gdot, outputLevel);
      std::cout << "Converged with nFEvals : " << nFEvals << std::endl;

      EXPECT_TRUE(nFEvals == 4) << "Not the expected number of function evaluations";
#ifdef ECMECH_DEBUG
      std::cout << "Updated hardness state : ";
      printVec<kinetics.nH>(hs_u, std::cout);
#endif
      // Our numbers are pretty large here, we should do a relative tolerance instead
      EXPECT_LT(fabs((hs_u[0] - hUpdtVal1)/hUpdtVal1), hUpdtTol) << "Did not get expected value";
      EXPECT_LT(fabs((hs_u[12] - hUpdtVal2)/hUpdtVal2), hUpdtTol) << "Did not get expected value";
   }
}

int main(int argc, char *argv[])
{
   ::testing::InitGoogleTest(&argc, argv);
   if (argc > 1) {
      outputLevel = atoi(argv[1]);
   }
   std::cout << "got outputLevel : " << outputLevel << std::endl;

   return RUN_ALL_TESTS();
}
