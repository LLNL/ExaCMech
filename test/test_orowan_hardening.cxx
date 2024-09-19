#include <gtest/gtest.h>

#include "SNLS_TrDLDenseG.h"

#include "kinetics/ECMech_kinetics.h"
#include "cases/ECMech_cases_fcc_defs.h"

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>

static int outputLevel = 1;

TEST(ecmech, hard_orowan_fcc)
{
   using namespace ecmech;
#ifdef LARGE_DD
   const double hUpdtVal1 = 1.000688247427146e+4;
   const double hUpdtVal2 = 4.001148109947428e+4;
   const int nevals = 2;
#else
   const double hUpdtVal1 = 1.992139383887259e-2;
   const double hUpdtVal2 = 5.896582642644239e-2;
   const int nevals = 4;
#endif

   const double hUpdtTol = 1e-11;

   const int nslip = 12;
   double dt = 0.001;
   double gdot[nslip] = { 1.0 };

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
      double hvals[12] = { 0.0 };
      double tkelv = 300;
      int nFEvals = kinetics.updateH(hs_u, &(init[0]), dt, gdot, hvals, tkelv, outputLevel);
      std::cout << "Converged with nFEvals : " << nFEvals << std::endl;

      EXPECT_TRUE(nFEvals == nevals) << "Not the expected number of function evaluations";
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
