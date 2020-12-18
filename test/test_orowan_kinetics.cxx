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
#ifdef LARGE_DD
   const double gdotVal1 = 11387.933989173 - 1.145e-10;
#else
   const double gdotVal1 = 0.011441126690517;
#endif

   const double hUpdtTol = 1.0e-11;

   const int nslip = 12;
   const double init_tau = 1.0e4;
   double dt = 1e-2;

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
      double kin_vals[kinetics.nVals];
      double _hdn_scale = kinetics.getVals(kin_vals, 0.0, 300.0, &(init[0]));

      double gdot[nslip] = {0.0};
      double dgdot_dtau[nslip] = {0.0};
      double dgdot_dg[nslip] = {0.0};
      double taua[nslip] = {0.0};
      for (int ig = 0; ig < nslip; ig++) {
         taua[ig] = init_tau; 
      }
      kinetics.evalGdots(gdot, dgdot_dtau, dgdot_dg, taua, kin_vals);

#ifdef ECMECH_DEBUG
      std::cout << "Gdot values : ";
      printVec<12>(gdot, std::cout);
#endif
      // Our numbers are pretty large here, we should do a relative tolerance instead
      EXPECT_LT(fabs(gdot[0] - gdotVal1), hUpdtTol) << "Did not get expected value";
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
