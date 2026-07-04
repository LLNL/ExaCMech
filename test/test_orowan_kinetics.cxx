/**
 * @file test_orowan_kinetics.cxx
 *
 * @brief Isolated unit test of `KineticsOrowanD::evalGdots` (the slip-rate law only --
 * no hardening-state ODE, contrast `test_orowan_hardening.cxx`): builds
 * `Kin_OroD_Iso_FCC` directly (not via a full `matModelBase`), evaluates the
 * dislocation-density-based CRSS/reference-rate values (`getVals`) from the initial
 * hardening state, then evaluates the resulting slip rate on every slip system for a
 * fixed, uniform resolved shear stress (`init_tau`) and checks slip system 0's slip
 * rate against a recorded reference value.
 *
 * `LARGE_DD` (see `test/CMakeLists.txt`) swaps in the `qM`/`qT` (initial mobile/total
 * dislocation density) values from `setup_kin_OroD_Iso_FCC.h`'s `LARGE_DD` branch,
 * giving a much higher initial CRSS and hence a much smaller slip rate at the same
 * resolved shear stress -- exercising the same code path at a very different point in
 * its range.
 */

#include <gtest/gtest.h>

#include "SNLS_TrDLDenseG.h"

#include "kinetics/ECMech_kinetics.h"
#include "cases/ECMech_cases_fcc_defs.h"

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>

static int outputLevel = 1;

TEST(ecmech, kin_orowan_fcc)
{
   using namespace ecmech;
#ifdef LARGE_DD
   const double gdotVal1 = 4.87738409778574465e-10;
#else
   const double gdotVal1 = 2.59181779318282111e+2;
#endif

   const double hUpdtTol = 1.0e-11;

   constexpr int nslip = 12;
   const double init_tau = 1.0e-2;

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
      double kin_vals[Kin_OroD_Iso_FCC::nVals];
      kinetics.getVals(kin_vals, 0.0, 300.0, &(init[0]));

      double gdot[nslip] = {0.0};
      double dgdot_dtau[nslip] = {0.0};
      double taua[nslip] = {0.0};
      for (int ig = 0; ig < nslip; ig++) {
         taua[ig] = init_tau; 
      }
      kinetics.evalGdots(gdot, dgdot_dtau, taua, kin_vals);

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
