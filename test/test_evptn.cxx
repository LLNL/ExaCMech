#include <gtest/gtest.h>

#include "SNLS_TrDLDenseG.h"

#include "ECMech_evptn.h"
#include "ECMech_cases.h"
#include "ECMech_kinetics.h"
#include "ECMech_slipgeom.h"
#include "ECMech_eosSimple.h"
#include "ECMech_util.h"

#ifndef KIN_TYPE
#define KIN_TYPE 1
#endif

static int outputLevel = 1;

#include "test_expectedVals.h"

TEST(ecmech, evptn_a)
{
#include "setup_base.h"

   // some convenience stuff
   using namespace ecmech;

#if KIN_TYPE == 3
   typedef ecmech::SlipGeom_BCC_A SlipGeom;
   typedef Kin_BCC_A Kinetics;
   typedef EvptnUpsdtProblem_BCC_A Prob;
   typedef EvptnSolver_BCC_A Solver;
   typedef evptn::ThermoElastNCubic ThermoElastN;
#elif KIN_TYPE == 2
   typedef ecmech::SlipGeom_HCP_A SlipGeom;
   typedef Kin_HCP_A Kinetics;
   typedef EvptnUpsdtProblem_HCP_A Prob;
   typedef EvptnSolver_HCP_A Solver;
   typedef evptn::ThermoElastNHexag ThermoElastN;
#elif KIN_TYPE == 1
   typedef ecmech::SlipGeomFCC SlipGeom;
   typedef Kin_FCC_B Kinetics;
   typedef EvptnUpsdtProblem_FCC_B Prob;
   typedef EvptnSolver_FCC_B Solver;
   typedef evptn::ThermoElastNCubic ThermoElastN;
#else
   typedef ecmech::SlipGeomFCC SlipGeom;
   typedef Kin_FCC_A Kinetics;
   typedef EvptnUpsdtProblem_FCC_A Prob;
   typedef EvptnSolver_FCC_A Solver;
   typedef evptn::ThermoElastNCubic ThermoElastN;
#endif

   SlipGeom slipGeom;
   Kinetics kinetics(slipGeom.nslip);
   ThermoElastN elastN;

#if KIN_TYPE == 3
#include "setup_slipGeom.h"
#include "setup_kin_KMBalD_TFF_BCC_A.h"
#include "setup_elastn.h"
   const int iGdotExpected = 1;
#elif KIN_TYPE == 2
#include "setup_slipGeom_HCP.h"
#include "setup_kin_KMBalD_TTT_HCP_A.h"
#include "setup_elastn_HCP.h"
   const int iGdotExpected = 12;
#elif KIN_TYPE == 1
#include "setup_slipGeom.h"
#include "setup_kin_KMBalD_FFF.h"
#include "setup_elastn.h"
   const int iGdotExpected = 1;
#else
#include "setup_slipGeom.h"
#include "setup_kin_VocePL.h"
#include "setup_elastn.h"
   const int iGdotExpected = 1;
#endif


   //////////////////////////////

   double p = 0.0, tK = 300.0;
   std::vector<double> h_state_vec;
   double* h_state;
   {
      std::vector<std::string> names;
      std::vector<bool>        plot;
      std::vector<bool>        state;
      kinetics.getHistInfo(names, h_state_vec, plot, state);
      h_state = &(h_state_vec[0]);
   }

#include "setup_conditions.h"
   double detV = 1.0;
   double eVref = 0.0;
   double e_vecd_n[ecmech::ntvec] = { 0.0 };
   double Cn_quat[ecmech::qdim] = { 1.0, 0.0, 0.0, 0.0 };

   Prob prob(slipGeom, kinetics, elastN,
             dt,
             detV, eVref, p, tK,
             h_state, e_vecd_n, Cn_quat,
             d_vecd_sm, w_veccp_sm);

   Solver solver(prob);

   snls::TrDeltaControl deltaControl;
   deltaControl._deltaInit = 1e0;
   {
      int maxIter = 100;
      solver.setupSolver(maxIter, tolerance, &deltaControl, outputLevel);
   }

   for (int iX = 0; iX < prob.nDimSys; ++iX) {
      solver._x[iX] = 0e0;
   }

   snls::SNLSStatus_t status = solver.solve( );
   if (status != snls::converged) {
      ECMECH_FAIL(__func__, "Solver failed to converge!");
   }
   std::cout << "Function evaluations: " << solver.getNFEvals() << std::endl;
   std::cout << "Last 'rho' in solver: " << solver.getRhoLast() << std::endl;
#ifdef ECMECH_DEBUG
   std::cout << "Slip system shearing rates : ";
   printVec<slipGeom.nslip>(prob.getGdot(), std::cout);
#endif
   EXPECT_TRUE(solver.getNFEvals() == expectedNFEvals) << "Not the expected number of function evaluations";
   {
      const double* gdot = prob.getGdot();
      EXPECT_LT(fabs(gdot[iGdotExpected] - expectedGdotVal), 1e-8) <<
         "Did not get expected value for gdot[iGdotExpected]";
   }
   EXPECT_LT(solver.getRhoLast() - 1.0, 1e-3) << "Final 'rho' from solver not as close to 1 as expected";

   //////////////////////////////////////////////////////////////////////
   //
   // this gives the same as above? does updateH internally, but with
   // beginning-of-step shearing rates which are all set to zero

   typedef ecmech::EosModelConst<false> EosModel;
   EosModel eos;
#include "setup_eos.h"

   double eInt[ecmech::ne] = { 0.0 };
   double stressSvecP[ecmech::nsvp] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                        0.0 };
   static const int iHistLbGdot = evptn::NumHist<SlipGeom, Kinetics, ThermoElastN, EosModel>::iHistLbGdot;
   static const int numHist = evptn::NumHist<SlipGeom, Kinetics, ThermoElastN, EosModel>::numHist;
   double hist[numHist] = { 0.0 };
   std::copy(Cn_quat, Cn_quat + ecmech::qdim, hist + evptn::iHistLbQ);
   std::copy(h_state, h_state + kinetics.nH, hist + evptn::iHistLbH);
   double* gdot = &(hist[iHistLbGdot]); // already zerod
   // do not bother with other stuff (like e_vecd_n) that is all zero above
   //
   double tkelv;
   double sdd[ecmech::nsdd];
   double mtanSD[ecmech::nsvec2];
   //
   evptn::getResponseSngl<SlipGeom, Kinetics, ThermoElastN, EosModel>
      (slipGeom, kinetics, elastN, eos,
      dt,
      tolerance,
      d_svec_kk_sm, w_veccp_sm, volRatio,
      eInt, stressSvecP, hist,
      tkelv, sdd, mtanSD);
   int nFEvals = hist[evptn::iHistA_nFEval];
   std::cout << "Function evaluations: " << nFEvals << std::endl;
#ifdef ECMECH_DEBUG
   std::cout << "Updated hist : ";
   printVec<numHist>(hist, std::cout);

   std::cout << "Slip system shearing rates : ";
   printVec<slipGeom.nslip>(gdot, std::cout);
#endif
   // add 1 to expectedNFEvals because asked for mtanSD
   EXPECT_TRUE(nFEvals == expectedNFEvals + 1) << "Not the expected number of function evaluations";
   EXPECT_LT(fabs(hist[evptn::iHistLbE + 1] - expectedE2), 1e-10) <<
      "Did not get expected value for lattice strain component";
   EXPECT_LT(fabs(hist[evptn::iHistLbQ] - expectedQ1), 1e-8) <<
      "Did not get expected value for quat_1";
   EXPECT_LT(fabs(gdot[iGdotExpected] - expectedGdotVal), 1e-8) <<
      "Did not get expected value for gdot[iGdotExpected]";
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

