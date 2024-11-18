#include <gtest/gtest.h>

#include "SNLS_TrDLDenseG.h"
#include "ECMech_util.h"

#include "cases/ECMech_cases_fcc_defs.h"
#include "cases/ECMech_cases_bcc_defs.h"
#include "cases/ECMech_cases_hcp_defs.h"

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
   using SlipGeom = SlipGeom_BCC_A ;
   using Kinetics = Kin_KMBalD_TFF;
   using ThermoElastN =  EVPTN_cubic;
#elif KIN_TYPE == 2
   using SlipGeom = SlipGeom_HCP_A;
   using Kinetics = Kin_HCP_A;
   using ThermoElastN =  EVPTN_hex;
#elif KIN_TYPE == 1
   using SlipGeom = SlipGeomFCC;
   using Kinetics = Kin_KMBalD_FFF;
   using ThermoElastN =  EVPTN_cubic;
#else
   using SlipGeom = SlipGeomFCC;
   using Kinetics = Kin_Voce;
   using ThermoElastN =  EVPTN_cubic;
#endif
   using ProblemState = evptn::ProblemState<SlipGeom, Kinetics, ThermoElastN, EosModelConst<false>>;

   using Prob = evptn::EvptnUpdstProblem<SlipGeom, Kinetics, ThermoElastN, ProblemState>;
   using Solver = snls::SNLSTrDlDenseG<Prob>;

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

   double tkelv = 300.0;
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

   constexpr int numHist1 = evptn::NumHist<SlipGeom, Kinetics, ThermoElastN, EosModelConst<false>>::numHist;
   double hist2[numHist1] = {};

   ProblemState prob_state(hist2, nullptr, tkelv, def_rate_d6v_sample, spin_vec_sample, rel_vol_ratios, dt);

   prob_state.quat_n[0] = 1.0;
   for (int iqdim = 1; iqdim < ecmech::qdim; iqdim++) {
      prob_state.quat_n[iqdim] = 0.0;
   }

   for (size_t iH = 0; iH < h_state_vec.size(); iH++) {
      prob_state.h_state_u[iH] = h_state[iH];
   }

   prob_state.energy_new = 0.0;
   prob_state.pressure_EOS = 0.0;

   Prob prob(slipGeom, kinetics, elastN, prob_state); 

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
   {
      double gdot[slipGeom.nslip] = {};
      double junk = 0.0;
      double junk_vec[ecmech::qdim] = {};
      double elast_strain_d5[ecmech::ntvec] = {};
      prob.stateFromX(elast_strain_d5, junk_vec, solver._x);
      prob.get_slip_contribution(junk, junk, gdot, elast_strain_d5);
      printVec<slipGeom.nslip>(gdot, std::cout);
   }
#endif
   EXPECT_TRUE(solver.getNFEvals() == expectedNFEvals) << "Not the expected number of function evaluations";
   {
      double gdot[slipGeom.nslip] = {};
      double junk = 0.0;
      double junk_vec[ecmech::qdim] = {};
      double elast_strain_d5[ecmech::ntvec] = {};
      prob.stateFromX(elast_strain_d5, junk_vec, solver._x);
      prob.get_slip_contribution(junk, junk, gdot, elast_strain_d5);
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

   double internal_energy[ecmech::ne] = { 0.0 };
   double cauchy_stress_d6p[ecmech::nsvp] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                        0.0 };
   static const int iHistLbGdot = evptn::NumHist<SlipGeom, Kinetics, ThermoElastN, EosModel>::iHistLbGdot;
   static const int numHist = evptn::NumHist<SlipGeom, Kinetics, ThermoElastN, EosModel>::numHist;
   double hist[numHist] = { 0.0 };
   std::copy(prob_state.quat_n, prob_state.quat_n + ecmech::qdim, hist + evptn::iHistLbQ);
   std::copy(h_state, h_state + kinetics.nH, hist + evptn::iHistLbH);
   double* gdot = &(hist[iHistLbGdot]); // already zerod
   // do not bother with other stuff (like e_vecd_n) that is all zero above
   //
   double tkelv2;
   double sdd[ecmech::nsdd];
   double mtanSD[ecmech::nsvec2];
   //
   evptn::getResponseSngl<SlipGeom, Kinetics, ThermoElastN, EosModel>
      (slipGeom, kinetics, elastN, eos,
      dt,
      tolerance,
      def_rate_d6v_sample, spin_vec_sample, rel_vol_ratios,
      internal_energy, cauchy_stress_d6p, hist,
      tkelv2, sdd, mtanSD);
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

