#include <gtest/gtest.h>

#include "SNLS_TrDLDenseG.h"

#include "cases/ECMech_cases_fcc_defs.h"

#define STACK_PARAMS

#ifndef DO_FD_CHECK_MTAN
#define DO_FD_CHECK_MTAN 0
#endif

#if DO_FD_CHECK_MTAN
// if doing finite-difference check of tangent stiffness, then change some other defaults
// it is useful to have a mushy rate sensitivity (more rate sensitive) for checking the tangent stiffness as the finite differencing gets a better result
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

static int outputLevel = 0;

#include "test_expectedVals.h"

TEST(ecmech, updst_a)
{
   using namespace ecmech;

#if KIN_TYPE
   using mat_model = matModelEvptn_FCC_B;
#else
   using mat_model = matModelEvptn_FCC_A;
#endif
   mat_model* mmodel = new mat_model();
   matModelBase* mmb = dynamic_cast<matModelBase*>(mmodel);

#include "setup_base.h"
   std::vector<int>           opts; // none
   std::vector<std::string>   strs; // none
   std::vector<double>         params { density0, cvav, tolerance };
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
   DUMPVEC("opts", opts);
   DUMPVEC("params", params);
   DUMPVEC("strs", strs);
   //
   mmb->setExecutionStrategy(ecmech::ExecutionStrategy::CPU);
   mmb->initFromParams(opts, params, strs);
   //
   mmb->complete();

   std::vector<double>       hist_vec;
   {
      std::vector<std::string> names;
      std::vector<bool>        plot;
      std::vector<bool>        state;
      mmb->getHistInfo(names, hist_vec, plot, state);
   }
   double* hist = &(hist_vec[0]);
#if NON_I_QUAT
   double* q_state = &(hist[ecmech::evptn::iHistLbQ]);
   {
      double th = 0.2;
      double et = 0.7;
      q_state[0] = cos(0.5 * th);
      q_state[1] = sin(0.5 * th) * cos(et);
      q_state[2] = sin(0.5 * th) * sin(et);
      q_state[3] = 0.0;
   }
#endif

   // int numHist = hist_vec.size() ; // should equal mmodel->numHist

   constexpr int nPassed = 1; // just do a single point here as a simple example

   mmodel->setOutputLevel(outputLevel); // would not normally do this in a production setting

   static const int iHistLbGdot = mmodel->iHistLbGdot;
   double* gdot = &(hist[iHistLbGdot]);
#if defined(ECMECH_DEBUG) && defined(__ecmech_host_only__)
   std::cout << "Initial hist : ";
   ecmech::printVec(hist, mmodel->numHist, std::cout);
#endif
#include "setup_conditions.h"
   {
      double internal_energy[ecmech::ne] = { 0.0 };
      double cauchy_stress_d6p[ecmech::nsvp] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                           0.0 };
      double tkelv[nPassed];
      double sdd[ecmech::nsdd * nPassed];

      mmb->getResponseECM(dt, def_rate_d6v_sample, spin_vec_sample, rel_vol_ratios,
                          internal_energy, cauchy_stress_d6p, hist, tkelv, sdd, nullptr,
                          nPassed);

      std::cout << "Function evaluations: " << hist[evptn::iHistA_nFEval] << std::endl;
   }
#if defined(ECMECH_DEBUG) && defined(__ecmech_host_only__)
   std::cout << "Updated hist : ";
   ecmech::printVec(hist, mmodel->numHist, std::cout);

   std::cout << "Hardness state : ";
   ecmech::printVec<mat_model::nH>(&(hist[ecmech::evptn::iHistLbH]), std::cout);

   std::cout << "Slip system shearing rates : ";
   ecmech::printVec<mat_model::nslip>(gdot, std::cout);
#endif
   EXPECT_TRUE(hist[evptn::iHistA_nFEval] == expectedNFEvals) << "Not the expected number of function evaluations";
   EXPECT_LT(fabs(hist[evptn::iHistLbE + 1] - expectedE2), 1e-10) <<
      "Did not get expected value for lattice strain component";
   EXPECT_LT(fabs(hist[evptn::iHistLbQ] - expectedQ1), 1e-8) <<
      "Did not get expected value for quat_1";
   EXPECT_LT(fabs(gdot[1] - expectedGdotVal), 1e-8) <<
      "Did not get expected value for gdot[1]";

   delete mmodel;
}

TEST(ecmech, driver_a)
{
   using namespace ecmech;

#if KIN_TYPE
   using mat_model = matModelEvptn_FCC_B;
#else
   using mat_model = matModelEvptn_FCC_A;
#endif
   mat_model* mmodel = new mat_model();
   matModelBase* mmb = dynamic_cast<matModelBase*>(mmodel);

#include "setup_base.h"
   std::vector<int>           opts; // none
   std::vector<std::string>   strs; // none
   std::vector<double>         params { density0, cvav, tolerance };
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
   DUMPVEC("opts", opts);
   DUMPVEC("params", params);
   DUMPVEC("strs", strs);
   //
   mmb->setExecutionStrategy(ecmech::ExecutionStrategy::CPU);
   mmb->initFromParams(opts, params, strs);
   //
   mmb->complete();

   std::vector<double>       hist_vec;
   {
      std::vector<std::string> names;
      std::vector<bool>        plot;
      std::vector<bool>        state;
      mmb->getHistInfo(names, hist_vec, plot, state);
   }
   double* hist = &(hist_vec[0]);
#if NON_I_QUAT
   double* q_state = &(hist[ecmech::evptn::iHistLbQ]);
   {
      double th = 0.2;
      double et = 0.7;
      q_state[0] = cos(0.5 * th);
      q_state[1] = sin(0.5 * th) * cos(et);
      q_state[2] = sin(0.5 * th) * sin(et);
      q_state[3] = 0.0;
   }
#endif

   // int numHist = hist_vec.size() ; // should equal mmodel->numHist

   const int nPassed = 1; // just do a single point here as a simple example

   mmodel->setOutputLevel(outputLevel); // would not normally do this in a production setting

   double relRate = 1e-6;
   double def_rate_d6v_sample[ecmech::nsvp] = { -0.5 * relRate, -0.5 * relRate, 1.0 * relRate,
                                         0.0, 0.0, 0.0,
                                         0.0 };
   // vecsVsa<ecmech::nsvp>(def_rate_d6v_sample, sqr2b3) ; // nope, choose not to do that here
   //
   double def_rate_d5_sample[ecmech::ntvec];
   svecToVecd(def_rate_d5_sample, def_rate_d6v_sample);

   // dt value in setup_conditions.h is meant to stress the implementation --
   // here go with a smaller value to be able to make a nicer curve
   double dt = 0.002 / relRate;
   int nStep = 100;

   double spin_vec_sample[ecmech::nwvec] = { 0.0, 0.0, 0.0 };

   double rel_vol_ratios[ecmech::nvr] = { 1.0, 1.0, 0.0, 0.0 };

   double internal_energy[ecmech::ne] = { 0.0 };
   double cauchy_stress_d6p[ecmech::nsvp] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                        0.0 };
   double tkelv[nPassed];
   double sdd[ecmech::nsdd * nPassed];

#if !(DO_FD_CHECK_MTAN)
   std::cout << "# time, Axial deviatoric stress, h[0], p : " << std::endl;
   double time = 0.0;
#endif
   //
   for (int iStep = 0; iStep<nStep; ++iStep) {

      // update current relative volume from the volumetric deformation rate
      //
      rel_vol_ratios[0] = rel_vol_ratios[1];
      rel_vol_ratios[1] = rel_vol_ratios[0] * exp(def_rate_d6v_sample[ecmech::iSvecP] * dt);
      rel_vol_ratios[3] = rel_vol_ratios[1] - rel_vol_ratios[0];
      rel_vol_ratios[2] = rel_vol_ratios[3] / (dt * 0.5 * (rel_vol_ratios[0] + rel_vol_ratios[1]) );

      mmb->getResponseECM(dt, def_rate_d6v_sample, spin_vec_sample, rel_vol_ratios,
                          internal_energy, cauchy_stress_d6p, hist, tkelv, sdd, nullptr,
                          nPassed);

#if !(DO_FD_CHECK_MTAN)
      time += dt;
      std::cout << time << " "
                << std::setprecision(14) << cauchy_stress_d6p[2] << " "
                << std::setprecision(14) << hist[ecmech::evptn::iHistLbH + 0] << " "
                << std::setprecision(14) << cauchy_stress_d6p[iSvecP] << " "
                << std::endl;

      // std::cout << "hist : " ;
      // printVec(hist, numHist, std::cout) ;
#endif
   }

#if KIN_TYPE && !(DO_FD_CHECK_MTAN)
   EXPECT_LT(fabs(cauchy_stress_d6p[2] - 0.006664208062085), 1e-10) <<
      "Did not get expected value for stress component";
   EXPECT_LT(fabs(hist[ecmech::evptn::iHistLbH + 0] - 87.96284116155), 1e-8) <<
      "Did not get expected value for history variable";
   EXPECT_LT(fabs(cauchy_stress_d6p[iSvecP] - 0.00332519207297), 1e-10) <<
      "Did not get expected value for stress component";
#endif

#if DO_FD_CHECK_MTAN
   {
      //
      // do another step, and do finite differencing to check mtanSD

      std::vector<double> hist_ref(hist, hist + mmodel->numHist);
      std::vector<double> internal_energy_ref(internal_energy, internal_energy + ecmech::ne);
      std::vector<double> cauchy_stress_d6p_ref(cauchy_stress_d6p, cauchy_stress_d6p + ecmech::nsvp);
      double v_ref = rel_vol_ratios[1];

      rel_vol_ratios[0] = v_ref;
      rel_vol_ratios[1] = rel_vol_ratios[0] * exp(def_rate_d6v_sample[ecmech::iSvecP] * dt);
      rel_vol_ratios[3] = rel_vol_ratios[1] - rel_vol_ratios[0];
      rel_vol_ratios[2] = rel_vol_ratios[3] / (dt * 0.5 * (rel_vol_ratios[0] + rel_vol_ratios[1]) );

      double mtanSD_an[ecmech::nsvec2];
      mmb->getResponseECM(dt, def_rate_d6v_sample, spin_vec_sample, rel_vol_ratios,
                          internal_energy, cauchy_stress_d6p, hist, tkelv, sdd, mtanSD_an,
                          nPassed);


      double cauchy_stress[ecmech::nsvec];
      svecpToSvec(cauchy_stress, cauchy_stress_d6p);
#if defined(ECMECH_DEBUG) && defined(__ecmech_host_only__)
      std::cout << "mtanSD_an : " << std::endl;
      printMat<ecmech::nsvec>(mtanSD_an, std::cout);
#endif
      double def_rate_d6v_sample_pert[ecmech::nsvp];
      const double pertVal = 1e-8 * relRate;
      double mtanSD_fd[ecmech::nsvec2];
      //
      double internal_energy_pert[ecmech::ne];
      double cauchy_stress_d6p_pert[ecmech::nsvp];
      //
      for (int jSvec = 0; jSvec<ecmech::nsvec; ++jSvec) {
         std::copy(def_rate_d6v_sample, def_rate_d6v_sample + ecmech::nsvp, def_rate_d6v_sample_pert);
         if (jSvec < 3) {
            def_rate_d6v_sample_pert[jSvec] += pertVal;
            double d_kk = def_rate_d6v_sample_pert[0] + def_rate_d6v_sample_pert[1] + def_rate_d6v_sample_pert[2];
            def_rate_d6v_sample_pert[ecmech::iSvecP] += d_kk;
            def_rate_d6v_sample_pert[0] += (-ecmech::onethird * d_kk);
            def_rate_d6v_sample_pert[1] += (-ecmech::onethird * d_kk);
            def_rate_d6v_sample_pert[2] += (-ecmech::onethird * d_kk);
         }
         else {
            // factor of 2 to go with l_ddsdde_gamma being true in call to mtan_conv_sd_svec ;
            def_rate_d6v_sample_pert[jSvec] += 0.5 * pertVal;
         }
         //
         rel_vol_ratios[0] = v_ref;
         rel_vol_ratios[1] = rel_vol_ratios[0] * exp(def_rate_d6v_sample_pert[ecmech::iSvecP] * dt);
         rel_vol_ratios[3] = rel_vol_ratios[1] - rel_vol_ratios[0];
         rel_vol_ratios[2] = rel_vol_ratios[3] / (dt * 0.5 * (rel_vol_ratios[0] + rel_vol_ratios[1]) );

         std::copy(internal_energy_ref.begin(), internal_energy_ref.end(), internal_energy_pert);
         std::copy(cauchy_stress_d6p_ref.begin(), cauchy_stress_d6p_ref.end(), cauchy_stress_d6p_pert);

         double tkelv_pert[nPassed];
         double sdd_pert[ecmech::nsdd * nPassed];

         std::copy(hist_ref.begin(), hist_ref.end(), hist); // make hist equal to hist_ref again

         mmb->getResponseECM(dt, def_rate_d6v_sample_pert, spin_vec_sample, rel_vol_ratios,
                             internal_energy_pert, cauchy_stress_d6p_pert, hist, tkelv_pert, sdd_pert, nullptr,
                             nPassed);

         double cauchy_stress_pert[ecmech::nsvec];
         svecpToSvec(cauchy_stress_pert, cauchy_stress_d6p_pert);
         //
         for (int iSvec = 0; iSvec<ecmech::nsvec; ++iSvec) {
            // divide by dt because tangent gets converted to a per-strain-increment type quantity
            mtanSD_fd[ECMECH_NN_INDX(iSvec, jSvec, ecmech::nsvec)] = (cauchy_stress_pert[iSvec] - cauchy_stress[iSvec]) / pertVal / dt;
         }
      }

#if defined(ECMECH_DEBUG) && defined(__ecmech_host_only__)
      std::cout << "mtanSD_fd : " << std::endl;
      printMat<ecmech::nsvec>(mtanSD_fd, std::cout);
#endif
      // do not bother restoring things to evaluation at non-perturbed condition

      for (int iiMtan = 0; iiMtan<ecmech::nsvec2; ++iiMtan) {
         EXPECT_LT(fabs(mtanSD_fd[iiMtan] - mtanSD_an[iiMtan]), 1e-3) <<
            "Analytic and finite-differenced mtan differ by more than expected";
      }
   }
#endif // DO_FD_CHECK_MTAN

   delete mmodel;
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

