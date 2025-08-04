#include <gtest/gtest.h>

#include "SNLS_TrDLDenseG.h"

#include "cases/ECMech_cases_fcc_defs.h"

#include <random>

#define STACK_PARAMS

#ifndef KIN_TYPE
#define KIN_TYPE 1
#endif

static int outputLevel = 0;

TEST(ecmech, px_a)
{
   // can adjust these to change the computational workload
   //
   constexpr int nPassed = 16;
   constexpr int nStep = 100;
   //
   const double weight = 1.0 / (double) (nPassed);

   using namespace ecmech;

   // it would be nice to try out calls like
   // matModelBase* mmb = makeMatModel("evptn_FCC_A");
   // here, but that does not play nicely with the parameter munging machinery used here
#if KIN_TYPE
   matModelEvptn_FCC_B* mmodel = new matModelEvptn_FCC_B();
#else
   matModelEvptn_FCC_A* mmodel = new matModelEvptn_FCC_A();
#endif
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

   mmb->setOutputLevel(outputLevel); // would not normally do this in a production setting

   std::vector<double> histInit_vec;
   {
      std::vector<std::string> names;
      std::vector<bool>        plot;
      std::vector<bool>        state;
      mmb->getHistInfo(names, histInit_vec, plot, state);
   }

   //
   // set up hist and other state information
   //
   const int numHist = mmb->getNumHist();
   std::vector<double> V_hist(numHist * nPassed, 0.0);
   {
      // Turns out with this is implementation dependent and on newer macos arm systems it doesn't follow our old rand distribution implementations :(
      // std::default_random_engine gen(1);
      std::minstd_rand0 gen(1);
      std::normal_distribution<double> distrib(0.0, 1.0);
      for (int iPassed = 0; iPassed < nPassed; iPassed++) {
         double* hist = &(V_hist[numHist * iPassed]);
         for (int iHist = 0; iHist < numHist; iHist++) {
            hist[iHist] = histInit_vec[iHist];
         }

         double* q_state = &(hist[ecmech::evptn::iHistLbQ]);
         q_state[0] = distrib(gen);
         q_state[1] = distrib(gen);
         q_state[2] = distrib(gen);
         q_state[3] = distrib(gen);
         vecsVNormalize<ecmech::qdim>(q_state);
      }
   }

   double relRate = 1e-6;
   double dt = 0.002 / relRate;

   double def_rate_d6v_sample[ecmech::nsvp] = { -0.5 * relRate, -0.5 * relRate, 1.0 * relRate,
                                         0.0, 0.0, 0.0,
                                         0.0 };
   double V_def_rate_d6v_sample[ecmech::nsvp * nPassed];
   for (int iPassed = 0; iPassed < nPassed; iPassed++) {
      int pOffsetSVP = ecmech::nsvp * iPassed;
      std::copy(def_rate_d6v_sample, def_rate_d6v_sample + ecmech::nsvp, &(V_def_rate_d6v_sample[pOffsetSVP]));
   }

   std::vector<double> spin_vec_sample_vec(ecmech::nwvec * nPassed, 0.0);
   double* V_spin_vec_sample = &(spin_vec_sample_vec[0]);

   std::vector<double> rel_vol_ratios_vec(ecmech::nvr * nPassed, 1.0); // not really 1 for all entries, but this works given what happens below
   double* V_rel_vol_ratios = &(rel_vol_ratios_vec[0]);

   std::vector<double> internal_energy_vec(ecmech::ne * nPassed, 0.0);
   double* V_internal_energy = &(internal_energy_vec[0]);

   std::vector<double> cauchy_stress_d6p_vec(ecmech::nsvp * nPassed, 0.0);
   double* V_cauchy_stress_d6p = &(cauchy_stress_d6p_vec[0]);

   double V_tkelv[nPassed];
   double V_sdd[ecmech::nsdd * nPassed];
   //
   std::cout << "# time, Axial deviatoric stress : " << std::endl;
   double time = 0.0;
   //
   double sAvg = 0.0;
   for (int iStep = 0; iStep<nStep; ++iStep) {
      //
      time += dt;

      // update current relative volume from the volumetric deformation rate
      //
      for (int iPassed = 0; iPassed < nPassed; iPassed++) {
         int pOffsetSVP = ecmech::nsvp * iPassed;
         int pOffsetVR = ecmech::nvr * iPassed;
         V_rel_vol_ratios[0 + pOffsetVR] = V_rel_vol_ratios[1 + pOffsetVR];
         V_rel_vol_ratios[1 + pOffsetVR] = V_rel_vol_ratios[0 + pOffsetVR] * exp(V_def_rate_d6v_sample[ecmech::iSvecP + pOffsetSVP] * dt);
         V_rel_vol_ratios[3 + pOffsetVR] = V_rel_vol_ratios[1 + pOffsetVR] - V_rel_vol_ratios[0 + pOffsetVR];
         V_rel_vol_ratios[2 + pOffsetVR] = V_rel_vol_ratios[3 + pOffsetVR] /
                                     (dt * 0.5 * (V_rel_vol_ratios[0 + pOffsetVR] + V_rel_vol_ratios[1 + pOffsetVR]) );
      }

      mmb->getResponseECM(dt,
                          V_def_rate_d6v_sample, V_spin_vec_sample, V_rel_vol_ratios,
                          V_internal_energy, V_cauchy_stress_d6p, V_hist.data(), V_tkelv, V_sdd, nullptr,
                          nPassed);

      sAvg = 0.0;
      for (int iPassed = 0; iPassed<nPassed; ++iPassed) {
         sAvg += weight * V_cauchy_stress_d6p[iPassed * ecmech::nsvp + 2];
      }

      std::cout << time << " "
                << std::setprecision(14) << sAvg
                << std::endl;
   }

#if KIN_TYPE
   EXPECT_LT(fabs(sAvg - 0.00816331652264), 1e-10) << "Did not get expected value";
#else
   EXPECT_LT(fabs(sAvg - 0.00344825801180), 1e-10) << "Did not get expected value";
#endif

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

