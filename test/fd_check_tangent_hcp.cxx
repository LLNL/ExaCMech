// Finite-difference verification of the consistent material tangent returned
// by evptn::getResponseSngl for an HCP (hexagonal) material at a general
// (non-identity) initial orientation, with both deviatoric and volumetric
// loading directions probed.
//
// The returned mtanSD is d(sigma_svec)/d(strain increment) with engineering
// (gamma) shear convention. Each FD column re-runs the full nonlinear solve
// from the same beginning-of-step state with a perturbed strain increment.

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <vector>
#include <string>

#include "ECMech_cases.h"
#include "cases/ECMech_cases_hcp_defs.h"

int main(int argc, char* argv[])
{
   using namespace ecmech;
   // loading scale: 1.0 reproduces the standard test conditions (deep plastic
   // flow); small values (e.g. 1e-6) probe the near-elastic limit where the
   // only legitimate tangent error is FD noise
   double load_scale = (argc > 1) ? atof(argv[1]) : 1.0;

   using SlipGeom = SlipGeom_HCP_A;
   using Kinetics = Kin_HCP_A;
   using ThermoElastN = EVPTN_hex;
   using EosModel = EosModelConst<false>;

   double density0 = 3.0, cvav = 2.0e-5;
   double tolerance = 1e-12;

   SlipGeom slipGeom;
   Kinetics kinetics(slipGeom.nslip);
   ThermoElastN elastN;
#include "setup_slipGeom_HCP.h"
#include "setup_kin_KMBalD_TTT_HCP_A.h"
#include "setup_elastn_HCP.h"
   EosModel eos;
#include "setup_eos.h"

   double tkelv_in = 300.0;
   std::vector<double> h_state_vec;
   {
      std::vector<std::string> names;
      std::vector<bool> plot, state;
      kinetics.getHistInfo(names, h_state_vec, plot, state);
   }

   const double dt = 1e-1;

   // base loading: unit-magnitude deviatoric stretch plus spin (as in test_evptn)
   double def_rate_base[ecmech::nsvp] = { -0.5, -0.5, 1.0, 0.0, 0.0, 0.0, 0.0 };
   vecsVsa<ecmech::nsvp>(def_rate_base, ecmech::sqr2b3 * load_scale);
   double spin_vec_sample[ecmech::nwvec] = { 0.0, 0.0, 0.5 * load_scale };

   // general orientation -- essential: frame errors in the volumetric coupling
   // column vanish when the c-axis aligns with a sample axis
   double q0[ecmech::qdim] = { 0.9, 0.3, -0.25, 0.15 };
   {
      double qn = 0.0;
      for (int i = 0; i < ecmech::qdim; ++i) { qn += q0[i] * q0[i]; }
      qn = 1.0 / sqrt(qn);
      for (int i = 0; i < ecmech::qdim; ++i) { q0[i] *= qn; }
   }

   constexpr int numHist = evptn::NumHist<SlipGeom, Kinetics, ThermoElastN, EosModel>::numHist;

   auto runModel = [&](const double* def_rate6, double vol_incr,
                       double* stress_svec_out, double* mtanSD) {
      double hist[numHist] = { 0.0 };
      std::copy(q0, q0 + ecmech::qdim, hist + evptn::iHistLbQ);
      std::copy(h_state_vec.begin(), h_state_vec.end(), hist + evptn::iHistLbH);
      double internal_energy[ecmech::ne] = { 0.0 };
      double cauchy_stress_d6p[ecmech::nsvp] = { 0.0 };
      double vNew = exp(vol_incr);
      double rel_vol_ratios[ecmech::nvr] = { 1.0, vNew, (vNew - 1.0) / dt, vNew - 1.0 };
      double tkelv2;
      double sdd[ecmech::nsdd];
      evptn::getResponseSngl<SlipGeom, Kinetics, ThermoElastN, EosModel>
         (slipGeom, kinetics, elastN, eos,
         dt, tolerance,
         def_rate6, spin_vec_sample, rel_vol_ratios,
         internal_energy, cauchy_stress_d6p, hist,
         tkelv2, sdd, mtanSD);
      // full Cauchy stress in svec: deviatoric part minus pressure on normals
      for (int i = 0; i < ecmech::nsvec; ++i) {
         stress_svec_out[i] = cauchy_stress_d6p[i];
      }
      for (int i = 0; i < 3; ++i) {
         stress_svec_out[i] -= cauchy_stress_d6p[6];
      }
   };

   double mtan[ecmech::nsvec2];
   double sig0[ecmech::nsvec];
   const double vol0 = (argc > 3) ? atof(argv[3]) : 0.0;
   runModel(def_rate_base, vol0, sig0, mtan);

   // finite differences: columns are svec strain-increment directions,
   // engineering shear convention for columns 3-5
   const double h = (argc > 2) ? atof(argv[2]) : 1e-6;
   double mtan_fd[ecmech::nsvec2];
   for (int j = 0; j < ecmech::nsvec; ++j) {
      double def_rate[ecmech::nsvp];
      std::copy(def_rate_base, def_rate_base + ecmech::nsvp, def_rate);
      double vol_incr = vol0;
      if (j < 3) {
         // normal strain increment h e_j x e_j : deviatoric part + trace part
         for (int k = 0; k < 3; ++k) {
            def_rate[k] += (((k == j) ? 1.0 : 0.0) - onethird) * h / dt;
         }
         vol_incr = vol0 + h;
      }
      else {
         // engineering shear increment gamma = h -> tensor component h/2
         def_rate[j] += 0.5 * h / dt;
      }
      double sigp[ecmech::nsvec];
      runModel(def_rate, vol_incr, sigp, nullptr);
      for (int i = 0; i < ecmech::nsvec; ++i) {
         mtan_fd[ECMECH_NN_INDX(i, j, ecmech::nsvec)] = (sigp[i] - sig0[i]) / h;
      }
   }

   double amax = 0.0, emax = 0.0, emax_vol = 0.0;
   for (int ij = 0; ij < ecmech::nsvec2; ++ij) {
      amax = fmax(amax, fabs(mtan_fd[ij]));
   }
   std::cout << std::setprecision(5);
   std::cout << "\nmtanSD (analytic) | mtan_fd (finite difference), row by row:\n";
   for (int i = 0; i < ecmech::nsvec; ++i) {
      for (int j = 0; j < ecmech::nsvec; ++j) {
         std::cout << std::setw(11) << mtan[ECMECH_NN_INDX(i, j, ecmech::nsvec)];
      }
      std::cout << "   |";
      for (int j = 0; j < ecmech::nsvec; ++j) {
         std::cout << std::setw(11) << mtan_fd[ECMECH_NN_INDX(i, j, ecmech::nsvec)];
      }
      std::cout << "\n";
      for (int j = 0; j < ecmech::nsvec; ++j) {
         int ij = ECMECH_NN_INDX(i, j, ecmech::nsvec);
         double e = fabs(mtan[ij] - mtan_fd[ij]) / amax;
         emax = fmax(emax, e);
      }
   }
   // volumetric response lives in the sum over the three normal-strain columns
   // (columns j<3 each carry trace h); isolate max error there for reporting
   for (int i = 0; i < ecmech::nsvec; ++i) {
      double a = 0.0, f = 0.0;
      for (int j = 0; j < 3; ++j) {
         a += mtan[ECMECH_NN_INDX(i, j, ecmech::nsvec)];
         f += mtan_fd[ECMECH_NN_INDX(i, j, ecmech::nsvec)];
      }
      emax_vol = fmax(emax_vol, fabs(a - f) / amax);
   }

   std::cout << "\nmax |mtan - mtan_fd| / max|mtan_fd|            : " << emax << "\n";
   std::cout << "same, row-sums over normal columns (vol probe) : " << emax_vol << "\n";
   return 0;
}
