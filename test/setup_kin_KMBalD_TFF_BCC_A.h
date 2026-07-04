/**
 * @file setup_kin_KMBalD_TFF_BCC_A.h
 *
 * @brief Fragment (see `setup_base.h` for the general `setup_*.h` inclusion pattern)
 * that sets parameters for `Kin_KMBalD_TFF` (`KineticsKMBalD<true, false, false, false,
 * 1>`, see `ECMech_cases.h` / `kinetics/ECMech_kinetics_KMBalD.h`), used by
 * `matModelEvptn_BCC_B` (`"evptn_BCC_B"`). The "TFF" name is the three leading template
 * bools: `withGAthermal = true` (CRSS ĝ and the MTS-normalizing stress `τ_a` are split
 * -- see `setup_kin_KMBalD_TTT_HCP_A.h`'s doc for how that changes the parameter list),
 * `pOne = false`, `qOne = false` (general `p`/`q` exponents, same values as the FFF
 * variant). `perSS = false` as well, so `c_1`/`go`/`s` remain single crystal-wide
 * values even though `SlipGeom_BCC_A` has 12 slip systems.
 *
 * Same parameter meanings and order as `setup_kin_KMBalD_FFF.h` (this file is
 * numerically identical to it, just grouped slightly differently in the initializer
 * list); see that file's doc for the full per-parameter breakdown.
 */
{
   double
      shear_modulus = 1.0,
      tkelv_ref = 300.,
      c_1 = 20000.,
      tau_a = 0.004,
      p = 0.28,
      q = 1.34,
      gam_wo = 20.,
      gam_ro = 1e3,
      wrD = 0.02,
      go = 10e-5,
      s = 5e-5;
   double
      k1 = 100.0,
      k2o = 10.0,
      ninv = 0.05,
      gamma_o = 1e-6;
   double
      rho_dd_init = 0.25;
   std::vector<double> paramsThese {
      shear_modulus, tkelv_ref,
      c_1,
      tau_a, p, q, gam_wo, gam_ro, wrD,
      go,
      s,
      k1, k2o, ninv, gamma_o,
      rho_dd_init
   };
#ifdef STACK_PARAMS
   params.insert(params.end(), paramsThese.begin(), paramsThese.end());
#else
   kinetics.setParams(paramsThese);
#endif
}
