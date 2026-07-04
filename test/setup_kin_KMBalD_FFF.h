/**
 * @file setup_kin_KMBalD_FFF.h
 *
 * @brief Fragment (see `setup_base.h` for the general `setup_*.h` inclusion pattern)
 * that sets parameters for `Kin_KMBalD_FFF` (`KineticsKMBalD<false, false, false,
 * false, 1>`, see `ECMech_cases.h` / `kinetics/ECMech_kinetics_KMBalD.h`), used by
 * `matModelEvptn_FCC_B` (`"evptn_FCC_B"`). The "FFF" name is the three leading template
 * bools in order: `withGAthermal = false` (no athermal-floor/MTS-normalizing-stress
 * split -- CRSS ĝ is used directly as both), `pOne = false` and `qOne = false` (the MTS
 * activation-energy exponents `p`/`q` are general, not fixed at 1), matching `p = 0.28`,
 * `q = 1.34` below. `perSS = false` here too, so `c_1`/`go`/`s` are single crystal-wide
 * values rather than per-slip-system ones (contrast `setup_kin_KMBalD_TTT_HCP_A.h`,
 * which repeats these per slip system).
 *
 * See `KineticsKMBalD::setParams`'s doc for the full parameter list/order; briefly:
 * `shear_modulus`/`tkelv_ref` (MTS reference shear modulus and temperature), `c_1`
 * (thermal-energy prefactor scale), `tau_a` (athermal/MTS reference stress), `p`/`q`
 * (MTS activation-energy exponents), `gam_wo`/`gam_ro` (forward/reverse reference shear
 * rates), `wrD` (phonon-drag weighting), `go`/`s` (initial CRSS and Voce-like
 * saturation-stress-law scale), `k1`/`k2o`/`ninv` (Kocks-Mecking hardening-law
 * coefficients), `gamma_o` (reference shear rate for the hardening law), and
 * `rho_dd_init` (initial relative dislocation density -- the single hardening state
 * variable this model tracks).
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
      shear_modulus, tkelv_ref, c_1, tau_a, p, q, gam_wo, gam_ro, wrD, go, s,
      k1, k2o, ninv, gamma_o,
      rho_dd_init
   };
#ifdef STACK_PARAMS
   params.insert(params.end(), paramsThese.begin(), paramsThese.end());
#else
   kinetics.setParams(paramsThese);
#endif
}
