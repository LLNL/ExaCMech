/**
 * @file setup_kin_KMBalD_TTT_HCP_A.h
 *
 * @brief Fragment (see `setup_base.h` for the general `setup_*.h` inclusion pattern)
 * that sets parameters for `Kin_HCP_A` (`KineticsKMBalD<true, true, true, true,
 * SlipGeom_HCP_A::nslip>`, see `ECMech_cases_hcp_defs.h` /
 * `kinetics/ECMech_kinetics_KMBalD.h`), used by `matModelEvptn_HCP_A`
 * (`"evptn_HCP_A"`). The "TTT" name is the three leading template bools, all true:
 * `withGAthermal = true` (CRSS ĝ and the MTS-normalizing stress `τ_a` are split, hence
 * both `go` (initial ĝ, athermal-floor-like) and `tau_a` appear below), `pOne = true`
 * and `qOne = true` (the MTS activation-energy exponents are fixed at 1 -- hence `p =
 * 1.0`, `q = 1.0` here, unlike the FFF/TFF variants' `p = 0.28`, `q = 1.34`).
 *
 * Unlike `setup_kin_KMBalD_FFF.h`/`_TFF_BCC_A.h`, this model also sets `perSS = true`
 * with `nVPer = SlipGeom_HCP_A::nslip = 24` -- so `c_1`, `go`, and `s` are each given as
 * 24 individual per-slip-system values below (grouped 3 + 3 + 6 + 12, matching HCP_A's
 * basal/prismatic/pyramidal-⟨a⟩/pyramidal-⟨c+a⟩ slip-system families) rather than one
 * shared crystal-wide value each. All 24 entries happen to be identical here (same
 * physical value repeated per system), but the *shape* of the parameter list itself is
 * what exercises the `perSS = true` code path.
 *
 * See `setup_kin_KMBalD_FFF.h`'s doc for what each parameter name physically means.
 */
{
   double
      shear_modulus = 1.0,
      tkelv_ref = 300.,
      c_1 = 20000.,
      tau_a = 0.004,
      p = 1.0,
      q = 1.0,
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
      c_1, c_1, c_1,
      c_1, c_1, c_1,
      c_1, c_1, c_1, c_1, c_1, c_1,
      c_1, c_1, c_1, c_1, c_1, c_1, c_1, c_1, c_1, c_1, c_1, c_1,
      tau_a, p, q, gam_wo, gam_ro, wrD,
      go, go, go,
      go, go, go,
      go, go, go, go, go, go,
      go, go, go, go, go, go, go, go, go, go, go, go,
      s, s, s,
      s, s, s,
      s, s, s, s, s, s,
      s, s, s, s, s, s, s, s, s, s, s, s,
      k1, k2o, ninv, gamma_o,
      rho_dd_init
   };
#ifdef STACK_PARAMS
   params.insert(params.end(), paramsThese.begin(), paramsThese.end());
#else
   kinetics.setParams(paramsThese);
#endif
}
