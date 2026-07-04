/**
 * @file setup_kin_VocePL.h
 *
 * @brief Fragment (see `setup_base.h` for the general `setup_*.h` inclusion pattern)
 * that sets parameters for linear-Voce power-law kinetics. Used two ways depending on
 * which test includes it: for `Kin_Voce` (`kinetics/ECMech_kinetics_VocePL.h`, used by
 * `matModelEvptn_FCC_A`/`"evptn_FCC_A"`) in `test_evptn.cxx`/`test_updst.cxx`/
 * `test_px.cxx`; and, when `ANISO_HARDENING` is defined, for `KineticsAnisoVocePL`
 * (declared in `test_aniso_kinetics_VocePL.h`, a per-slip-system-hardening test variant)
 * in `test_aniso_hardening.cxx`. See `kinetics/ECMech_kinetics_VocePL.h`'s `setParams`
 * doc for what each shared name means; briefly: `shear_modulus`/`xm`/`gam_w`
 * (power-law slip-rate parameters), `h0`/`tausi`/`taus0`/`xms`/`gamss0` (Voce
 * hardening-law parameters -- initial hardening rate, initial CRSS, saturation-stress
 * reference, saturation-stress rate-sensitivity exponent, and saturation-stress
 * reference rate), and a trailing `hdn_init` (initial hardening state, here just set
 * equal to `tausi`).
 *
 * Three build-time macros toggle which variant gets built:
 * - `XM_MUSHY`: increases the rate-sensitivity exponent `xm` from `0.01` to `0.1` (a
 *   "mushier"/more rate-sensitive material), used by `test_updst.cxx`'s finite-difference
 *   tangent-stiffness check since a more rate-sensitive response finite-differences more
 *   accurately.
 * - `KIN_NONLINEAR`: switches to the nonlinear-Voce parameter list, which adds an extra
 *   `xmprime` parameter between `taus0` and `xms`.
 * - `ANISO_HARDENING`: repeats `tausi` 12 times (one per FCC slip system) instead of
 *   once, for `KineticsAnisoVocePL`'s per-slip-system initial CRSS. `hdn_init` stays a
 *   single trailing value regardless -- `test_aniso_hardening.cxx` broadcasts it across
 *   all 12 slip systems itself (via `std::fill`) rather than needing 12 separate values
 *   here.
 */
#if XM_MUSHY
#define XM_VAL 0.1
#else
#define XM_VAL 0.01
#endif

{
   double shear_modulus = 1.0, xm = XM_VAL, gam_w = 1.0;
   double h0 = 200e-5, tausi = 100e-5, taus0 = 400e-5, xms = 0.05, gamss0 = 1.0e-6;
#ifdef KIN_NONLINEAR
   double xmprime = 2.0;
   #ifdef ANISO_HARDENING
   std::vector<double> paramsThese {
      shear_modulus, xm, gam_w,
      h0,
      tausi, tausi, tausi, tausi,
      tausi, tausi, tausi, tausi,
      tausi, tausi, tausi, tausi,
      taus0, xmprime, xms, gamss0,
      tausi // hdn_init
   };
   #else
   std::vector<double> paramsThese {
      shear_modulus, xm, gam_w,
      h0, tausi, taus0, xmprime, xms, gamss0,
      tausi // hdn_init
   };
   #endif
#else
   #ifdef ANISO_HARDENING
   std::vector<double> paramsThese {
      shear_modulus, xm, gam_w,
      h0,
      tausi, tausi, tausi, tausi,
      tausi, tausi, tausi, tausi,
      tausi, tausi, tausi, tausi,
      taus0, xms, gamss0,
      tausi // hdn_init
   };
   #else
   std::vector<double> paramsThese {
      shear_modulus, xm, gam_w,
      h0, tausi, taus0, xms, gamss0,
      tausi // hdn_init
   };
   #endif
#endif
#ifdef STACK_PARAMS
   params.insert(params.end(), paramsThese.begin(), paramsThese.end());
#else
   kinetics.setParams(paramsThese);
#endif
}
