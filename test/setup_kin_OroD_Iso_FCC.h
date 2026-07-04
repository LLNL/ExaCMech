/**
 * @file setup_kin_OroD_Iso_FCC.h
 *
 * @brief Fragment (see `setup_base.h` for the general `setup_*.h` inclusion pattern)
 * that sets parameters for `Kin_OroD_Iso_FCC` (`KineticsOrowanD<false, false, false,
 * true, false, 1, SlipGeomFCC, true>`, see `ECMech_cases_fcc_defs.h` /
 * `kinetics/ECMech_kinetics_OrowanD.h`), used by `matModelEvptn_FCC_C`
 * (`"evptn_FCC_C"`). The 4th template argument, `isotropic = true`, is why `inter_mat`
 * below is a single scalar rather than the full `nslip × nslip` forest-interaction
 * matrix -- contrast `setup_kin_OroD_Iso_FCC_ns.h`, which supplies the full matrix (with
 * every entry equal, so numerically equivalent to this file, but exercising the
 * `isotropic = false` code path instead).
 *
 * See `KineticsOrowanD::setParams`'s doc for the full parameter list; briefly:
 * `shear_modulus_ref`/`tkelv_ref` (MTS reference shear modulus/temperature),
 * `berg_mag`/`lbar` (Burgers-vector magnitude and mean free path), `fD` (phonon-drag
 * reference frequency), `c_1` (thermal-energy prefactor scale), `tau_a`/`p`/`q` (MTS
 * athermal reference stress and activation-energy exponents), `c2` (drag-stress
 * scale), `gam_ro`/`wrD` (reverse reference shear rate and phonon-drag weighting),
 * `inter_mat` (forest-interaction coefficient), `c_ann`/`d_ann`/`c_trap`/`c_mult`
 * (dislocation annihilation/trapping/multiplication coefficients), and `qM`/`qT`
 * (initial mobile/total dislocation densities, one pair per slip system).
 */
{
   double
      shear_modulus_ref = 1.0, //MBar
      tkelv_ref = 300., //K
      berg_mag = 1.0e-4, // microns
      lbar = 10.0 * berg_mag, // microns
#ifdef KIN_BCC
      fD = 1.0e5, //micro-sec
#else
      fD = 1.0e5, //micro-sec
#endif
      // Getting the units right here can be a pain...
      // c_1 = 0.65e-4 * berg_mag * berg_mag * berg_mag / 1.3806504e-23, // g_0 * b^3 / kB
      // Taken from the other example
      c_1 = 20000.,
      tau_a = 0.004, //MBar
      p = 0.28, //unitless
      q = 1.34, //unitless
      c2 = shear_modulus_ref * berg_mag, // MBar * microns
      gam_ro = 1e3, //1/micro-sec /(1/micron^2) / (micron)
      wrD = 0.02, // MBar?
      inter_mat = 1.0;// unitless;
   double
      c_trap = 1.0e-3, // unitless
      c_mult = 2.5e-3, // unitless
      c_ann = 2.0e-4, //(c_mult - c_trap),//2.0e-4, // unitless
      d_ann = 6.0 * berg_mag; // microns or a variation of
#ifdef LARGE_DD
   double
      qM = 1.0e4, // 1 / microns^2
      qT = 4.0e4; // 1 / micron^2
#else
   double
      qM = 1.0e-2, // 1 / microns^2
      qT = 4.0e-2; // 1 / micron^2
#endif

   gam_ro *= 1.0 / qM;
   // This should really be fD *= sqrt(qM) / berg_mag
   // However, we're going to keep it as below so our test suite
   // stays the same as before. However, we can think of this as if
   // we scaled fD by berg_mag and then scaled things by:
   // sqrt(qM)/berg_mag
   fD *= sqrt(qM);
   std::vector<double> paramsThese {
      shear_modulus_ref, tkelv_ref, berg_mag, lbar,
      gam_ro, wrD,
      fD, c_1, tau_a, p, q, c2, inter_mat,
      c_ann, d_ann, c_trap, c_mult,
      qM, qM, qM, qM, qM, qM,
      qM, qM, qM, qM, qM, qM,
      qT, qT, qT, qT, qT, qT,
      qT, qT, qT, qT, qT, qT
   };
#ifdef STACK_PARAMS
   params.insert(params.end(), paramsThese.begin(), paramsThese.end());
#else
   kinetics.setParams(paramsThese);
#endif
}
