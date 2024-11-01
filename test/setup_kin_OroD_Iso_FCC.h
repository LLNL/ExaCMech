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
