{
   double
      mu_ref = 1.0, //MBar
      tK_ref = 300., //K
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
      c2 = mu_ref * berg_mag, // MBar * microns
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

   gam_ro *= 1.0 / (berg_mag * qM);
   fD *= sqrt(berg_mag * qM);

   std::vector<double> paramsThese {
      mu_ref, tK_ref, berg_mag, lbar,
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
