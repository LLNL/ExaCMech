{

   double
      mu_ref = 2.7e10, // mbar
      tK_ref = 300., // k
      berg_mag = 1.0e-10, //microns
      lbar = 10.0 * berg_mag, // microns
      shear_speed = std::sqrt(mu_ref / 2703), //m/s -> sqrt(mu/rho)
      fD = 4.0e12, //Hz or variation of
      c_1 = 0.65 * berg_mag * berg_mag * berg_mag / 1.3806504e-23, // g_0 * b^3 / kB
      tau_a = 1.0e3,// Pa or a variation?
      p = 0.5,//unitless
      q = 2.0,//unitless
      tau_0 = 1.0e7,//Pa or a variation?
      c2 = mu_ref * berg_mag, //Pa * length
      inter_mat = 1.0;// unitless
   double beta_0 = (3.0 * 1.3806504e-23 * tK_ref * 4.0) / (20.0 * shear_speed * berg_mag * berg_mag);
   double c_3 = beta_0 * shear_speed / berg_mag; //?
   double
      c_ann = 1.5, // unitless
      d_ann = 6.0 * berg_mag, // microns or a variation of
      c_trap = 0.5, // unitless
      c_mult = 1.2; // unitless
#ifdef LARGE_DD
   double
      qM = 1.0e16, // 1 / microns^2
      qT = 4.0e16; // 1 / micron^2
#else
   double
      qM = 1.0e10, // 1 / microns^2
      qT = 4.0e10; // 1 / micron^2
#endif
   std::vector<double> paramsThese {
      mu_ref, tK_ref, berg_mag, lbar,
      shear_speed, c_3,
      fD, c_1, tau_a, p, q, tau_0, c2, inter_mat,
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
