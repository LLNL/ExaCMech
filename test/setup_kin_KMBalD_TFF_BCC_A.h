{
   double
      mu = 1.0,
      tK_ref = 300.,
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
      mu, tK_ref,
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
