#if XM_MUSHY
#define XM_VAL 0.1
#else
#define XM_VAL 0.01
#endif

{
   double mu = 1.0, xm = XM_VAL, gam_w = 1.0;
   double h0 = 200e-5, tausi = 100e-5, taus0 = 400e-5, xms = 0.05, gamss0 = 1.0e-6;
#ifdef KIN_NONLINEAR
   double xmprime = 2.0;
   std::vector<double> paramsThese {
      mu, xm, gam_w,
      h0, tausi, taus0, xmprime, xms, gamss0,
      tausi // hdn_init
   };
#else
   std::vector<double> paramsThese {
      mu, xm, gam_w,
      h0, tausi, taus0, xms, gamss0,
      tausi // hdn_init
   };
#endif
#ifdef STACK_PARAMS
   params.insert(params.end(), paramsThese.begin(), paramsThese.end());
#else
   kinetics.setParams(paramsThese);
#endif
}
