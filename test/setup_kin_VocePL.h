      {
         real8 mu = 1.0, xm = 0.01, gam_w = 1.0 ;
         real8 h0 = 200e-5, tausi = 100e-5, taus0 = 400e-5, xms = 0.05, gamss0 = 1.0e-6 ;
         std::vector<real8> paramsThese{
            mu, xm, gam_w,
            h0, tausi, taus0, xms, gamss0,
            tausi // hdn_init
         } ;
#ifdef STACK_PARAMS
         params.insert(params.end(), paramsThese.begin(), paramsThese.end()) ;
#else
         kinetics.setParams( paramsThese ) ;
#endif
      }
