   {
      real8
         mu     = 1.0,
         tK_ref = 300.,
         c_1    = 20000.,
         tau_a  = 0.004,
         p      = 0.28,
         q      = 1.34,
         gam_wo = 20.,
         gam_ro = 1e3,
         wrD    = 0.02,
         go     = 10e-5,
         s      = 5e-5 ;
      real8
         k1       = 100.0,
         k2o      =  10.0,
         ninv     =  20.0,
         gamma_o  =  1e-6 ;
      real8
         rho_dd_init = 0.25 ;
      std::vector<real8> params{
         mu, tK_ref, c_1, tau_a, p, q, gam_wo, gam_ro, wrD, go, s,
         k1, k2o, ninv, gamma_o,
         rho_dd_init
      } ;
      kinetics.setParams( params ) ;
   }
