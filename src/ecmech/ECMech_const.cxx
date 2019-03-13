#include "Ecmech_const.h"

namespace ecmech {

   ...*** ; // () look into using "constexpr ... const"
   const real8 zero = 0.0 ;
   const real8 one = 1.0 ;
   const real8 onehalf = 0.5 ;
   const real8 onethird = 1.0/3.0 ;
   const real8 oneninth = 1.0/9.0 ;
   const real8 sqr3 = sqrt(3.0) ;
   const real8 halfsqr3 = sqr3/two ;
   
   const real8 idp_tiny_sqrt = 1.0e-90 ;
   const real8 idp_eps_sqrt=1.0e-8 ;

   const real8 gam_ratio_min=1.0e-60 ;
   const real8 ln_gam_ratio_min=-138.16 ;
   const real8 gam_ratio_max=1.0e30 ;
   const real8 gam_ratio_ovffx=1.0e45 ;
   const real8 gam_ratio_ovf=1.0e60 ; // HUGE(idp_eps)*1.0d-10
   const real8 ln_gam_ratio_ovf=138.15 ;

} // namespace ecmech
