#ifndef __ECMECH_CONST_H__
#define __ECMECH_CONST_H__

#include <math.h>
#include "ECMech_port.h"

/* [J/deg. K] */
#ifndef M_BOLTZ_J_K
#define M_BOLTZ_J_K 1.3806504e-23 
#endif

/* per mole */
#ifndef M_AVOGADRO
#define M_AVOGADRO 6.02214179e23
#endif

#define MORE_DERIVS 0

#define ECMECH_INFO_LEVEL_INFO       0
#define ECMECH_INFO_LEVEL_WARN       1
#define ECMECH_INFO_LEVEL_IMPORTANT  2
#define ECMECH_INFO_LEVEL_CRITICAL   3

namespace ecmech
{

const int nsvp         =  7  ;
const int ndim         =  3  ;
const int ne           =  2  ;
const int nsvec        =  6  ;
const int nsvec2       = 36  ;
const int nvr          =  4  ;

const int ntvec        =  5  ;
const int nwvec        =  3  ;
const int qdim         =  4  ;
const int invdim       =  4  ;
const int emapdim      =  3  ;

const int iSvecS       =  nsvec-1  ; // index like SVEC in F90 coding

// indexing into array of outputs 
const int i_sdd_bulk   = 0  ;
const int i_sdd_gmod   = 1  ;
const int nsdd         = 2  ;

const int i_ne_total   = 0 ;
const int i_ne_cold    = 1 ;
// const int i_ne_melt    = I_NE_MELT ; // no -- do not require this as an output -- generally not needed

const int info_level_info      = ECMECH_INFO_LEVEL_INFO     ;
const int info_level_warn      = ECMECH_INFO_LEVEL_WARN     ;
const int info_level_important = ECMECH_INFO_LEVEL_IMPORTANT;
const int info_level_critical  = ECMECH_INFO_LEVEL_CRITICAL ;

const real8 zero = 0.0 ;
const real8 one = 1.0 ;
const real8 two = 2.0 ;
const real8 three = 3.0 ;
const real8 six = 6.0 ;
const real8 onehalf = 0.5 ;
const real8 onethird = 1.0/3.0 ;
const real8 oneninth = 1.0/9.0 ;
const real8 oneqrtr = 0.25 ;
const real8 thrhalf = 1.5 ;
const real8 fourthirds = 4.0/3.0 ;
const real8 twothird   = 2.0/3.0 ;

const real8 sqr2       = sqrt(two) ;
const real8 sqr3       = sqrt(three) ;
const real8 sqr3b2     = sqrt(thrhalf) ;
const real8 sqr2i      = one/sqr2 ;
const real8 sqr3i      = one/sqr3 ;
const real8 sqr6       = sqrt(six) ;
const real8 sqr6i      = one/sqr6 ;
const real8 sqr2b3     = one/sqr3b2 ;
const real8 halfsqr3   = sqr3/two ;
const real8 twosqr3    = two*sqr3 ;

const real8 idp_tiny_sqrt = 1.0e-90 ;
const real8 idp_eps_sqrt=1.0e-8 ;

const real8 gam_ratio_min=1.0e-60 ;
const real8 ln_gam_ratio_min=-138.16 ;
const real8 gam_ratio_max=1.0e30 ;
const real8 gam_ratio_ovffx=1.0e45 ;
const real8 gam_ratio_ovf=1.0e60 ; // HUGE(idp_eps)*1.0d-10
const real8 ln_gam_ratio_ovf=138.15 ;

}  // (namespace ecmech)

#endif  // __ECMECH_CONST_H__
