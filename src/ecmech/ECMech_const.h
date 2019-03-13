#ifndef __ECMECH_CONST_H__
#define __ECMECH_CONST_H__

#include <math.h>

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

}  // (namespace ecmech)

#endif  // __ECMECH_CONST_H__
