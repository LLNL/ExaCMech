
#if KIN_TYPE
static const int   expectedNFEvals = 23 ;
static const real8 expectedGdotVal = 0.2398180885495 ;
static const real8 expectedE2      = 0.00407276458021 ;
static const real8 expectedQ1 = 0.999687516276 ;
#else

#if XM_MUSHY
static const int   expectedNFEvals = 13 ;
static const real8 expectedGdotVal = 0.4607285834886 ;
static const real8 expectedE2      = 0.00082808734102797 ;
static const real8 expectedQ1 = 0.9956165282270 ;
#else
static const int   expectedNFEvals = 18 ;
static const real8 expectedGdotVal = 0.2475346625929 ;
static const real8 expectedE2      = 0.0009861349707681 ;
static const real8 expectedQ1 = 0.999687516276 ;
#endif

#endif
