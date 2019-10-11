#ifndef __ECMECH_CONST_H__
#define __ECMECH_CONST_H__

// #include <math.h>
#include "ECMech_port.h"
#include "ECMech_cuda_portability.h"

/* [J/deg. K] */
#ifndef M_BOLTZ_J_K
#define M_BOLTZ_J_K 1.3806504e-23
#endif

/* per mole */
#ifndef M_AVOGADRO
#define M_AVOGADRO 6.02214179e23
#endif

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880168872421
#endif

#ifndef M_SQRT3
#define M_SQRT3 1.73205080756887729352744634151
#endif

#ifndef M_SQRT3B2
#define M_SQRT3B2 1.22474487139158904909864203735
#endif

#ifndef M_SQRT2I
#define M_SQRT2I 0.707106781186547524400844362105
#endif

#ifndef M_SQRT3I
#define M_SQRT3I 0.577350269189625764509148780501
#endif

#ifndef M_SQRT6
#define M_SQRT6 2.44948974278317809819728407471
#endif

#ifndef M_SQRT6I
#define M_SQRT6I 0.408248290463863016366214012450
#endif

#ifndef M_SQRT2B3
#define M_SQRT2B3 0.816496580927726032732428024904
#endif

#ifndef M_HALFSQRT3
#define M_HALFSQRT3 0.866025403784438646763723170755
#endif

#ifndef M_TWOSQRT3
#define M_TWOSQRT3 3.46410161513775458705489268302
#endif

#define MORE_DERIVS 0

// #include <limits>

namespace ecmech
{
   // We're going to use this to determine what RAJA code to run for our
   // kernels.
   enum class Accelerator { CPU, CUDA, OPENMP };

   const int nsvp = 7;
   const int ndim = 3;
   const int ne = 1;
   const int nsvec = 6;
   const int nsvec2 = 36;
   const int nvr = 4;

   const int ntvec = 5;
   const int nwvec = 3;
   const int qdim = 4;
   const int invdim = 4;
   const int emapdim = 3;

   const int iSvecS = nsvec - 1; // index like SVEC in F90 coding
   const int iSvecP = nsvec;

   // indexing into array of outputs
   const int i_sdd_bulk = 0;
   const int i_sdd_gmod = 1;
   const int nsdd = 2;

   const int i_ne_total = 0;
   // const int i_ne_cold    = ? ; // generally not needed
   // const int i_ne_melt    = ? ; // generally not needed

   const double zero = 0.0;
   const double one = 1.0;
   const double two = 2.0;
   const double three = 3.0;
   const double six = 6.0;
   const double onehalf = 0.5;
   const double onethird = 1.0 / 3.0;
   const double oneninth = 1.0 / 9.0;
   const double oneqrtr = 0.25;
   const double thrhalf = 1.5;
   const double fourthirds = 4.0 / 3.0;
   const double twothird = 2.0 / 3.0;

   const double sqr2 = M_SQRT2;
   const double sqr3 = M_SQRT3;
   const double sqr3b2 = M_SQRT3B2;
   const double sqr2i = M_SQRT2I;
   const double sqr3i = M_SQRT3I;
   const double sqr6 = M_SQRT6;
   const double sqr6i = M_SQRT6I;
   const double sqr2b3 = M_SQRT2B3;
   const double halfsqr3 = M_HALFSQRT3;
   const double twosqr3 = M_TWOSQRT3;

   const double idp_tiny_sqrt = 1.0e-90;
   const double idp_eps_sqrt = 1.0e-8;

   const double gam_ratio_min = 1.0e-60;
   const double ln_gam_ratio_min = -138.16;
   const double gam_ratio_max = 1.0e30;
   const double gam_ratio_ovffx = 1.0e45;
   const double gam_ratio_ovf = 1.0e60; // HUGE(idp_eps)*1.0d-10
   const double ln_gam_ratio_ovf = 138.15;

   // as in evptn and evptnconst
   const double st_toler = 1.0e-11;
   const double epsdot_scl_nzeff = idp_eps_sqrt;
   const double e_scale = 5e-4;
   const double r_scale = 0.01;
   const int st_max_iter = 200;
} // (namespace ecmech)

#endif // __ECMECH_CONST_H__
