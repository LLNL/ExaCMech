#ifndef __ECMECH_CONST_H__
#define __ECMECH_CONST_H__

#include "ECMech_port.h"
#include "ECMech_gpu_portability.h"

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

#define ISTRIDE_DEF_RATE 0
#define ISTRIDE_SPIN_V 1
#define ISTRIDE_VOL_RATIO 2
#define ISTRIDE_INT_ENG 3
#define ISTRIDE_STRESS 4
#define ISTRIDE_HISTORY 5
#define ISTRIDE_TKELV 6
#define ISTRIDE_SDD 7
#define ECMECH_NSTRIDE 8

#define MORE_DERIVS 0

namespace ecmech
{
   // We're going to use this to determine what RAJA code to run for our
   // kernels.
   enum class ExecutionStrategy { CPU, GPU, OPENMP };

   constexpr int nsvp = 7;
   constexpr int ndim = 3;
   constexpr int ne = 1;
   constexpr int nsvec = 6;
   constexpr int nsvec2 = 36;
   constexpr int nvr = 4;

   constexpr int nMiller = 4;

   // Provide indices for the matModel stride array so codes outside of the library
   // can use them.
   constexpr int istride_def_rate = ISTRIDE_DEF_RATE;
   constexpr int istride_spin_v = ISTRIDE_SPIN_V;
   constexpr int istride_vol_ratio = ISTRIDE_VOL_RATIO;
   constexpr int istride_int_eng = ISTRIDE_INT_ENG;
   constexpr int istride_stress = ISTRIDE_STRESS;
   constexpr int istride_history = ISTRIDE_HISTORY;
   constexpr int istride_tkelv = ISTRIDE_TKELV;
   constexpr int istride_sdd = ISTRIDE_SDD;
   constexpr int nstride = ECMECH_NSTRIDE;

   constexpr int ntvec = 5;
   constexpr int nwvec = 3;
   constexpr int qdim = 4;
   constexpr int invdim = 4;
   constexpr int emapdim = 3;

   constexpr int iSvecS = nsvec - 1; // index like SVEC in F90 coding
   constexpr int iSvecP = nsvec;

   // indexing into array of outputs
   constexpr int i_sdd_bulk = 0;
   constexpr int i_sdd_gmod = 1;
   constexpr int nsdd = 2;

   constexpr int i_ne_total = 0;

   constexpr double zero = 0.0;
   constexpr double one = 1.0;
   constexpr double two = 2.0;
   constexpr double three = 3.0;
   constexpr double six = 6.0;
   constexpr double onehalf = 0.5;
   constexpr double onethird = 1.0 / 3.0;
   constexpr double oneninth = 1.0 / 9.0;
   constexpr double oneqrtr = 0.25;
   constexpr double thrhalf = 1.5;
   constexpr double fourthirds = 4.0 / 3.0;
   constexpr double twothird = 2.0 / 3.0;

   constexpr double sqr2 = M_SQRT2;
   constexpr double sqr3 = M_SQRT3;
   constexpr double sqr3b2 = M_SQRT3B2;
   constexpr double sqr2i = M_SQRT2I;
   constexpr double sqr3i = M_SQRT3I;
   constexpr double sqr6 = M_SQRT6;
   constexpr double sqr6i = M_SQRT6I;
   constexpr double sqr2b3 = M_SQRT2B3;
   constexpr double halfsqr3 = M_HALFSQRT3;
   constexpr double twosqr3 = M_TWOSQRT3;

   constexpr double idp_tiny_sqrt = 1.0e-90;
   constexpr double idp_eps_sqrt = 1.0e-8;
   constexpr double idp_eps = 2.0e-16;


   constexpr double gam_ratio_min = 1.0e-60;
   constexpr double ln_gam_ratio_min = -138.16;
   constexpr double gam_ratio_max = 1.0e30;
   constexpr double gam_ratio_ovffx = 1.0e45;
   constexpr double gam_ratio_ovf = 1.0e60; // HUGE(idp_eps)*1.0d-10
   constexpr double ln_gam_ratio_ovf = 138.15;

   // as in evptn and evptnconst
   constexpr double st_toler = 1.0e-11;
   constexpr double epsdot_scl_nzeff = idp_eps_sqrt;
   constexpr double e_scale = 5e-4;
   constexpr double r_scale = 0.01;
   constexpr int st_max_iter = 200;
} // (namespace ecmech)

#define ECM_EXEC_STRAT_CPU    ecmech::ExecutionStrategy::CPU
#define ECM_EXEC_STRAT_GPU    ecmech::ExecutionStrategy::GPU
#define ECM_EXEC_STRAT_OPENMP ecmech::ExecutionStrategy::OPENMP

#endif // __ECMECH_CONST_H__
