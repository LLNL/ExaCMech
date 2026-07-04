/**
 * @file ECMech_const.h
 * @brief Central header for physical/mathematical constants, indexing conventions, and
 * solver tolerances used throughout ExaCMech.
 *
 * This header has no dependencies on any crystal-plasticity-specific data structures; it
 * is included (directly or transitively via ECMech_core.h) by nearly every other header
 * in the library. It provides:
 *
 * - **Physical constants**: Boltzmann's constant, Avogadro's number
 * - **Math constants**: precomputed square roots and related irrational factors, used so
 *   that expressions like `1/sqrt(3)` are not recomputed at every call site
 * - **Stride/indexing conventions**: the fixed layout of the flattened arrays that
 *   `matModelBase::getResponse` (see ECMech_matModelBase.h) uses to communicate strain
 *   rate, spin, stress, history, temperature, and derived-derivative data to and from a
 *   calling code
 * - **Dimensionality constants**: fixed sizes for the various tensor representations used
 *   across the library (see ECMech_util.h for the representations themselves)
 * - **Named numeric literals**: `zero`, `one`, `two`, ... used in place of bare literals
 *   for readability and to guarantee consistent `double` typing in mixed-precision
 *   expressions
 * - **Solver tolerances**: cutoffs used by the slip-kinetics power-law evaluations (see
 *   the `kinetics/` model headers) and by the nonlinear state-update solvers to avoid
 *   overflow/underflow and to bound iteration counts
 *
 * @see ECMech_core.h for the umbrella header that pulls this in along with
 * ECMech_gpu_portability.h and ECMech_port.h
 * @see ECMech_matModelBase.h for how the `istride_*`/`nstride` values are used
 */

#ifndef __ECMECH_CONST_H__
#define __ECMECH_CONST_H__

#include "ECMech_port.h"
#include "ECMech_gpu_portability.h"

/** @brief Boltzmann's constant, in J/K. */
/* [J/deg. K] */
#ifndef M_BOLTZ_J_K
#define M_BOLTZ_J_K 1.3806504e-23
#endif

/** @brief Avogadro's number, per mole. */
/* per mole */
#ifndef M_AVOGADRO
#define M_AVOGADRO 6.02214179e23
#endif

/**
 * @brief Precomputed irrational math constants.
 *
 * These are `#define`-guarded so that a value already provided by `<cmath>` or another
 * library (e.g. `M_SQRT2`) is not clobbered, and so that ExaCMech still has a definition
 * available on platforms where the C library does not provide the corresponding POSIX
 * math constant.
 *
 * Naming convention: `M_SQRTn` is √n, `M_SQRTnI` is 1/√n, `M_SQRTnBm` is √(n/m).
 */
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

/**
 * @brief Indices into the flattened stride array passed to/from
 * `matModelBase::getResponse`, one per field it communicates.
 *
 * `matModelBase::getResponse` (ECMech_matModelBase.h) communicates its inputs and
 * outputs as a set of raw pointers plus a stride array of length #ECMECH_NSTRIDE. Each
 * `ISTRIDE_*` macro is the index into that stride array for one field:
 *
 * - `ISTRIDE_DEF_RATE`: deformation rate (rate of deformation tensor, deviatoric 5-vector)
 * - `ISTRIDE_SPIN_V`: spin vector (3-component axial vector)
 * - `ISTRIDE_VOL_RATIO`: relative volume ratios
 * - `ISTRIDE_INT_ENG`: internal energy
 * - `ISTRIDE_STRESS`: Cauchy stress (deviatoric + pressure, 7-vector)
 * - `ISTRIDE_HISTORY`: per-point history/state variables
 * - `ISTRIDE_TKELV`: temperature in Kelvin
 * - `ISTRIDE_SDD`: "state-dependent data" outputs (e.g. bulk modulus, shear modulus)
 *
 * The `ecmech::istride_*` / `ecmech::nstride` constants below mirror these macros as
 * proper `constexpr int`s for use in C++ code; the macros exist mainly so external codes
 * that only need the raw values can pick them up without linking against the `ecmech`
 * namespace.
 */
#define ISTRIDE_DEF_RATE 0
#define ISTRIDE_SPIN_V 1
#define ISTRIDE_VOL_RATIO 2
#define ISTRIDE_INT_ENG 3
#define ISTRIDE_STRESS 4
#define ISTRIDE_HISTORY 5
#define ISTRIDE_TKELV 6
#define ISTRIDE_SDD 7
/** @brief Total number of fields in the `matModelBase::getResponse` stride array. */
#define ECMECH_NSTRIDE 8

/**
 * @brief Reserved for requesting additional derivative output beyond what is currently
 * implemented; not presently consumed by any code path.
 */
#define MORE_DERIVS 0

namespace ecmech
{
   /**
    * @brief Selects which RAJA execution back end a kernel dispatch should target.
    *
    * Used by call sites that need to pick a RAJA policy (sequential/CPU, CUDA/HIP GPU,
    * or OpenMP) at runtime rather than compile time.
    */
   enum class ExecutionStrategy {
      CPU,   ///< Run using RAJA's sequential/CPU execution policy.
      GPU,   ///< Run using RAJA's CUDA or HIP execution policy.
      OPENMP ///< Run using RAJA's OpenMP execution policy.
   };

   /**
    * @brief Fixed dimensionality constants for the tensor representations used across
    * ExaCMech (see ECMech_util.h for the representations themselves).
    *
    * - `nsvp` : size of the extended stress-like 7-vector (deviatoric 6 + pressure), as
    *   used for some legacy/alternate stress storage
    * - `ndim` : spatial dimension (always 3)
    * - `ne`   : number of internal-energy values carried per point (currently always 1)
    * - `nsvec` : size of the Voigt symmetric 6-vector representation
    * - `nsvec2` : `nsvec * nsvec`, size of a flattened 6x6 matrix (e.g. a stiffness
    *   tensor in Voigt notation)
    * - `nvr` : number of relative-volume values carried per point: `[rel_vol_n, rel_vol_{n+1},
    *   (rel_vol_{n+1}-rel_vol_n)/dt, rel_vol_{n+1}-rel_vol_n]`, where `rel_vol` is the relative volume (volume ratio
    *   to the reference state)
    */
   constexpr int nsvp = 7;
   constexpr int ndim = 3;
   constexpr int ne = 1;
   constexpr int nsvec = 6;
   constexpr int nsvec2 = 36;
   constexpr int nvr = 4;

   /** @brief Number of Miller indices used to specify a crystallographic plane/direction (h,k,i,l for hexagonal, or truncated for cubic). */
   constexpr int nMiller = 4;

   // Provide indices for the matModel stride array so codes outside of the library
   // can use them.
   /** @brief Index of the deformation-rate field in the `getResponse` stride array; see #ISTRIDE_DEF_RATE. */
   constexpr int istride_def_rate = ISTRIDE_DEF_RATE;
   /** @brief Index of the spin-vector field in the `getResponse` stride array; see #ISTRIDE_SPIN_V. */
   constexpr int istride_spin_v = ISTRIDE_SPIN_V;
   /** @brief Index of the relative-volume field in the `getResponse` stride array; see #ISTRIDE_VOL_RATIO. */
   constexpr int istride_vol_ratio = ISTRIDE_VOL_RATIO;
   /** @brief Index of the internal-energy field in the `getResponse` stride array; see #ISTRIDE_INT_ENG. */
   constexpr int istride_int_eng = ISTRIDE_INT_ENG;
   /** @brief Index of the stress field in the `getResponse` stride array; see #ISTRIDE_STRESS. */
   constexpr int istride_stress = ISTRIDE_STRESS;
   /** @brief Index of the history/state field in the `getResponse` stride array; see #ISTRIDE_HISTORY. */
   constexpr int istride_history = ISTRIDE_HISTORY;
   /** @brief Index of the temperature field in the `getResponse` stride array; see #ISTRIDE_TKELV. */
   constexpr int istride_tkelv = ISTRIDE_TKELV;
   /** @brief Index of the state-dependent-data field in the `getResponse` stride array; see #ISTRIDE_SDD. */
   constexpr int istride_sdd = ISTRIDE_SDD;
   /** @brief Total number of fields in the `getResponse` stride array; see #ECMECH_NSTRIDE. */
   constexpr int nstride = ECMECH_NSTRIDE;

   /** @brief Size of the deviatoric 5-vector tensor representation (traceless symmetric tensor). */
   constexpr int ntvec = 5;
   /** @brief Size of the spin/axial 3-vector representation of a skew-symmetric tensor. */
   constexpr int nwvec = 3;
   /** @brief Size of a unit quaternion (4 components) used for crystal orientation. */
   constexpr int qdim = 4;
   /** @brief Dimension of the inverse/derivative bookkeeping arrays associated with a quaternion (matches #qdim). */
   constexpr int invdim = 4;
   /** @brief Dimension used when mapping between exponential-map and quaternion orientation representations. */
   constexpr int emapdim = 3;

   /** @brief Index of the last deviatoric component within a #nsvec-sized Voigt array (i.e. `nsvec - 1`); mirrors the Fortran `SVEC` convention. */
   constexpr int iSvecS = nsvec - 1; // index like SVEC in F90 coding
   /** @brief Index of the pressure/trace component appended after a #nsvec-sized Voigt array. */
   constexpr int iSvecP = nsvec;

   // indexing into array of outputs
   /** @brief Index of the bulk modulus within the "state-dependent data" (SDD) output array. */
   constexpr int i_sdd_bulk = 0;
   /** @brief Index of the shear modulus within the "state-dependent data" (SDD) output array. */
   constexpr int i_sdd_gmod = 1;
   /** @brief Total number of entries in the "state-dependent data" (SDD) output array. */
   constexpr int nsdd = 2;

   /** @brief Index of the total-strain-energy entry when accumulating per-element energy output. */
   constexpr int i_ne_total = 0;

   /**
    * @brief Named floating-point literals used in place of bare numeric literals.
    *
    * These exist purely for readability (`onethird` reads better than `1.0 / 3.0`
    * scattered through tensor-math expressions) and to guarantee the literals are typed
    * as `double` and computed once rather than re-derived at every use site.
    */
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

   /**
    * @brief `constexpr double` mirrors of the `M_SQRT*` macros above, for use directly
    * in C++ expressions within the `ecmech` namespace.
    */
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

   /** @brief Floor below which a squared quantity is treated as exactly zero, to avoid `sqrt`/divide-by-zero issues from roundoff. */
   constexpr double idp_tiny_sqrt = 1.0e-90;
   /** @brief Small tolerance (~sqrt of double-precision epsilon) used to switch between an exact and a linearized/small-argument approximation of an expression. */
   constexpr double idp_eps_sqrt = 1.0e-8;
   /** @brief Double-precision machine epsilon used directly (as opposed to its square root, #idp_eps_sqrt) for tight tolerance checks. */
   constexpr double idp_eps = 2.0e-16;


   /**
    * @brief Bounds on the ratio of resolved shear stress to flow strength (or similar
    * dimensionless "how far into the power law" ratios) used by the slip-kinetics power
    * laws in `kinetics/` to keep `pow()`/`exp()`/`log()` evaluations well-behaved.
    *
    * - `gam_ratio_min` / `ln_gam_ratio_min`: below this ratio (or its log), the slip rate
    *   contribution is treated as exactly zero rather than evaluating a power law on a
    *   near-zero base
    * - `gam_ratio_max` / `gam_ratio_ovf` / `ln_gam_ratio_ovf`: above these ratios, the
    *   power-law kinetics would overflow, so the slip rate is instead capped using
    *   `gam_ratio_ovffx`
    * - `gam_ratio_ovffx`: the large-but-finite multiplier used to report an "essentially
    *   infinite" slip rate once the ratio exceeds `gam_ratio_ovf`, so that downstream
    *   solvers see a large finite number (signaling the need to cut the time step) rather
    *   than `inf`/`nan`
    */
   constexpr double gam_ratio_min = 1.0e-60;
   constexpr double ln_gam_ratio_min = -138.16;
   constexpr double gam_ratio_max = 1.0e30;
   constexpr double gam_ratio_ovffx = 1.0e45;
   constexpr double gam_ratio_ovf = 1.0e60; // HUGE(idp_eps)*1.0d-10
   constexpr double ln_gam_ratio_ovf = 138.15;

   // as in evptn and evptnconst
   /** @brief Convergence tolerance for the implicit elastic-strain/lattice-rotation (evptn) nonlinear solve. */
   constexpr double st_toler = 1.0e-11;
   /** @brief Scale factor (set to #idp_eps_sqrt) used to non-dimensionalize the effective strain rate when it is used as a solver scaling quantity. */
   constexpr double epsdot_scl_nzeff = idp_eps_sqrt;
   /** @brief Characteristic strain increment used to scale the elastic-strain residual in the evptn nonlinear solve. */
   constexpr double e_scale = 5e-4;
   /** @brief Characteristic rotation increment used to scale the lattice-orientation residual in the evptn nonlinear solve. */
   constexpr double r_scale = 0.01;
   /** @brief Maximum iteration count allowed for the evptn nonlinear solve before it is treated as a failed step. */
   constexpr int st_max_iter = 200;
} // (namespace ecmech)

/** @brief Convenience macro form of ecmech::ExecutionStrategy::CPU for use outside the `ecmech` namespace. */
#define ECM_EXEC_STRAT_CPU    ecmech::ExecutionStrategy::CPU
/** @brief Convenience macro form of ecmech::ExecutionStrategy::GPU for use outside the `ecmech` namespace. */
#define ECM_EXEC_STRAT_GPU    ecmech::ExecutionStrategy::GPU
/** @brief Convenience macro form of ecmech::ExecutionStrategy::OPENMP for use outside the `ecmech` namespace. */
#define ECM_EXEC_STRAT_OPENMP ecmech::ExecutionStrategy::OPENMP

#endif // __ECMECH_CONST_H__
