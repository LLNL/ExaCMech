/**
 * @file ECMech_util.h
 * @brief Core mathematical utility library for crystal plasticity tensor operations.
 * 
 * This header provides a comprehensive suite of mathematical utilities for crystal plasticity
 * finite deformation mechanics, including vector/matrix operations, tensor transformations,
 * quaternion manipulations, frame rotations, and various specialized conversions between
 * different tensor representations.
 * 
 * Functionality categories:
 * - **Basic vector operations**: Addition, scaling, dot products, norms, normalization
 * - **Matrix-vector operations**: Matrix-vector products, transposed operations
 * - **Matrix-matrix operations**: Products, transposes, outer products
 * - **Tensor decompositions**: Symmetric/skew decomposition, deviatoric extraction
 * - **Tensor conversions**: Between full 3×3, deviatoric 5-vector, Voigt 6-vector, extended 7-vector
 * - **Quaternion operations**: Conversions, products, inverses, exponential map
 * - **Rotation tensors**: Quaternion to rotation matrix, 5×5 rotation matrices for deviatoric tensors
 * - **Rotation derivatives**: Sensitivities for implicit Jacobians in crystal plasticity
 * - **Miller index conversions**: Crystallographic notation to Cartesian coordinates
 * - **Debug utilities**: Formatted printing for vectors and matrices
 * 
 * Tensor representation conventions:
 * 
 * 1. **Full 3×3 tensors**: Stored in row-major order using ECMECH_NN_INDX macro
 * 2. **Deviatoric 5-vector** (traceless symmetric): Uses √2 and √3 normalization for norm preservation
 * 3. **Voigt 6-vector** (general symmetric): Standard [11, 22, 33, 23, 31, 12] ordering
 * 4. **Extended 7-vector** (deviatoric + pressure): [dev11', dev22', dev33', dev23, dev31, dev12, -p]
 * 5. **Spin 3-vector** (skew-symmetric): Axial vector form [w1, w2, w3]
 * 6. **Quaternions**: 4-component unit vectors [q0, q1, q2, q3] representing rotations
 * 7. **Exponential map**: 3-component angle-axis representation
 * 
 * Index ordering:
 * - Matrices: Row-major storage via ECMECH_NN_INDX and ECMECH_NM_INDX macros
 * - Voigt stress/strain: [11, 22, 33, 23=yz, 31=zx, 12=xy]
 * - Deviatoric: Component 0-4 (5 DOF), optional component 5 for trace/pressure
 * 
 * Frame conventions:
 * - **Sample frame**: Global/laboratory/spatial, fixed in space
 * - **Lattice frame**: Crystal/material, rotates with crystal orientation
 * - Transformations via rotation matrices C or quaternions q
 * 
 * Mathematical notation (in comments):
 * - Vectors: lowercase bold or component form
 * - Matrices/tensors: uppercase bold or index form
 * - Deviatoric: prime superscript (e.g., σ')
 * - Sample frame: subscript "sm" or no subscript
 * - Lattice frame: subscript "lat" or "xtal"
 * 
 * Design philosophy:
 * - Template parameters for dimensions enable compile-time optimization
 * - Inline functions minimize call overhead
 * - __ecmech_hdev__ decorator enables GPU execution where supported
 * - Minimal dynamic memory allocation (stack-based operations)
 * - Separate read-only inputs (const) from writable outputs
 * - Output arguments first in parameter lists for consistency
 * 
 * Performance considerations:
 * - All operations are inlined for optimization
 * - Template dimensions allow loop unrolling
 * - Row-major storage matches SNLS solver conventions
 * - RAJA::View used for multi-dimensional array indexing in complex functions
 * 
 * Usage with util headers:
 * Many rotation derivative functions include util/.h headers (mc_vars_set.h, 
 * vad_vars_set.h, etc.) to extract tensor components into named variables,
 * improving readability of complex mathematical expressions.
 * 
 * @ingroup ECMech_core_utilities
 * 
 * @see ECMech_const.h for mathematical constants (√2, √3, etc.)
 * @see util/mc_vars_set.h for rotation matrix component extraction
 * @see util/vad_vars_set.h for deviatoric vector components (sample frame)
 * @see util/vadl_vars_set.h for deviatoric vector components (lattice frame)
 */
// -*-c++-*-

#ifndef ECMECH_UTIL_H
#define ECMECH_UTIL_H

#include "ECMech_core.h"

#include <cmath>

#if defined(ECMECH_DEBUG) && defined(__ecmech_host_only__)
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#endif

#include "RAJA/RAJA.hpp"

/**
 * @defgroup array_indexing Array Indexing Macros
 * @brief Macros for row-major array indexing of matrices.
 * 
 * ECMech uses row-major storage for all matrices to maintain consistency with
 * the SNLS solver library and RAJA execution framework. These macros compute
 * linear array indices from (row, column) pairs for square and non-square matrices.
 * 
 * Row-major storage pattern for 3×3 matrix:
 * ```
 * A = | A00  A01  A02 |    Stored as: [A00, A01, A02, A10, A11, A12, A20, A21, A22]
 *     | A10  A11  A12 |    Index:      [ 0,   1,   2,   3,   4,   5,   6,   7,   8 ]
 *     | A20  A21  A22 |
 * ```
 * 
 * @note Column-major storage (Fortran-style) would be: index = row + col*nrow
 * @note These macros are candidates for replacement with RAJA::View in future updates
 * 
 * @{
 */

/**
 * @brief Compute linear index for element (row, col) of square n×n matrix in row-major order.
 * @param p Row index (0-based)
 * @param q Column index (0-based)
 * @param nDim Matrix dimension (number of rows = number of columns)
 * @return Linear array index
 * 
 * Formula: index = row × ncols + col = p × nDim + q
 * 
 * Example: For 3×3 matrix element A[1][2] → ECMECH_NN_INDX(1,2,3) = 1×3+2 = 5
 */
#define ECMECH_NN_INDX(p, q, nDim) (p) * (nDim) + (q)

/**
 * @brief Compute linear index for element (row, col) of non-square p×q matrix in row-major order.
 * @param p Row index (0-based)
 * @param q Column index (0-based)
 * @param pDim Number of rows
 * @param qDim Number of columns
 * @return Linear array index
 * 
 * Formula: index = row × ncols + col = p × qDim + q
 * 
 * Example: For 5×3 matrix element A[2][1] → ECMECH_NM_INDX(2,1,5,3) = 2×3+1 = 7
 * 
 * Common usage: Indexing into Schmid tensor arrays (ntvec × nslip matrices)
 */
#define ECMECH_NM_INDX(p, q, pDim, qDim) (p) * (qDim) + (q)

/** @} */ // end of array_indexing group

namespace ecmech {

/**
 * @defgroup vector_operations Vector Operations
 * @brief Fundamental vector arithmetic and properties.
 * 
 * Template-based vector operations with compile-time dimensions for optimal performance.
 * All functions operate on contiguous arrays of doubles with specified compile-time or
 * runtime lengths.
 * 
 * Conventions:
 * - Output parameters first (non-const pointers)
 * - Input parameters second (const pointers)
 * - Template parameter `n` specifies vector dimension
 * - Runtime variant `_n` suffix for dynamic dimensions
 * 
 * @{
 */

   /**
    * @brief Vector addition: v = a + b (element-wise).
    * @tparam n Vector dimension
    * @param[out] v Result vector (v[i] = a[i] + b[i])
    * @param[in] a First input vector
    * @param[in] b Second input vector
    * 
    * Computes element-wise sum of two vectors.
    * 
    * Mathematical operation: v_i = a_i + b_i for i = 0, ..., n-1
    */
   template<int n>
   __ecmech_hdev__
   inline void vecsVapb(double* const v,
                        const double* const a,
                        const double* const b) {
      for (int i = 0; i<n; ++i) {
         v[i] = a[i] + b[i];
      }
   }

   /**
    * @brief Element-wise (Hadamard) product: v = a ⊙ b.
    * @tparam n Vector dimension
    * @param[out] v Result vector (v[i] = a[i] × b[i])
    * @param[in] a First input vector
    * @param[in] b Second input vector
    * 
    * Computes element-wise product (not dot product or outer product).
    * Also known as Hadamard product or Schur product.
    * 
    * Mathematical operation: v_i = a_i · b_i for i = 0, ..., n-1
    * 
    * @note "Adiag" in function name suggests diagonal matrix interpretation:
    *       If A = diag(a), then v = A·b
    */
   template<int n>
   __ecmech_hdev__
   inline void vecsVAdiagB(double* const v,
                           const double* const a,
                           const double* const b) {
      for (int i = 0; i<n; ++i) {
         v[i] = a[i] * b[i];
      }
   }

   /**
    * @brief Dot product (inner product): y = a · b.
    * @tparam n Vector dimension
    * @param[in] a First vector
    * @param[in] b Second vector
    * @return Scalar dot product value
    * 
    * Computes standard Euclidean inner product.
    * 
    * Mathematical operation: y = Σ_{i=0}^{n-1} a_i · b_i
    * 
    * Common usage: Resolving stress onto slip systems, computing strain energy
    */
   template<int n>
   __ecmech_hdev__
   inline double vecsyadotb(const double* const a,
                            const double* const b) {
      double y = 0.0;
      for (int i = 0; i<n; ++i) {
         y += a[i] * b[i];
      }

      return y;
   }

   /**
    * @brief Sum of absolute values: s = Σ|a_i|.
    * @tparam n Vector dimension
    * @param[in] a Input vector
    * @return Sum of absolute values (L1 norm)
    * 
    * Computes the L1 (Manhattan) norm of the vector.
    * 
    * Mathematical operation: s = Σ_{i=0}^{n-1} |a_i|
    * 
    * Common usage: Checking convergence criteria, measuring slip rate totals
    */
   template<int n>
   __ecmech_hdev__
   inline double vecsssumabs(const double* const a) {
      double s = 0.0;
      for (int i = 0; i<n; ++i) {
         s += fabs(a[i]);
      }

      return s;
   }

   /**
    * @brief Sum of vector components: s = Σa_i.
    * @tparam n Vector dimension
    * @param[in] a Input vector
    * @return Sum of all components
    * 
    * Mathematical operation: s = Σ_{i=0}^{n-1} a_i
    * 
    * @note For symmetric tensors in Voigt form, this does NOT give the trace
    *       due to the √2 factors on off-diagonal components
    */
   template<int n>
   __ecmech_hdev__
   inline double vecsssum(const double* const a) {
      double s = 0.0;
      for (int i = 0; i<n; ++i) {
         s += a[i];
      }

      return s;
   }

   /**
    * @brief Sum of absolute values with runtime dimension.
    * @param[in] a Input vector
    * @param[in] n Vector dimension (runtime parameter)
    * @return Sum of absolute values
    * 
    * Runtime variant of vecsssumabs<n> for dynamic vector sizes.
    * 
    * Mathematical operation: s = Σ_{i=0}^{n-1} |a_i|
    */
   __ecmech_hdev__
   inline double vecsssumabs_n(const double* const a, int n) {
      double s = 0.0;
      for (int i = 0; i<n; ++i) {
         s += fabs(a[i]);
      }

      return s;
   }

   /**
    * @brief Euclidean norm (L2 norm, magnitude): ‖v‖ = √(Σv_i²).
    * @tparam n Vector dimension
    * @param[in] v Input vector
    * @return Euclidean norm
    * 
    * Computes the L2 (Euclidean) norm of the vector.
    * 
    * Mathematical operation: ‖v‖ = √(Σ_{i=0}^{n-1} v_i²)
    * 
    * Common usage: Normalization, checking convergence, computing tensor magnitudes
    */
   template<int n>
   __ecmech_hdev__
   inline double vecNorm(const double* const v){
      double retval = 0.0;
      for (int i = 0; i<n; ++i) {
         retval += v[i] * v[i];
      }

      retval = sqrt(retval);
      return retval;
   }

   /**
    * @brief Scalar multiplication: v = x · a.
    * @tparam n Vector dimension
    * @param[out] v Result vector
    * @param[in] x Scalar multiplier
    * @param[in] a Input vector
    * 
    * Scales vector by scalar value.
    * 
    * Mathematical operation: v_i = x · a_i for i = 0, ..., n-1
    */
   template<int n>
   __ecmech_hdev__
   inline void vecsVxa(double* const v,
                       double x,
                       const double* const a) {
      for (int i = 0; i<n; ++i) {
         v[i] = x * a[i];
      }
   }

   /**
    * @brief Scale vector in place: a ← s · a.
    * @tparam n Vector dimension
    * @param[in,out] a Vector to scale (modified in place)
    * @param[in] s Scalar multiplier
    * 
    * Multiplies all vector components by scalar, modifying the original vector.
    * 
    * Mathematical operation: a_i ← s · a_i for i = 0, ..., n-1
    * 
    * @note Modifies input array in place (destructive operation)
    */
   template<int n>
   __ecmech_hdev__
   inline void vecsVsa(double* const a,
                       double s) {
      for (int i = 0; i<n; ++i) {
         a[i] *= s;
      }
   }

   /**
    * @brief Normalize vector in place: v ← v / ‖v‖.
    * @tparam n Vector dimension
    * @param[in,out] v Vector to normalize (modified in place)
    * 
    * Scales vector to unit length. If vector norm is below tolerance (idp_eps),
    * scales by large value (1/idp_eps) instead to avoid division by near-zero.
    * 
    * Mathematical operation:
    * - If ‖v‖ > ε: v ← v / ‖v‖
    * - If ‖v‖ ≤ ε: v ← v × (1/ε)  [effectively sets to large value]
    * 
    * Safeguard: Uses idp_eps threshold to prevent division by zero for
    * nearly-zero vectors.
    * 
    * @note Modifies input array in place
    * @note For ‖v‖ ≈ 0, does not produce unit vector but scales by 1/ε
    * 
    * Common usage: Normalizing quaternions, slip plane normals, direction vectors
    */
   template<int n>
   __ecmech_hdev__
   inline void vecsVNormalize(double* const v){
      const double norm = vecNorm<n>(v);
      const double s = (fabs(norm) > idp_eps) ? 1.0 / norm : 1.0 / idp_eps;
      vecsVsa<n>(v, s);
   }

   /**
    * @brief Vector cross product: v = a × b (3D only).
    * 
    * Computes the standard right-handed cross product of two 3D vectors.
    * 
    * **Mathematical operation**:
    * ```
    * v[0] = a[1]·b[2] - a[2]·b[1]
    * v[1] = a[2]·b[0] - a[0]·b[2]
    * v[2] = a[0]·b[1] - a[1]·b[0]
    * ```
    * 
    * @param[out] v Result vector (length 3), perpendicular to a and b
    * @param[in] a First input vector (length 3)
    * @param[in] b Second input vector (length 3)
    * 
    * @note ‖v‖ = ‖a‖ · ‖b‖ · sin(θ) where θ is angle between a and b
    * @note Direction follows right-hand rule
    */
   __ecmech_hdev__
   inline
   void
   vecCrossProd( // vec_cross_prod
      double* const val,
      const double* const v1,
      const double* const v2) {
      val[0] = v1[1] * v2[2] - v1[2] * v2[1];
      val[1] = -v1[0] * v2[2] + v1[2] * v2[0];
      val[2] = v1[0] * v2[1] - v1[1] * v2[0];
   }

   /** @} */ // end of vector_operations group

   /**
    * @defgroup matrix_operations Matrix and Matrix-Vector Operations
    * @brief Linear algebra operations for square and non-square matrices.
    * 
    * All matrices stored in row-major order using ECMECH_NN_INDX and ECMECH_NM_INDX.
    * 
    * Naming conventions:
    * - vecsV*: Output is vector
    * - vecsM*: Output is matrix
    * - *Ma*: Matrix times vector (a)
    * - *MTa*: Matrix transpose times vector (a)
    * - *aTM*: Vector transpose times matrix (a^T · M)
    * - *ABT*: Matrix times matrix transpose (A · B^T)
    * - *AB*: Matrix times matrix (A · B)
    * - *aTb*: Outer product (a ⊗ b)
    * 
    * @{
    */

   /**
    * @brief Matrix-vector product: p = M^T · v for square matrix.
    * 
    * Computes the product of a transposed square matrix with a vector.
    * This is equivalent to computing the linear combination of matrix columns
    * weighted by vector components.
    * 
    * **Mathematical operation**: p_i = Σⱼ M[j,i] · v[j] for i = 0, ..., n-1
    * 
    * @tparam n Matrix dimension (n×n) and vector length
    * @param[out] p Result vector (length n)
    * @param[in] M Square matrix in row-major order (n×n)
    * @param[in] v Input vector (length n)
    * 
    * @note Unlike FORTRAN matt_x_vec_5, output is first argument
    */
   template<int n>
   __ecmech_hdev__
   inline void vecsVMTa(double* const p,
                        const double* const M,
                        const double* const v) {
      for (int i = 0; i<n; ++i) {
         p[i] = 0.0;
         for (int j = 0; j<n; ++j) {
            p[i] += M[ECMECH_NN_INDX(j, i, n)] * v[j];
         }
      }
   }

   /**
    * @brief Vector-matrix product: p = a^T · M for non-square matrix (n×q).
    * 
    * Computes the product of a row vector with a non-square matrix,
    * resulting in a row vector of dimension q.
    * 
    * **Mathematical operation**: p[j] = Σᵢ a[i] · M[i,j] for j = 0, ..., q-1
    * 
    * @tparam n Number of rows in M, length of vector a
    * @tparam q Number of columns in M, length of result p
    * @param[out] p Result vector (length q)
    * @param[in] a Input row vector (length n)
    * @param[in] M Matrix in row-major order (n×q)
    */
   template<int n, int q>
   __ecmech_hdev__
   inline void vecsVaTM(double* const p,
                        const double* const a,
                        const double* const M) {
      for (int iQ = 0; iQ<q; ++iQ) {
         p[iQ] = 0.0;
         for (int iN = 0; iN<n; ++iN) {
            p[iQ] += a[iN] * M[ECMECH_NM_INDX(iN, iQ, n, q)];
         }
      }
   }

   /**
    * @brief Matrix-vector product: p = M · a for non-square matrix (n×q).
    * 
    * Computes the standard matrix-vector multiplication where M is non-square.
    * 
    * **Mathematical operation**: p[i] = Σⱼ M[i,j] · a[j] for i = 0, ..., n-1
    * 
    * @tparam n Number of rows in M, length of result p
    * @tparam q Number of columns in M, length of vector a
    * @param[out] p Result vector (length n)
    * @param[in] M Matrix in row-major order (n×q)
    * @param[in] a Input vector (length q)
    */
   template<int n, int q>
   __ecmech_hdev__
   inline void vecsVMa(double* const p,
                       const double* const M,
                       const double* const a) {
      for (int iN = 0; iN<n; ++iN) {
         p[iN] = 0.0;
         for (int iQ = 0; iQ<q; ++iQ) {
            p[iN] += M[ECMECH_NM_INDX(iN, iQ, n, q)] * a[iQ];
         }
      }
   }

   /**
    * @brief Matrix-vector product: p = M · a for square matrix (n×n).
    * 
    * Computes the standard matrix-vector multiplication where M is square.
    * 
    * **Mathematical operation**: p[i] = Σⱼ M[i,j] · a[j] for i = 0, ..., n-1
    * 
    * @tparam n Matrix dimension (n×n) and vector length
    * @param[out] p Result vector (length n)
    * @param[in] M Square matrix in row-major order (n×n)
    * @param[in] a Input vector (length n)
    */
   template<int n>
   __ecmech_hdev__
   inline void vecsVMa(double* const p,
                       const double* const M,
                       const double* const a) {
      for (int iN = 0; iN<n; ++iN) {
         p[iN] = 0.0;
         for (int jN = 0; jN<n; ++jN) {
            p[iN] += M[ECMECH_NN_INDX(iN, jN, n)] * a[jN];
         }
      }
   }

   /**
    * @brief Matrix product P = A · B^T where A and B are (n×q), result is (n×n).
    * 
    * Computes the product of matrix A with the transpose of matrix B.
    * Both A and B have the same dimensions (n×q), producing a square result.
    * Commonly used for computing Gram matrices or outer product sums.
    * 
    * **Mathematical operation**: P[i,j] = Σₖ A[i,k] · B[j,k] for i,j = 0, ..., n-1
    * 
    * @tparam n Number of rows in A and B, dimension of square result P
    * @tparam q Number of columns in A and B (inner product dimension)
    * @param[out] P Result matrix (n×n) in row-major order
    * @param[in] A First matrix (n×q) in row-major order
    * @param[in] B Second matrix (n×q) in row-major order
    * 
    * @note Output P is initialized to zero before accumulation
    */
   template<int n, int q>
   __ecmech_hdev__
   inline void vecsMABT(double* const P,
                        const double* const A,
                        const double* const B)
   {
      for (int ij = 0; ij<n * n; ++ij) {
         P[ij] = 0.0;
      }

      for (int iN = 0; iN < n; ++iN) {
         for (int jN = 0; jN < n; ++jN) {
            for (int iQ = 0; iQ < q; ++iQ) {
               P[ECMECH_NM_INDX(iN, jN, n, n)] += A[ECMECH_NM_INDX(iN, iQ, n, q)] * B[ECMECH_NM_INDX(jN, iQ, n, q)];
            }
         }
      }
   }

   /**
    * @brief Matrix product P = A · B^T where A is (n×q), B is (m×q), result is (n×m).
    * 
    * Computes the product of matrix A with the transpose of matrix B.
    * A and B can have different numbers of rows but must have the same
    * number of columns (q).
    * 
    * **Mathematical operation**: P[i,j] = Σₖ A[i,k] · B[j,k] for i=0,...,n-1, j=0,...,m-1
    * 
    * @tparam n Number of rows in A, first dimension of result P
    * @tparam m Number of rows in B, second dimension of result P
    * @tparam q Number of columns in A and B (inner product dimension)
    * @param[out] P Result matrix (n×m) in row-major order
    * @param[in] A First matrix (n×q) in row-major order
    * @param[in] B Second matrix (m×q) in row-major order
    * 
    * @note Output P is initialized to zero before accumulation
    */
   template<int n, int m, int q>
   __ecmech_hdev__
   inline void vecsMABT(double* const P,
                        const double* const A,
                        const double* const B) {
      for (int ij = 0; ij<n * m; ++ij) {
         P[ij] = 0.0;
      }

      for (int iN = 0; iN < n; ++iN) {
         for (int jM = 0; jM < m; ++jM) {
            for (int iQ = 0; iQ < q; ++iQ) {
               P[ECMECH_NM_INDX(iN, jM, n, m)] += A[ECMECH_NM_INDX(iN, iQ, n, q)] * B[ECMECH_NM_INDX(jM, iQ, m, q)];
            }
         }
      }
   }

   /**
    * @brief Matrix product P = A · B where A is (n×q), B is (q×m), result is (n×m).
    * 
    * Computes the standard matrix-matrix product. This is the most general
    * rectangular matrix multiplication.
    * 
    * **Mathematical operation**: P[i,j] = Σₖ A[i,k] · B[k,j] for i=0,...,n-1, j=0,...,m-1
    * 
    * @tparam n Number of rows in A, first dimension of result P
    * @tparam m Number of columns in B, second dimension of result P
    * @tparam q Number of columns in A / rows in B (inner dimension)
    * @param[out] P Result matrix (n×m) in row-major order
    * @param[in] A First matrix (n×q) in row-major order
    * @param[in] B Second matrix (q×m) in row-major order
    * 
    * @note Output P is initialized to zero before accumulation
    */
   template<int n, int m, int q>
   __ecmech_hdev__
   inline void vecsMAB(double* const P,
                       const double* const A,
                       const double* const B) {
      for (int ij = 0; ij<n * m; ++ij) {
         P[ij] = 0.0;
      }

      for (int iN = 0; iN < n; ++iN) {
         for (int iQ = 0; iQ < q; ++iQ) {
            for (int jM = 0; jM < m; ++jM) {
               P[ECMECH_NM_INDX(iN, jM, n, m)] += A[ECMECH_NM_INDX(iN, iQ, n, q)] * B[ECMECH_NM_INDX(iQ, jM, q, m)];
            }
         }
      }
   }

   /**
    * @brief Outer product P = a ⊗ b where a and b are n-vectors, result is (n×n).
    * 
    * Computes the tensor (outer) product of two vectors, producing a matrix.
    * Each element P[i,j] is simply the product of a[i] and b[j].
    * 
    * **Mathematical operation**: P[i,j] = a[i] · b[j] for i,j = 0, ..., n-1
    * 
    * @tparam n Vector length and matrix dimension
    * @param[out] P Result matrix (n×n) in row-major order
    * @param[in] a First input vector (length n)
    * @param[in] b Second input vector (length n)
    * 
    * @note For symmetric outer product (a ⊗ a), use same vector for both arguments
    */
   template<int n>
   __ecmech_hdev__
   inline void vecsMaTb(double* const P,
                        const double* const a,
                        const double* const b)
   {
      for (int iN = 0; iN < n; ++iN) {
         for (int jN = 0; jN < n; ++jN) {
            P[ECMECH_NN_INDX(iN, jN, n)] = a[iN] * b[jN];
         }
      }
   }

   /**
    * @brief Compute symmetric part of matrix: P = (A + A^T) / 2.
    * 
    * Extracts the symmetric part of a matrix by averaging it with its transpose.
    * This is the orthogonal projection onto the space of symmetric matrices.
    * 
    * **Mathematical operation**: P[i,j] = (A[i,j] + A[j,i]) / 2 for i,j = 0, ..., n-1
    * 
    * @tparam n Matrix dimension (n×n)
    * @param[out] P Symmetric part of A (n×n) in row-major order
    * @param[in] A Input matrix (n×n) in row-major order
    * 
    * @note P is guaranteed symmetric: P[i,j] = P[j,i]
    * @note Can be used in-place: P and A can be the same array
    */
   template<int n>
   __ecmech_hdev__
   inline void vecsMsymm(double* const P,
                         const double* const A)
   {
      for (int iN = 0; iN < n; ++iN) {
         for (int jN = 0; jN < n; ++jN) {
            P[ECMECH_NN_INDX(iN, jN, n)] = 0.5 * (A[ECMECH_NN_INDX(iN, jN, n)] + A[ECMECH_NN_INDX(jN, iN, n)]);
         }
      }
   }

   /**
    * @brief Compute skew-symmetric part of matrix: Q = (A - A^T) / 2.
    * 
    * Extracts the skew-symmetric (antisymmetric) part of a matrix by computing
    * the difference with its transpose. This is the orthogonal projection onto
    * the space of skew-symmetric matrices.
    * 
    * **Mathematical operation**: Q[i,j] = (A[i,j] - A[j,i]) / 2 for i,j = 0, ..., n-1
    * 
    * @tparam n Matrix dimension (n×n)
    * @param[out] Q Skew-symmetric part of A (n×n) in row-major order
    * @param[in] A Input matrix (n×n) in row-major order
    * 
    * @note Q is guaranteed skew-symmetric: Q[i,j] = -Q[j,i], Q[i,i] = 0
    * @note Can be used in-place: Q and A can be the same array
    */
   template<int n>
   __ecmech_hdev__
   inline void vecsMskew(double* const Q,
                         const double* const A)
   {
      for (int iN = 0; iN < n; ++iN) {
         for (int jN = 0; jN < n; ++jN) {
            Q[ECMECH_NN_INDX(iN, jN, n)] = 0.5 * (A[ECMECH_NN_INDX(iN, jN, n)] - A[ECMECH_NN_INDX(jN, iN, n)]);
         }
      }
   }

   /** @} */ // end of matrix_operations group


   /**
    * @defgroup tensor_conversions Tensor Representation Conversions
    * @brief Convert between different tensor storage formats.
    * 
    * ECMech uses multiple representations for symmetric and skew-symmetric tensors:
    * - Full 3×3 matrix (9 components, row-major)
    * - Deviatoric 5-vector (traceless symmetric, normalized)
    * - Voigt 6-vector (general symmetric: 11, 22, 33, 23, 31, 12)
    * - Extended 7-vector (deviatoric + pressure: vecd + p)
    * - Spin 3-vector (skew-symmetric, axial vector)
    * 
    * Normalization preserves inner products between representations.
    * 
    * @{
    */

   /**
    * @brief Compute trace of 3×3 matrix: tr(A) = A[0,0] + A[1,1] + A[2,2].
    * 
    * **Mathematical operation**: tr(A) = Σᵢ A[i,i] = A₀₀ + A₁₁ + A₂₂
    * 
    * @param[in] A Matrix in row-major order (3×3 = 9 components)
    * @return Trace of the matrix (sum of diagonal elements)
    */
   __ecmech_hdev__
   inline double trace3(const double* const A // (DIMS,DIMS)
                        ) {
      double trace =
         A[ECMECH_NN_INDX(0, 0, ecmech::ndim)] +
         A[ECMECH_NN_INDX(1, 1, ecmech::ndim)] +
         A[ECMECH_NN_INDX(2, 2, ecmech::ndim)];

      return trace;
   }

   /**
    * @brief Convert trace to deviatoric-volumetric (vecds) scalar component.
    * 
    * Computes the normalized spherical (volumetric) component from a trace value.
    * Used in deviatoric-spherical decomposition of tensors.
    * 
    * **Mathematical operation**: vecds_s = dkk / √3
    * 
    * @param[in] dkk Trace value (sum of diagonal components)
    * @return Spherical component in vecds representation
    * 
    * @note The √3 normalization preserves inner product norms
    */
   __ecmech_hdev__
   inline double traceToVecdsS(double dkk) {
      double vecds_s = sqr3i * dkk;
      return vecds_s;
   }

   /**
    * @brief Compute effective strain rate from deviatoric components.
    * 
    * Computes the von Mises equivalent (effective) strain rate from the
    * 5-component deviatoric representation. This is the scalar measure of
    * deformation rate intensity.
    * 
    * **Mathematical operation**: D_eff = √(2/3 · D_dev : D_dev)
    * 
    * @param[in] vecd Deviatoric strain rate (5 components)
    * @return Effective (von Mises equivalent) strain rate
    * 
    * @note Always non-negative by definition
    * @note Equivalent to sqrt(2/3) times the deviatoric tensor norm
    */
   __ecmech_hdev__
   inline double vecd_Deff(const double* const vecd) {
      double retval = vecNorm<ecmech::ntvec>(vecd);
      retval = sqr2b3 * retval;
      return retval;
   }

   /**
    * @brief Inner product of deviatoric parts of two extended 7-vectors.
    * 
    * Computes the Frobenius inner product of the deviatoric parts only,
    * ignoring the pressure (component 6) in both vectors.
    * 
    * **Mathematical operation**:
    * ```
    * result = Σ_{i=0}^{4} stress[i] · d[i]
    * ```
    * 
    * @param[in] stressSvec First extended 7-vector (typically stress)
    * @param[in] dSvec Second extended 7-vector (typically strain rate)
    * @return Inner product of deviatoric parts
    * 
    * @note Only components 0-4 contribute (deviatoric)
    * @note Component 5 is skipped, component 6 (pressure) is ignored
    * @note Common usage: Computing deviatoric stress power
    */
   __ecmech_hdev__
   inline double vecsInnerSvecDev(const double* const stressSvec,
                                  const double* const dSvec) {
      double retval =
         stressSvec[0] * dSvec[0] +
         stressSvec[1] * dSvec[1] +
         stressSvec[2] * dSvec[2] +
         two * (
            stressSvec[3] * dSvec[3] +
            stressSvec[4] * dSvec[4] +
            stressSvec[5] * dSvec[5]);
      return retval;
   }

   /**
    * @brief Convert skew-symmetric 3×3 matrix to axial 3-vector (spin vector).
    * 
    * Represents a skew-symmetric tensor by its dual axial vector.
    * For a skew-symmetric matrix W, the axial vector w satisfies:
    * W·v = w × v for any vector v.
    * 
    * **Component mapping**:
    * ```
    * veccp[0] = W[2,1] = -W[1,2]  (rotation about x₁-axis)
    * veccp[1] = W[0,2] = -W[2,0]  (rotation about x₂-axis)
    * veccp[2] = W[1,0] = -W[0,1]  (rotation about x₃-axis)
    * ```
    * 
    * @param[out] veccp Axial 3-vector (spin components)
    * @param[in] W Skew-symmetric 3×3 matrix in row-major order
    * 
    * @note Only 3 independent components in skew-symmetric 3×3 matrix
    * @note Magnitude ‖w‖ equals the rotation rate
    */
   __ecmech_hdev__
   inline void skewToVeccp(double* const veccp, // (WVEC)
                           const double* const W // (DIMS,DIMS)
                           )
   {
      veccp[0] = W[ECMECH_NN_INDX(2, 1, ecmech::ndim)];
      veccp[1] = W[ECMECH_NN_INDX(0, 2, ecmech::ndim)];
      veccp[2] = W[ECMECH_NN_INDX(1, 0, ecmech::ndim)];
   }

   /**
    * @brief Convert symmetric 3×3 tensor to deviatoric 5-vector.
    * 
    * Extracts the traceless (deviatoric) part of a symmetric tensor and
    * represents it in normalized 5-component form.
    * 
    * **Component mapping**:
    * ```
    * vecd[0] = (T₁₁ - T₂₂) / √2
    * vecd[1] = (2·T₃₃ - T₁₁ - T₂₂) / √6
    * vecd[2] = √2 · T₂₃
    * vecd[3] = √2 · T₁₃
    * vecd[4] = √2 · T₁₂
    * ```
    * 
    * @param[out] vecd Deviatoric 5-vector (traceless symmetric representation)
    * @param[in] A Symmetric 3×3 tensor in row-major order (9 components)
    * 
    * @note The √2 and √6 factors preserve the Frobenius norm
    * @note Result is guaranteed traceless: vecd represents T' = T - (tr(T)/3)·I
    */
   __ecmech_hdev__
   inline void symmToVecd(double* const vecd, // (TVEC)
                          const double* const A // (DIMS,DIMS)
                          )
   {
      vecd[0] = sqr2i * (A[ECMECH_NN_INDX(0, 0, ecmech::ndim)] - A[ECMECH_NN_INDX(1, 1, ecmech::ndim)]);
      vecd[1] = sqr6i *
                (two *
                 A[ECMECH_NN_INDX(2, 2,
                                  ecmech::ndim)] - A[ECMECH_NN_INDX(0, 0, ecmech::ndim)] - A[ECMECH_NN_INDX(1, 1, ecmech::ndim)]); // = sqr6i * (3*A33 - Akk) = sqr3b2 * (A33 - Akk/3) = sqr3b2 * Adev33
      vecd[2] = sqr2 * A[ECMECH_NN_INDX(1, 0, ecmech::ndim)];
      vecd[3] = sqr2 * A[ECMECH_NN_INDX(2, 0, ecmech::ndim)];
      vecd[4] = sqr2 * A[ECMECH_NN_INDX(2, 1, ecmech::ndim)];
   }

   /**
    * @brief Convert symmetric 3×3 tensor to deviatoric-spherical 6-vector.
    * 
    * Decomposes a symmetric tensor into deviatoric (5 components) plus
    * spherical/volumetric (1 component) parts.
    * 
    * **Component mapping**:
    * ```
    * vecds[0:4] = deviatoric part (as in symmToVecd)
    * vecds[5] = tr(T) / √3  (spherical part)
    * ```
    * 
    * @param[out] vecds Deviatoric-spherical 6-vector
    * @param[in] A Symmetric 3×3 tensor in row-major order (9 components)
    * 
    * @note This is a complete decomposition: T = T' + (tr(T)/3)·I
    * @note Preserves norm: ‖vecds‖² = ‖T'‖² + (tr(T))²/3
    */
   __ecmech_hdev__
   inline void symmToVecds(double* const vecds, // (SVEC)
                           const double* const A // (DIMS,DIMS)
                           )
   {
      symmToVecd(vecds, A);
      double Akk = trace3(A);
      vecds[iSvecS] = traceToVecdsS(Akk);
   }


   /**
    * @brief Convert Voigt 6-vector (symmetric) to deviatoric 5-vector.
    * 
    * Extracts deviatoric part from Voigt representation and converts to
    * normalized 5-component form.
    * 
    * **Input format** (Voigt): svec_kk = [T₁₁, T₂₂, T₃₃, T₂₃, T₁₃, T₁₂]
    * 
    * **Output format**: vecd[0:4] = deviatoric part with √2, √6 normalization
    * 
    * @param[out] vecd Deviatoric 5-vector
    * @param[in] svec_kk Symmetric tensor in Voigt 6-vector form
    * 
    * @note Voigt indexing: [11, 22, 33, 23, 13, 12] (engineering notation)
    */
   __ecmech_hdev__
   inline void svecToVecd(double* const vecd, // (TVEC)
                          const double* const svec_kk // (SVEC[+1])
                          )
   {
      vecd[0] = sqr2i * (svec_kk[0] - svec_kk[1]);
      vecd[1] = sqr3b2 * svec_kk[2];
      vecd[2] = sqr2 * svec_kk[5];
      vecd[3] = sqr2 * svec_kk[4];
      vecd[4] = sqr2 * svec_kk[3];
      // vecds[5] = traceToVecdsS( svec_kk[6] ) ;
   }

   /**
    * @brief Convert deviatoric-spherical 6-vector to extended 7-vector with pressure.
    * 
    * Transforms vecds (deviatoric + spherical) to svecp format used in stress updates.
    * 
    * **Component mapping**:
    * ```
    * svecp[0:4] = vecds[0:4]  (deviatoric part, unchanged)
    * svecp[5] = 0.0           (unused/placeholder)
    * svecp[6] = -vecds[5]     (negative of spherical part = pressure)
    * ```
    * 
    * @param[out] svecp Extended 7-vector (deviatoric stress + pressure)
    * @param[in] vecds Deviatoric-spherical 6-vector
    * 
    * @note svecp[6] = -p where p is the mean normal stress (pressure)
    * @note svecp[5] is typically unused but may be repurposed in some contexts
    */
   __ecmech_hdev__
   inline void vecdsToSvecP(double* const svecp, // (SVEC+1)
                            const double* const vecds // (SVEC)
                            )
   {
      svecp[iSvecP] = -sqr3i * vecds[iSvecS]; // -Akk_by_3

      double t1 = sqr2i * vecds[0];
      double t2 = sqr6i * vecds[1];
      //
      svecp[0] = t1 - t2; // 11'
      svecp[1] = -t1 - t2; // 22'
      svecp[2] = sqr2b3 * vecds[1]; // 33'
      svecp[3] = sqr2i * vecds[4]; // 23
      svecp[4] = sqr2i * vecds[3]; // 31
      svecp[5] = sqr2i * vecds[2]; // 12
   }

   /**
    * @brief Convert extended 7-vector to standard Voigt 6-vector.
    * 
    * Reconstructs full symmetric tensor in Voigt form from deviatoric + pressure.
    * 
    * **Component reconstruction**:
    * ```
    * T₁₁ = dev₁₁ - p  (svecp[6] = -p)
    * T₂₂ = dev₂₂ - p
    * T₃₃ = dev₃₃ - p
    * T₂₃, T₁₃, T₁₂ = off-diagonals from svecp
    * ```
    * 
    * @param[out] a_svec Symmetric tensor in Voigt 6-vector form [11, 22, 33, 23, 13, 12]
    * @param[in] a_svec_p Extended 7-vector (deviatoric + pressure)
    * 
    * @note Adds back the spherical part to get the full tensor
    */
   __ecmech_hdev__
   inline void svecpToSvec(double* const a_svec,
                           const double* const a_svec_p
                           )
   {
      double a_mean = -a_svec_p[iSvecP];

      for (int i_svec = 0; i_svec < ecmech::nsvec; i_svec++) {
         a_svec[i_svec] = a_svec_p[i_svec];
      }

      a_svec[0] = a_svec[0] + a_mean;
      a_svec[1] = a_svec[1] + a_mean;
      a_svec[2] = a_svec[2] + a_mean;
   }

   /**
    * @brief Decompose 3×3 matrix into deviatoric symmetric and skew parts.
    * 
    * Decomposes an arbitrary 3×3 tensor into:
    * - P_vecd: Deviatoric symmetric part (5 components, traceless)
    * - Q_veccp: Skew-symmetric part (3 components, axial vector form)
    * 
    * **Mathematical operation**:
    * ```
    * Symmetric part: S = (T + T^T) / 2
    * Deviatoric: P = S - (tr(S)/3)·I
    * Skew part: Q = (T - T^T) / 2
    * ```
    * 
    * @param[out] P_vecd Deviatoric symmetric part (5 components)
    * @param[out] Q_veccp Skew-symmetric part in axial form (3 components)
    * @param[in] T Full 3×3 tensor in row-major order
    * 
    * @note Recovers T via: T = P + (tr(T)/3)·I + Q
    * @note Commonly used in slip system kinematics: T_ref = P + Q
    */
   __ecmech_hdev__
   inline
   void matToPQ(double* const P_vecd, // ntvec
                double* const Q_veccp, // nwvec
                const double* const T // ndim*ndim
                )
   {
      // CALL mat_to_symm_3(crys%p_ref(:,:,is), crys%t_ref(:,:,is))
      // CALL symm_to_vecds(P_ref_svec, crys%p_ref(:,:,is))
      // crys%P_ref_vec(:, is) = P_ref_svec(1:TVEC)
      //
      double P[ ecmech::ndim * ecmech::ndim ];
      vecsMsymm<ndim>(P, T);
      symmToVecd(P_vecd, P);

      // CALL mat_to_skew_3(crys%q_ref(:,:,is), crys%t_ref(:,:,is))
      // CALL skew_to_veccp(crys%q_ref_vec(:,is), crys%q_ref(:,:,is))
      double Q[ ecmech::ndim * ecmech::ndim ];
      vecsMskew<ndim>(Q, T);
      skewToVeccp(Q_veccp, Q);
   }

   /** @} */ // end of tensor_conversions group

   /**
    * @defgroup quaternion_operations Quaternion Operations
    * @brief Unit quaternion algebra for rotation representation.
    * 
    * Quaternions q = [q0, q1, q2, q3] where q0 is scalar, [q1,q2,q3] is vector.
    * For rotations: ‖q‖ = 1 (unit quaternions).
    * 
    * **Conventions**:
    * - Ordering: [q0, q1, q2, q3] = [scalar, vector_x, vector_y, vector_z]
    * - Unit constraint: q0² + q1² + q2² + q3² = 1
    * - Composition: q = q1 ⊗ q2 means "rotate by q2, then by q1"
    * - Inverse: q⁻¹ = [q0, -q1, -q2, -q3] for unit quaternions
    * 
    * **Relationship to rotations**:
    * - Quaternion q represents rotation by angle θ about axis n̂:
    *   q = [cos(θ/2), sin(θ/2)·n̂]
    * - Exponential map: e = θ·n̂ (angle-axis, 3 components)
    * 
    * @{
    */


   /**
    * @brief Convert inverse exponential map (angle-axis in 4-component form) to quaternion.
    * 
    * Takes an "inverse" representation which is a 4-component angle-axis vector:
    * inv = [θ, n̂_x, n̂_y, n̂_z] where θ is the rotation angle and n̂ is the unit axis.
    * 
    * This is an intermediate representation used internally by emap_to_quat:
    * - emap_to_quat computes θ = ‖emap‖ and n̂ = emap/θ
    * - Then calls inv_to_quat with inv = [θ, n̂_x, n̂_y, n̂_z]
    * 
    * **Mathematical operation**:
    * ```
    * θ = inv[0]  (total rotation angle in radians)
    * n̂ = [inv[1], inv[2], inv[3]]  (unit rotation axis)
    * quat = [cos(θ/2), sin(θ/2)·n̂]
    * ```
    * 
    * **Implementation**:
    * ```
    * a = θ / 2
    * quat[0] = cos(a)
    * quat[1:3] = sin(a) · n̂
    * ```
    * 
    * @param[out] quat Quaternion [q0, q1, q2, q3] (length 4, unit norm)
    * @param[in] inv Inverse representation [θ, n̂_x, n̂_y, n̂_z] (length 4)
    *            - inv[0]: Rotation angle θ in radians
    *            - inv[1:3]: Unit rotation axis n̂ (should be normalized)
    * 
    * @note Name "inv" is historical/legacy - not related to quaternion inverse
    * @note This is a helper function, typically called via emap_to_quat
    * @note Input axis n̂ should be unit length for correct quaternion norm
    */
   __ecmech_hdev__
   inline void inv_to_quat(double* const quat,
                           const double* const inv) {
      double a = inv[0] * 0.5;
      quat[0] = cos(a);
      a = sin(a);
      vecsVxa<nwvec>(&(quat[1]), a, &(inv[1]));
   }

   /**
    * @brief Convert exponential map (angle-axis) to quaternion.
    * 
    * Converts a rotation represented as angle-axis (exponential map) to
    * unit quaternion form.
    * 
    * **Mathematical operation**:
    * ```
    * θ = ‖emap‖  (rotation angle in radians)
    * n̂ = emap / θ  (rotation axis, unit vector)
    * quat = [cos(θ/2), sin(θ/2)·n̂]
    * ```
    * 
    * @param[out] quat Quaternion [q0, q1, q2, q3] (length 4, unit norm)
    * @param[in] emap Exponential map θ·n̂ (length 3, angle-axis representation)
    * 
    * @note Special case: ‖emap‖ = 0 → quat = [1, 0, 0, 0] (identity rotation)
    * @note ‖quat‖ = 1 (unit quaternion) upon return
    */
   __ecmech_hdev__
   inline void emap_to_quat(double* const quat,
                            const double* const emap) {
      double inv[invdim] = { 0.0, 1.0, 0.0, 0.0 };
      inv[0] = vecNorm<emapdim>(emap);
      if (inv[0] > idp_tiny_sqrt) {
         double invInv = 1.0 / inv[0];
         vecsVxa<emapdim>(&(inv[1]), invInv, emap);
      } // else, emap is effectively zero, so axis does not matter
      inv_to_quat(quat, inv);
   }

   /**
    * @brief Convert quaternion to exponential map (angle-axis).
    * 
    * Converts a unit quaternion to angle-axis (exponential map) representation.
    * 
    * **Mathematical operation**:
    * ```
    * θ = 2·acos(q0)  (rotation angle, 0 ≤ θ ≤ π)
    * n̂ = [q1, q2, q3] / sin(θ/2)  (rotation axis)
    * emap = θ·n̂
    * ```
    * 
    * @param[out] emap Exponential map θ·n̂ (length 3)
    * @param[in] quat Quaternion [q0, q1, q2, q3] (length 4, should be unit norm)
    * 
    * @note Special case: q0 ≈ ±1 → emap ≈ [0, 0, 0] (near identity)
    * @note Sign convention: Always produces θ ∈ [0, π]
    */
   __ecmech_hdev__
   inline void quat_to_emap(double* const emap,
                            const double* const quat) {

      constexpr auto tol = std::numeric_limits<double>::epsilon();
      const auto phi = 2.0 * acos(quat[0]);

      if (fabs(quat[0]) < tol) {
         emap[0] = quat[1] * M_PI;
         emap[1] = quat[2] * M_PI;
         emap[2] = quat[3] * M_PI;
      } else {
         const double sign = (quat[0] < 0.0) ? -1.0 : 1.0; 
         const double s = sign / sqrt(quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]);
         emap[0] = s * quat[1] * phi;
         emap[1] = s * quat[2] * phi;
         emap[2] = s * quat[3] * phi; 
      }
   }

   /**
    * @brief Quaternion multiplication: q = q1 ⊗ q2 (Hamilton product).
    * 
    * Computes the Hamilton product of two quaternions. For unit quaternions
    * representing rotations, this composition corresponds to sequential rotations:
    * first q2, then q1.
    * 
    * **Quaternion convention**: q = [q0, q1, q2, q3] where q0 is scalar part
    * 
    * **Mathematical operation**:
    * ```
    * q0 = q1_0·q2_0 - q1_1·q2_1 - q1_2·q2_2 - q1_3·q2_3
    * q1 = q1_0·q2_1 + q1_1·q2_0 + q1_2·q2_3 - q1_3·q2_2
    * q2 = q1_0·q2_2 - q1_1·q2_3 + q1_2·q2_0 + q1_3·q2_1
    * q3 = q1_0·q2_3 + q1_1·q2_2 - q1_2·q2_1 + q1_3·q2_0
    * ```
    * 
    * @param[out] q Result quaternion (length 4)
    * @param[in] q1 First quaternion (length 4)
    * @param[in] q2 Second quaternion (length 4)
    * 
    * @note For rotations: q = q1 ⊗ q2 means "rotate by q2, then by q1"
    * @note Result is NOT automatically normalized
    * @note Quaternion multiplication is NON-commutative: q1 ⊗ q2 ≠ q2 ⊗ q1
    */
   __ecmech_hdev__
   inline void quat_prod(double* const q,
                         const double* const a,
                         const double* const b) {
      q[0] = a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3];

      q[1] = a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2];
      q[2] = a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1];
      q[3] = a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0];
   }

   /**
    * @brief Compute quaternion inverse: q_inv = [q0, -q1, -q2, -q3].
    * 
    * For unit quaternions (‖q‖ = 1), the inverse is simply the conjugate.
    * Represents the opposite rotation: if q rotates by θ about n̂, then
    * q⁻¹ rotates by -θ about n̂ (or equivalently θ about -n̂).
    * 
    * **Mathematical operation**:
    * ```
    * q_inv = [q0, -q1, -q2, -q3]  (for unit quaternions)
    * ```
    * 
    * @param[out] inv_quat Inverse quaternion (length 4)
    * @param[in] quat Input quaternion (length 4, should be unit norm)
    * 
    * @note For non-unit quaternions: q⁻¹ = conj(q) / ‖q‖²
    * @note Property: q ⊗ q⁻¹ = q⁻¹ ⊗ q = [1, 0, 0, 0] (identity)
    */
   __ecmech_hdev__
   inline void quat_inverse(double* const inv_quat,
                            const double* const quat) {
      // I mean this should be equal to 1 as we're dealing with unit quats...
      const double inv_quat_norm = 1.0 / vecNorm<ecmech::qdim>(quat);
      inv_quat[0] = inv_quat_norm * quat[0];
      inv_quat[1] = -inv_quat_norm * quat[1];
      inv_quat[2] = -inv_quat_norm * quat[2];
      inv_quat[3] = -inv_quat_norm * quat[3];
   }

   /**
    * @brief Compute relative rotation quaternion: qprime = q_new ⊗ q_old⁻¹.
    * 
    * Computes the incremental rotation from orientation q_old to q_new.
    * The result qprime represents the rotation that, when applied to q_old,
    * produces q_new.
    * 
    * **Mathematical operation**:
    * ```
    * qprime = q_new ⊗ q_old⁻¹
    * ```
    * where q_old⁻¹ = [q_old[0], -q_old[1], -q_old[2], -q_old[3]]
    * 
    * @param[out] qprime Relative rotation quaternion (length 4)
    * @param[in] quat_new Final orientation (length 4)
    * @param[in] quat_old Initial orientation (length 4)
    * 
    * @note Applying qprime to q_old: q_new = qprime ⊗ q_old
    * @note Useful for computing rotation increments in time integration
    */
   __ecmech_hdev__
   inline void quat_rel_rotation(double* const qprime,
                                 const double* const q1,
                                 const double* const q2) {
      double q1_inv[4] = {1.0, 0.0, 0.0, 0.0};
      quat_inverse(q1_inv, q1);
      quat_prod(qprime, q1_inv, q2);
   }

   /**
    * @brief Update crystal-to-sample quaternion using rotation increment.
    * 
    * Computes the current time step crystal-to-sample orientation quaternion
    * by applying a rotation increment to the previous time step orientation.
    * 
    * **Physical interpretation**:
    * - cn_quat: Crystal-to-sample rotation at previous time step (t_n)
    * - dr_quat: Incremental rotation over time step Δt (rotation increment)
    * - c_quat: Crystal-to-sample rotation at current time step (t_{n+1})
    * 
    * **Mathematical operation**:
    * ```
    * c_quat = cn_quat ⊗ dr_quat
    * ```
    * where ⊗ is quaternion multiplication (Hamilton product)
    * 
    * **Time integration context**:
    * In crystal plasticity, the lattice orientation evolves due to plastic spin:
    * ```
    * dR/dt = W^p · R
    * ```
    * 
    * The discrete update is:
    * ```
    * R_{n+1} = ΔR · R_n
    * ```
    * where ΔR is the rotation over time step Δt.
    * 
    * In quaternion form:
    * ```
    * q_{n+1} = Δq ⊗ q_n
    * ```
    * 
    * **Usage in implicit solvers**:
    * - dr_quat is computed from the exponential map: dr_quat = exp(Δω)
    * - Δω is the rotation vector increment (often an implicit solve variable)
    * - This function updates the orientation as part of the Newton-Raphson iteration
    * 
    * @param[out] c_quat Crystal-to-sample quaternion at current step t_{n+1} (length 4)
    * @param[in] dr_quat Incremental rotation quaternion over Δt (length 4, unit norm)
    *                    - Represents rotation increment: Δq = exp(Δω)
    *                    - Typically computed via emap_to_quat(dr_quat, delta_omega)
    * @param[in] cn_quat Crystal-to-sample quaternion at previous step t_n (length 4, unit norm)
    *                    - Previous orientation (from beginning of time step)
    * 
    * @note Quaternion multiplication order: c = cn ⊗ dr (increment applied "from the right")
    * @note Result is NOT automatically normalized (may need normalization in long simulations)
    * @note Function name "c_quat" historically refers to crystal-to-sample rotation
    * 
    * @see emap_to_quat for converting rotation vector to quaternion
    * @see quat_prod for underlying quaternion multiplication
    * @see quat_rel_rotation for computing incremental rotations
    */
   __ecmech_hdev__
   inline void get_c_quat(double* const c_quat,
                          const double* const dr_quat,
                          const double* const cn_quat) {
      // Compute : c = c_n * dr
      quat_prod(c_quat, cn_quat, dr_quat);
   }

   /**
    * @brief Convert quaternion to 3×3 rotation matrix.
    * 
    * Converts a unit quaternion to its corresponding orthogonal rotation matrix using \cite{kri-etal-94a}.
    * 
    * **Mathematical operation** (for q = [q0, q1, q2, q3]):
    * ```
    * C = | 1-2(q2²+q3²)   2(q1q2-q0q3)   2(q1q3+q0q2) |
    *     | 2(q1q2+q0q3)   1-2(q1²+q3²)   2(q2q3-q0q1) |
    *     | 2(q1q3-q0q2)   2(q2q3+q0q1)   1-2(q1²+q2²) |
    * ```
    * 
    * @param[out] c Rotation matrix (3×3 = 9 components, row-major order)
    * @param[in] quat Unit quaternion [q0, q1, q2, q3]
    * 
    * @note Result is orthogonal: C^T · C = I
    * @note det(C) = +1 (proper rotation, right-handed)
    * @note quat should be normalized: ‖quat‖ = 1
    */
   __ecmech_hdev__
   inline void quat_to_tensor(double* const c, // ndim * ndim
                              const double* const quat // qdim
                              ) {
      double x0sq = quat[0] * quat[0];
      double x1sq = quat[1] * quat[1];
      double x2sq = quat[2] * quat[2];
      double x3sq = quat[3] * quat[3];

      double x0x1 = quat[0] * quat[1];
      double x0x2 = quat[0] * quat[2];
      double x0x3 = quat[0] * quat[3];

      double x1x2 = quat[1] * quat[2];
      double x1x3 = quat[1] * quat[3];

      double x2x3 = quat[2] * quat[3];

      c[ECMECH_NN_INDX(0, 0, ndim)] = x0sq + x1sq - x2sq - x3sq;
      c[ECMECH_NN_INDX(0, 1, ndim)] = two * (x1x2 - x0x3);
      c[ECMECH_NN_INDX(0, 2, ndim)] = two * (x1x3 + x0x2);
      c[ECMECH_NN_INDX(1, 0, ndim)] = two * (x1x2 + x0x3);
      c[ECMECH_NN_INDX(1, 1, ndim)] = x0sq - x1sq + x2sq - x3sq;
      c[ECMECH_NN_INDX(1, 2, ndim)] = two * (x2x3 - x0x1);
      c[ECMECH_NN_INDX(2, 0, ndim)] = two * (x1x3 - x0x2);
      c[ECMECH_NN_INDX(2, 1, ndim)] = two * (x2x3 + x0x1);
      c[ECMECH_NN_INDX(2, 2, ndim)] = x0sq - x1sq - x2sq + x3sq;
   }

   /**
    * @brief Compute 5×5 rotation matrix for deviatoric tensors from 3×3 rotation matrix.
    * 
    * Constructs the 5×5 rotation matrix that transforms deviatoric 5-vectors
    * between lattice and sample frames. The matrix is derived from the 3×3
    * rotation matrix C(quat).
    * 
    * **Mathematical foundation**:
    * For deviatoric tensor D': D'_sample = C · D'_lattice · C^T
    * In 5-vector form: vecd_sample = R(5×5) · vecd_lattice
    * 
    * @param[out] qr5x5_raw Rotation matrix (5×5 = 25 components, row-major)
    * @param[in] c Rotation matrix (3×3 = 9 components, row-major order)
    * 
    * @note R is orthogonal: R^T · R = I (5×5 identity)
    * @note Preserves norm: ‖vecd_sample‖ = ‖vecd_lattice‖
    * @note Used in lattice rotation updates during crystal plasticity evolution
    */
   __ecmech_hdev__
   inline void get_rot_mat_vecd(double* const qr5x5_raw, // ntvec * ntvec
                                const double* const c // ndim * ndim
                                ) {
      // include "mc_vars.f90"
      // include "set_mc.f90"
#include "util/mc_vars_set.h"


      // IF ((UBOUND(c,DIM=1) /= DIMS) .OR. (UBOUND(c,DIM=2) /= DIMS)) &
      // & CALL consider_ierr(1,location,CIERR_DIMS_p,IERR_FATAL_p)
      // IF ((UBOUND(qr5x5,DIM=1) /= TVEC) .OR. (UBOUND(qr5x5,DIM=2) /= TVEC)) &
      // & CALL consider_ierr(1,location,CIERR_DIMS_p,IERR_FATAL_p)

      RAJA::View<double, RAJA::Layout<2> > qr5x5(qr5x5_raw, ecmech::ntvec, ecmech::ntvec);

      // ! if do not want to assume (c31**2+c32**2+c33**2=1)
      // qr5x5(1, 1)  =  c33 * c33 - onehalf * (c31 * c31 + c32 * c32)

      qr5x5(0, 0) = onehalf * (c11 * c11 - c12 * c12 - c21 * c21 + c22 * c22);
      qr5x5(0, 1) = sqr3 * onehalf * (c13 * c13 - c23 * c23);
      qr5x5(0, 2) = c11 * c12 - c21 * c22;
      qr5x5(0, 3) = c11 * c13 - c21 * c23;
      qr5x5(0, 4) = c12 * c13 - c22 * c23;
      qr5x5(1, 0) = sqr3 * onehalf * (c31 * c31 - c32 * c32);
      qr5x5(1, 1) = thrhalf * c33 * c33 - onehalf;
      qr5x5(1, 2) = sqr3 * c31 * c32;
      qr5x5(1, 3) = sqr3 * c31 * c33;
      qr5x5(1, 4) = sqr3 * c32 * c33;
      qr5x5(2, 0) = c11 * c21 - c12 * c22;
      qr5x5(2, 1) = sqr3 * c13 * c23;
      qr5x5(2, 2) = c11 * c22 + c12 * c21;
      qr5x5(2, 3) = c11 * c23 + c13 * c21;
      qr5x5(2, 4) = c12 * c23 + c13 * c22;
      qr5x5(3, 0) = c11 * c31 - c12 * c32;
      qr5x5(3, 1) = sqr3 * c13 * c33;
      qr5x5(3, 2) = c11 * c32 + c12 * c31;
      qr5x5(3, 3) = c11 * c33 + c13 * c31;
      qr5x5(3, 4) = c12 * c33 + c13 * c32;
      qr5x5(4, 0) = c21 * c31 - c22 * c32;
      qr5x5(4, 1) = sqr3 * c23 * c33;
      qr5x5(4, 2) = c21 * c32 + c22 * c31;
      qr5x5(4, 3) = c21 * c33 + c23 * c31;
      qr5x5(4, 4) = c22 * c33 + c23 * c32;
   }

   /** @} */ // end of quaternion_operations group

   /**
    * @brief Construct 3×5 matrix operator for matrix commutator [A,B] = A·B - B·A.
    * 
    * Computes the linear operator M such that the axial vector representation of
    * the commutator V = A·B - B·A can be obtained via matrix-vector multiplication:
    * 
    * **Mathematical operation**:
    * ```
    * w = M · a
    * ```
    * where:
    * - w = veccp(V) = axial vector of skew matrix V (3 components)
    * - V = A · B - B · A (matrix commutator, skew-symmetric result)
    * - a = vecds(A) = deviatoric representation of symmetric matrix A (5 components)
    * - b = vecds(B) = deviatoric representation of symmetric matrix B (5 components, input parameter)
    * - M is the 3×5 linear operator (depends on B)
    * 
    * **Physical interpretation**: In crystal plasticity, this computes how the plastic
    * spin W^p changes with respect to a symmetric tensor A (e.g., elastic strain),
    * given a fixed symmetric tensor B (e.g., elastic strain rate or another state variable).
    * 
    * **Key mathematical facts**:
    * - A and B are symmetric → V = A·B - B·A is **skew-symmetric**
    * - Skew-symmetric 3×3 matrix has 3 independent components
    * - These are represented as axial vector: w_i = ½ε_ijk V_jk
    * - The operation is LINEAR in A: V(αA + βA', B) = αV(A,B) + βV(A',B)
    * 
    * **Derivative interpretation**: This is effectively ∂(commutator)/∂A
    * - M represents the Jacobian of the commutator operation w.r.t. changes in A
    * - Used in implicit Jacobian assembly for spin-related derivatives
    * 
    * **Reversing the operation** (computing derivative w.r.t. B instead of A):
    * - To get M for V = B·A - A·B (negative commutator), pass cmv6a and negate result
    * - This gives: w' = -M'·b where M' is computed from A
    * 
    * @param[out] M35 Linear operator matrix (3×5 = 15 components, row-major)
    *                 - M35[i,j] at index ECMECH_NM_INDX(i, j, nwvec, ntvec)
    *                 - Maps 5-component deviatoric vector to 3-component spin vector
    * @param[in] cmv6b Deviatoric representation of symmetric matrix B (5 or 6 components)
    *                  - cmv6b[0:4] are the 5 deviatoric components
    *                  - cmv6b[iSvecS] (spherical/trace component) is NOT accessed
    *                  - Can pass 6-component svec; only first 5 used
    * 
    * @note A and B are symmetric but NOT necessarily deviatoric (can have trace)
    * @note Result V is always traceless (skew-symmetric → tr(V) = 0)
    * @note The commutator [A,B] measures non-commutativity of matrix multiplication
    * @note Used in computing elastic-elastic interaction contributions to plastic spin
    * 
    * **Typical usage context**:
    * - Implicit Jacobian: ∂W^p/∂ε^e where W^p involves commutators
    * - Elastic-elastic spin coupling: W^ee = ½[ε^e, dε^e/dt]
    * - Tangent stiffness assembly for crystal plasticity
    * 
    * @see skewToVeccp for converting skew matrix to axial vector
    * @see vecsMaTb for computing outer product (used to build full commutator)
    * @see matToPQ for decomposing matrices into symmetric/skew parts
    */
   __ecmech_hdev__
   inline void M35_d_AAoB_dA(double* const M35, // nwvec * ntvec
                             const double* const cmv6b // nsvec or ntvec -- cmv6b[iSvecS] not accessed
                             ) {
#include "util/vb_d_vars_set.h"

      M35[ECMECH_NM_INDX(0, 0, nwvec, ntvec)] = vb5 * onehalf;
      M35[ECMECH_NM_INDX(1, 0, nwvec, ntvec)] = vb4 * onehalf;
      M35[ECMECH_NM_INDX(2, 0, nwvec, ntvec)] = -vb3;
      M35[ECMECH_NM_INDX(0, 1, nwvec, ntvec)] = vb5 * halfsqr3;
      M35[ECMECH_NM_INDX(1, 1, nwvec, ntvec)] = -vb4 * halfsqr3;
      M35[ECMECH_NM_INDX(2, 1, nwvec, ntvec)] = zero;
      M35[ECMECH_NM_INDX(0, 2, nwvec, ntvec)] = -vb4 * onehalf;
      M35[ECMECH_NM_INDX(1, 2, nwvec, ntvec)] = vb5 * onehalf;
      M35[ECMECH_NM_INDX(2, 2, nwvec, ntvec)] = vb1;
      M35[ECMECH_NM_INDX(0, 3, nwvec, ntvec)] = vb3 * onehalf;
      M35[ECMECH_NM_INDX(1, 3, nwvec, ntvec)] = halfsqr3 * vb2 - vb1 * onehalf;
      M35[ECMECH_NM_INDX(2, 3, nwvec, ntvec)] = -vb5 * onehalf;
      M35[ECMECH_NM_INDX(0, 4, nwvec, ntvec)] = -vb1 * onehalf - halfsqr3 * vb2;
      M35[ECMECH_NM_INDX(1, 4, nwvec, ntvec)] = -vb3 * onehalf;
      M35[ECMECH_NM_INDX(2, 4, nwvec, ntvec)] = vb4 * onehalf;
      // M36[ECMECH_NM_INDX(0,5,nwvec,nsvec)] = zero ;
      // M36[ECMECH_NM_INDX(1,5,nwvec,nsvec)] = zero ;
      // M36[ECMECH_NM_INDX(2,5,nwvec,nsvec)] = zero ;
   }

   // included for documentation purposes:
   /*
       SUBROUTINE rot_mat_wveccp(& ! NOT same as rot_mat_skew
          &   c, qr3x3&
          &   )

       ! NOTE: cross-product notation

       IMPLICIT NONE

   #ifdef DO_CHECKS_ALL
       REAL(idp), INTENT(in) :: c(:,:)
       REAL(idp) :: qr3x3(:,:)
       IF ((UBOUND(c,DIM=1) /= DIMS) .OR. (UBOUND(c,DIM=2) /= DIMS)) &
            & CALL consider_ierr(1,location,CIERR_DIMS_p,IERR_FATAL_p)
       IF ((UBOUND(qr3x3,DIM=1) /= DIMS) .OR. (UBOUND(qr3x3,DIM=2) /= DIMS)) &
            & CALL consider_ierr(1,location,CIERR_DIMS_p,IERR_FATAL_p)
   #else
       REAL(idp), INTENT(in) :: c(DIMS, DIMS)
       REAL(idp) :: qr3x3(DIMS, DIMS)
   #endif
       !
       !     Construct 3X3 rotation matrix for skew 2nd order tensors
       !     [W]_sm = [c] [W]_lat [c]'  <=>  {W}_sm = [qr3x3] {W}_lat
       !
       ! for notation w_i = epsilon_jik W_jk, [qr3x3] = [c];

       qr3x3 = c

     END SUBROUTINE rot_mat_wveccp
    */

   /**
    * @brief Derivative of quaternion w.r.t. exponential map (transposed).
    * 
    * Computes ∂quat/∂emap in transposed form for Jacobian assembly.
    * 
    * **Mathematical operation**:
    * ```
    * ∂quat[i] / ∂emap[j]  for i ∈ {0,1,2,3}, j ∈ {0,1,2}
    * ```
    * 
    * @param[out] dq_demap_T_raw Jacobian matrix (4×3 = 12 components, transposed)
    * @param[in] emap Exponential map (angle-axis, length 3)
    * 
    * @note Used in chain rule: ∂f/∂emap = (∂f/∂quat) · (∂quat/∂emap)
    */
   __ecmech_hdev__
   inline
   void dquat_demap_T(double* const dqdeT_raw, // (EMAPDIM_p,QDIM_p)
                      const double* const emap // (EMAPDIM_p)
                      )
   {
      const double theta_sm_a = 1e-9;
      const double oo48 = 1.0 / 48.0;


      double theta = vecNorm<emapdim>(emap);

      double theta_inv, sthhbyth, halfsthh, na, nb, nc;
      if (fabs(theta) < theta_sm_a) {
         sthhbyth = onehalf - theta * theta * oo48; // truncated Taylor seriers; probably safe to just use onehalf and be done with it
         halfsthh = theta * oneqrtr; // truncated Taylor seriers
         if (fabs(theta) < idp_tiny_sqrt) {
            // n is arbitrary, as theta is effectively zero
            na = one; nb = zero; nc = zero;
         }
         else {
            theta_inv = one / theta;
            na = emap[0] * theta_inv; nb = emap[1] * theta_inv; nc = emap[2] * theta_inv;
         }
      }
      else {
         halfsthh = sin(theta * onehalf);
         sthhbyth = halfsthh / theta;
         halfsthh = halfsthh * onehalf;
         theta_inv = one / theta;
         na = emap[0] * theta_inv; nb = emap[1] * theta_inv; nc = emap[2] * theta_inv;
      }
      //
      double halfcthh = cos(theta * onehalf) * onehalf;
      //
      // now have: halfsthh, sthhbyth, halfcthh, theta, na, nb, nc

      RAJA::View<double, RAJA::Layout<2> > dqdeT(dqdeT_raw, ecmech::emapdim, ecmech::qdim);

      dqdeT(0, 0) = -halfsthh * na;
      dqdeT(1, 0) = -halfsthh * nb;
      dqdeT(2, 0) = -halfsthh * nc;

      double temp = na * na;
      dqdeT(0, 1) = halfcthh * temp + sthhbyth * (one - temp);
      //
      temp = nb * nb;
      dqdeT(1, 2) = halfcthh * temp + sthhbyth * (one - temp);
      //
      temp = nc * nc;
      dqdeT(2, 3) = halfcthh * temp + sthhbyth * (one - temp);

      temp = halfcthh - sthhbyth;
      //
      double tempb = temp * na * nb;
      dqdeT(1, 1) = tempb;
      dqdeT(0, 2) = tempb;
      //
      tempb = temp * na * nc;
      dqdeT(2, 1) = tempb;
      dqdeT(0, 3) = tempb;
      //
      tempb = temp * nb * nc;
      dqdeT(2, 2) = tempb;
      dqdeT(1, 3) = tempb;
   }

   /**
    * @brief Derivative of rotation matrix (3×3) w.r.t. quaternion components.
    * 
    * Computes ∂C/∂quat where C is the 3×3 rotation matrix from quat_to_tensor.
    * 
    * **Mathematical operation**:
    * ```
    * dCmatx_dq[i,j,k] = ∂C[i,j] / ∂quat[k]  for i,j,k ∈ {0,1,2,3}
    * ```
    * 
    * @param[out] dCmatx_dq Derivative tensor (3×3×4 = 36 components)
    * @param[in] quat Unit quaternion [q0, q1, q2, q3]
    * 
    * @note Less commonly used than d_rot_mat_vecd variants (those work on 5-vectors)
    */
   __ecmech_hdev__
   inline void d_quat_to_tensor(double* const dcdq_raw, // (DIMS,DIMS,QDIM_p)
                                const double* const quat // (QDIM_p)
                                )
   {
      double tqa = two * quat[0];
      double tqb = two * quat[1];
      double tqc = two * quat[2];
      double tqd = two * quat[3];

      RAJA::View<double, RAJA::Layout<3> > dcdq(dcdq_raw, ecmech::ndim, ecmech::ndim, ecmech::qdim);

      // c(1,1) = x1sq+x2sq-x3sq-x4sq
      dcdq(0, 0, 0) = tqa;
      dcdq(0, 0, 1) = tqb;
      dcdq(0, 0, 2) = -tqc;
      dcdq(0, 0, 3) = -tqd;

      // c(0,1) = two*(x1x2-x0x3)
      dcdq(0, 1, 0) = -tqd;
      dcdq(0, 1, 1) = tqc;
      dcdq(0, 1, 2) = tqb;
      dcdq(0, 1, 3) = -tqa;

      // c(0,2) = two*(x1x3+x0x2)
      dcdq(0, 2, 0) = tqc;
      dcdq(0, 2, 1) = tqd;
      dcdq(0, 2, 2) = tqa;
      dcdq(0, 2, 3) = tqb;

      // c(1,0) = two*(x1x2+x0x3)
      dcdq(1, 0, 0) = tqd;
      dcdq(1, 0, 1) = tqc;
      dcdq(1, 0, 2) = tqb;
      dcdq(1, 0, 3) = tqa;

      // c(1,1) = x0sq-x1sq+x2sq-x3sq
      dcdq(1, 1, 0) = tqa;
      dcdq(1, 1, 1) = -tqb;
      dcdq(1, 1, 2) = tqc;
      dcdq(1, 1, 3) = -tqd;

      // c(1,2) = two*(x2x3-x0x1)
      dcdq(1, 2, 0) = -tqb;
      dcdq(1, 2, 1) = -tqa;
      dcdq(1, 2, 2) = tqd;
      dcdq(1, 2, 3) = tqc;

      // c(2,0) = two*(x1x3-x0x2)
      dcdq(2, 0, 0) = -tqc;
      dcdq(2, 0, 1) = tqd;
      dcdq(2, 0, 2) = -tqa;
      dcdq(2, 0, 3) = tqb;

      // c(2,1) = two*(x2x3+x0x1)
      dcdq(2, 1, 0) = tqb;
      dcdq(2, 1, 1) = tqa;
      dcdq(2, 1, 2) = tqd;
      dcdq(2, 1, 3) = tqc;

      // c(2,2) = x0sq-x1sq-x2sq+x3sq
      dcdq(2, 2, 0) = tqa;
      dcdq(2, 2, 1) = -tqb;
      dcdq(2, 2, 2) = -tqc;
      dcdq(2, 2, 3) = tqd;
   }

   /**
    * @brief Derivative of lattice-frame vector w.r.t. rotation matrix components.
    * 
    * Computes the Jacobian of the transformation from sample frame to lattice frame
    * with respect to the 3×3 rotation matrix C. Given a sample-frame deviatoric
    * vector, this computes how the corresponding lattice-frame vector changes as
    * the rotation matrix C varies.
    * 
    * **Transformation**: vec_lat = Q^T · vec_sm where Q is the 5×5 rotation matrix
    * derived from the 3×3 rotation matrix C.
    * 
    * **Mathematical operation**:
    * ```
    * dvdc[i,k,l] = ∂(vec_lat[i]) / ∂(C[k,l])
    * ```
    * for i ∈ {0,1,2,3,4} (deviatoric components)
    *     k,l ∈ {0,1,2} (rotation matrix indices)
    * 
    * **Usage context**: Used in implicit Jacobian assembly to compute sensitivities
    * of lattice-frame quantities (stress, strain) with respect to orientation changes.
    * 
    * @param[out] dvdc_raw Jacobian tensor (5×3×3 = 45 components)
    *                      - dvdc[i,k,l] at index i*9 + k*3 + l
    * @param[in] c Rotation matrix C (3×3 = 9 components, row-major)
    *              - Crystal-to-sample rotation matrix
    * @param[in] vec_sm Sample-frame deviatoric vector (5 components)
    *                   - The vector being transformed to lattice frame
    * 
    * @note This is "latop" = lattice operation (computing lattice-frame result)
    * @note Pairs with d_rot_mat_vecd_smop (sample operation)
    * @note Uses util/mc_vars_set.h and util/vad_vars_set.h for component extraction
    * 
    * @see d_rot_mat_vecd_smop for the inverse transformation (lattice to sample)
    * @see get_rot_mat_vecd for computing the 5×5 rotation matrix itself
    */
   __ecmech_hdev__
   inline void d_rot_mat_vecd_latop(double* const dvdc_raw, // (TVEC,DIMS,DIMS)
                                    const double* const c, // (DIMS, DIMS)
                                    const double* const vec_sm // (TVEC)
                                    )
   {
#include "util/mc_vars_set.h"
#include "util/vad_vars_set.h"

      // include "d_Alat_dC.f90"
      RAJA::View<double, RAJA::Layout<3> > dvdc(dvdc_raw, ecmech::ntvec, ecmech::ndim, ecmech::ndim);
      dvdc(0, 0, 0) = c11 * va1 - c11 * sqr3 * va2 * onethird + c21 * va3 + c31 * va4;
      dvdc(1, 0, 0) = -sqr3 * (3 * c11 * va1 - c11 * sqr3 * va2 + 3 * c21 * va3 + 3 * c31 * va4) * oneninth;
      dvdc(2, 0, 0) = c12 * va1 - c12 * sqr3 * va2 * onethird + c22 * va3 + c32 * va4;
      dvdc(3, 0, 0) = c13 * va1 - c13 * sqr3 * va2 * onethird + c23 * va3 + c33 * va4;
      dvdc(4, 0, 0) = zero;
      dvdc(0, 1, 0) = c11 * va3 - c21 * va1 - c21 * sqr3 * va2 * onethird + c31 * va5;
      dvdc(1, 1, 0) = sqr3 * (-3 * c11 * va3 + 3 * c21 * va1 + c21 * sqr3 * va2 - 3 * c31 * va5) * oneninth;
      dvdc(2, 1, 0) = c12 * va3 - c22 * va1 - c22 * sqr3 * va2 * onethird + c32 * va5;
      dvdc(3, 1, 0) = c13 * va3 - c23 * va1 - c23 * sqr3 * va2 * onethird + c33 * va5;
      dvdc(4, 1, 0) = zero;
      dvdc(0, 2, 0) = c11 * va4 + c21 * va5 + twothird * c31 * sqr3 * va2;
      dvdc(1, 2, 0) = -sqr3 * c11 * va4 * onethird - sqr3 * c21 * va5 * onethird - twothird * c31 * va2;
      dvdc(2, 2, 0) = c12 * va4 + c22 * va5 + twothird * c32 * sqr3 * va2;
      dvdc(3, 2, 0) = c13 * va4 + c23 * va5 + twothird * c33 * sqr3 * va2;
      dvdc(4, 2, 0) = zero;
      dvdc(0, 0, 1) = -c12 * va1 + c12 * sqr3 * va2 * onethird - c22 * va3 - c32 * va4;
      dvdc(1, 0, 1) = sqr3 * (-3 * c12 * va1 + c12 * sqr3 * va2 - 3 * c22 * va3 - 3 * c32 * va4) * oneninth;
      dvdc(2, 0, 1) = c11 * va1 - c11 * sqr3 * va2 * onethird + c21 * va3 + c31 * va4;
      dvdc(3, 0, 1) = zero;
      dvdc(4, 0, 1) = c13 * va1 - c13 * sqr3 * va2 * onethird + c23 * va3 + c33 * va4;
      dvdc(0, 1, 1) = -c12 * va3 + c22 * va1 + c22 * sqr3 * va2 * onethird - c32 * va5;
      dvdc(1, 1, 1) = sqr3 * (-3 * c12 * va3 + 3 * c22 * va1 + c22 * sqr3 * va2 - 3 * c32 * va5) * oneninth;
      dvdc(2, 1, 1) = c11 * va3 - c21 * va1 - c21 * sqr3 * va2 * onethird + c31 * va5;
      dvdc(3, 1, 1) = zero;
      dvdc(4, 1, 1) = c13 * va3 - c23 * va1 - c23 * sqr3 * va2 * onethird + c33 * va5;
      dvdc(0, 2, 1) = -c12 * va4 - c22 * va5 - twothird * c32 * sqr3 * va2;
      dvdc(1, 2, 1) = -sqr3 * c12 * va4 * onethird - sqr3 * c22 * va5 * onethird - twothird * c32 * va2;
      dvdc(2, 2, 1) = c11 * va4 + c21 * va5 + twothird * c31 * sqr3 * va2;
      dvdc(3, 2, 1) = zero;
      dvdc(4, 2, 1) = c13 * va4 + c23 * va5 + twothird * c33 * sqr3 * va2;
      dvdc(0, 0, 2) = zero;
      dvdc(1, 0, 2) = -twothird * sqr3i * (-3 * c13 * va1 + sqr3 * c13 * va2 - 3 * c23 * va3 - 3 * c33 * va4);
      dvdc(2, 0, 2) = zero;
      dvdc(3, 0, 2) = c11 * va1 - c11 * sqr3 * va2 * onethird + c21 * va3 + c31 * va4;
      dvdc(4, 0, 2) = c12 * va1 - c12 * sqr3 * va2 * onethird + c22 * va3 + c32 * va4;
      dvdc(0, 1, 2) = zero;
      dvdc(1, 1, 2) = -twothird * sqr3i * (-3 * va3 * c13 + 3 * va1 * c23 + sqr3 * c23 * va2 - 3 * va5 * c33);
      dvdc(2, 1, 2) = zero;
      dvdc(3, 1, 2) = c11 * va3 - c21 * va1 - c21 * sqr3 * va2 * onethird + c31 * va5;
      dvdc(4, 1, 2) = c12 * va3 - c22 * va1 - c22 * sqr3 * va2 * onethird + c32 * va5;
      dvdc(0, 2, 2) = zero;
      dvdc(1, 2, 2) = twothird * sqr3 * c13 * va4 + twothird * sqr3 * c23 * va5 + fourthirds * c33 * va2;
      dvdc(2, 2, 2) = zero;
      dvdc(3, 2, 2) = c11 * va4 + c21 * va5 + twothird * c31 * sqr3 * va2;
      dvdc(4, 2, 2) = c12 * va4 + c22 * va5 + twothird * c32 * sqr3 * va2;
   } // d_rot_mat_vecd_latop

   /**
    * @brief Derivative of lattice-frame spin vector w.r.t. rotation matrix components.
    * 
    * Computes the Jacobian of the transformation of a spin (skew-symmetric tensor)
    * from sample frame to lattice frame with respect to the 3×3 rotation matrix C.
    * Spin vectors use 3-component axial representation (cross-product notation).
    * 
    * **Transformation**: vec_lat = C^T · vec_sm for spin vectors in axial form
    * 
    * For skew-symmetric tensors W: W_lat = C^T · W_sm · C
    * In axial vector form: w_lat = C^T · w_sm
    * 
    * **Mathematical operation**:
    * ```
    * dvdc[i,k,l] = ∂(vec_lat[i]) / ∂(C[k,l])
    * ```
    * for i ∈ {0,1,2} (spin vector components)
    *     k,l ∈ {0,1,2} (rotation matrix indices)
    * 
    * **Implementation detail**: The rotation matrix parameter `c` is commented out
    * in the signature because the derivative formula turns out to NOT depend on
    * the actual values of C - only on the structure of the transformation.
    * 
    * @param[out] dvdc_raw Jacobian tensor (3×3×3 = 27 components)
    *                      - dvdc[i,k,l] at index i*9 + k*3 + l
    * @param[in] cmv3w Sample-frame spin vector (3 components)
    *                  - Axial vector representation of skew-symmetric tensor
    *                  - Comment indicates this is vec_sm(WVEC)
    * 
    * @note Parameter name "cmv3w" is legacy notation
    * @note The rotation matrix C is NOT needed for this derivative (commented out)
    * @note "wveccp" refers to spin vector in cross-product (axial) notation
    * @note "latop" = lattice operation (computing lattice-frame result)
    * @note Uses util/vw_vars_set.h for spin component extraction
    * 
    * @see d_rot_mat_vecd_latop for deviatoric tensor version
    * @see skewToVeccp for converting skew matrix to axial vector
    */
   __ecmech_hdev__
   inline void d_rot_mat_wveccp_latop(double* const dvdc_raw, // (WVEC,DIMS,DIMS)
                                      // const double* const c, // (DIMS, DIMS) // not used
                                      const double* const cmv3w // (WVEC) // vec_sm(WVEC)
                                      )
   {
#include "util/vw_vars_set.h"

      // include "d_Wlat_dC.f90"
      RAJA::View<double, RAJA::Layout<3> > dvdc(dvdc_raw, ecmech::nwvec, ecmech::ndim, ecmech::ndim);

      dvdc(0, 0, 0) = vw1;
      dvdc(1, 0, 0) = zero;
      dvdc(2, 0, 0) = zero;
      dvdc(0, 1, 0) = vw2;
      dvdc(1, 1, 0) = zero;
      dvdc(2, 1, 0) = zero;
      dvdc(0, 2, 0) = vw3;
      dvdc(1, 2, 0) = zero;
      dvdc(2, 2, 0) = zero;
      dvdc(0, 0, 1) = zero;
      dvdc(1, 0, 1) = vw1;
      dvdc(2, 0, 1) = zero;
      dvdc(0, 1, 1) = zero;
      dvdc(1, 1, 1) = vw2;
      dvdc(2, 1, 1) = zero;
      dvdc(0, 2, 1) = zero;
      dvdc(1, 2, 1) = vw3;
      dvdc(2, 2, 1) = zero;
      dvdc(0, 0, 2) = zero;
      dvdc(1, 0, 2) = zero;
      dvdc(2, 0, 2) = vw1;
      dvdc(0, 1, 2) = zero;
      dvdc(1, 1, 2) = zero;
      dvdc(2, 1, 2) = vw2;
      dvdc(0, 2, 2) = zero;
      dvdc(1, 2, 2) = zero;
      dvdc(2, 2, 2) = vw3;
   } // d_rot_mat_wveccp_latop

   /**
    * @brief Compute Jacobian derivatives for implicit quaternion-based crystal plasticity.
    * 
    * Evaluates all derivative terms needed for the implicit Jacobian when solving
    * the coupled system of elastic strain and lattice orientation using quaternion
    * representation. Computes how the applied (apparent) deformation rate and spin
    * change with respect to the rotation increment variable xi (exponential map).
    * 
    * **Physical context**: In implicit crystal plasticity, the solve variables include
    * the rotation increment xi (exponential map). This function provides the Jacobian
    * contributions relating changes in xi to changes in the velocity gradient quantities.
    * 
    * **Mathematical operations**:
    * ```
    * 1. dxtal_ori_quat_dxi_T = ∂quat / ∂xi  (transposed)
    * 2. dDapp_dxi = ∂(D_apparent) / ∂xi  (applied deformation rate)
    * 3. dWapp_dxi = ∂(W_apparent) / ∂xi  (applied spin)
    * ```
    * 
    * **Derivative chain**:
    * ```
    * xi (rotation increment) 
    *   → Δquat (via emap_to_quat)
    *   → quat_{n+1} (via get_c_quat with quat_n)
    *   → C (rotation matrix via quat_to_tensor)
    *   → Q (5×5 rotation matrix via get_rot_mat_vecd)
    *   → transformed velocity gradient quantities
    * ```
    * 
    * **Implementation approach**:
    * - Computes dquat/dxi using dquat_demap_T
    * - Uses chain rule through quaternion product
    * - Transforms to rotation matrix derivatives
    * - Applies to deformation rate and spin
    * 
    * @param[out] dxtal_ori_quat_dxi_T Derivative of orientation quaternion w.r.t. xi (3×4 = 12 components, transposed)
    *                                   - dquat[j]/dxi[i] at index i*4 + j
    *                                   - Transposed for efficient matrix multiplication
    * @param[out] dDapp_dxi Derivative of applied deformation rate w.r.t. xi (5×3 = 15 components)
    *                       - dD[i]/dxi[j] at index i*3 + j
    *                       - Only deviatoric part (5 components), trace unchanged by rotation
    * @param[out] dWapp_dxi Derivative of applied spin w.r.t. xi (3×3 = 9 components)
    *                       - dW[i]/dxi[j] at index i*3 + j
    * @param[in] def_rate_d5_sample Deformation rate in sample frame (5 components)
    *                               - Can be TVEC (deviatoric) or SVEC (with trace) - only deviatoric used
    * @param[in] w_vec_sm Spin vector in sample frame (3 components)
    *                     - Axial representation of skew-symmetric spin tensor
    * @param[in] xi Rotation increment in exponential map form (3 components)
    *               - Current value of solve variable (angle-axis representation)
    * @param[in] xtal_ori_quat_n Crystal orientation quaternion at previous step t_n (4 components)
    *                            - Orientation at beginning of time step
    * @param[in] xtal_rmat Rotation matrix at current iteration (3×3 = 9 components)
    *                      - Computed from current quaternion estimate
    * @param[in] xtal_ori_quat Crystal orientation quaternion at current iteration (4 components)
    *                          - Current estimate: quat = quat_n ⊗ exp(xi)
    * 
    * @note "dxi" refers to derivative with respect to exponential map xi (rotation increment)
    * @note All rotations and derivatives are in the reference (constituent) frame
    * @note Only deviatoric part of D has derivatives (trace invariant under rotation)
    * @note Uses util headers (vad_vars_set.h, mc_vars_set.h) for component extraction
    * @note Commented-out parameter dxtal_rmat_dxi can be computed internally if needed
    * @note Commented-out parameter A_quat (Δquat) is computed internally
    * 
    * @see dquat_demap_T for quaternion-to-exponential map derivative
    * @see d_rot_mat_vecd_latop for rotation matrix derivatives
    * @see d_rot_mat_wveccp_latop for spin rotation derivatives
    * @see get_c_quat for quaternion composition
    */
   __ecmech_hdev__
   inline
   void eval_d_dxi_impl_quat(double* const dxtal_ori_quat_dxi_T, // (WVEC,QDIM_p)
                             // double* const dxtal_rmat_dxi, // (DIMS,DIMS,WVEC)
                             double* const dDapp_dxi, // dDapp_dxi(TVEC, WVEC)
                             double* const dWapp_dxi, // dWapp_dxi(WVEC, WVEC)
                             const double* const def_rate_d5_sample, // (TVEC), or (SVEC) is fine too
                             const double* const w_vec_sm, // (WVEC)
                             const double* const xi, // (WVEC)
                             const double* const xtal_ori_quat_n, // (QDIM_p)
                             const double* const xtal_rmat, // (DIMS,DIMS)
                             const double* const xtal_ori_quat // (QDIM_p)
                             // const double* const A_quat // (QDIM_p) // not used
                             ) {
      // working with quats, so do not call eval_d_cA_dxi(dc_dxi, dA_dxi, xi, c_n)
      //
      {
         double dA_quat_dxi_T[ ecmech::ndim * ecmech::qdim ]; // (QDIM_p,DIMS)^T
         dquat_demap_T(dA_quat_dxi_T, xi);

         // can get away with these three calls as quat_prod is bilinear in the input arguments
         //
         quat_prod(&(dxtal_ori_quat_dxi_T[ecmech::qdim * 0]), xtal_ori_quat_n, &(dA_quat_dxi_T[ecmech::qdim * 0]) );
         quat_prod(&(dxtal_ori_quat_dxi_T[ecmech::qdim * 1]), xtal_ori_quat_n, &(dA_quat_dxi_T[ecmech::qdim * 1]) );
         quat_prod(&(dxtal_ori_quat_dxi_T[ecmech::qdim * 2]), xtal_ori_quat_n, &(dA_quat_dxi_T[ecmech::qdim * 2]) );
      }
      // now have dxtal_ori_quat_dxi

      double dxtal_rmat_dxi[ (ecmech::ndim * ecmech::ndim) *ecmech::nwvec ]; // (DIMS,DIMS,WVEC)
      {
         double dCmatx_dq[ (ecmech::ndim * ecmech::ndim) *ecmech::qdim ]; // (DIMS,DIMS,QDIM_p)
         // get dxtal_rmat_dxi
         d_quat_to_tensor(dCmatx_dq, xtal_ori_quat);
         vecsMABT<ndim*ndim, nwvec, qdim>(dxtal_rmat_dxi, dCmatx_dq, dxtal_ori_quat_dxi_T); // vecsMABT because _T on dxtal_ori_quat_dxi_T
      }

      {
         double dD_dxtal_rmat[ ecmech::ntvec * (ecmech::ndim * ecmech::ndim) ];
         d_rot_mat_vecd_latop(dD_dxtal_rmat, xtal_rmat, def_rate_d5_sample);
         //
         vecsMAB<ntvec, nwvec, ndim*ndim>(dDapp_dxi, dD_dxtal_rmat, dxtal_rmat_dxi);
         // dDapp_dxi(SVEC,:) = zero
      }

      {
         double dW_dxtal_rmat[ ecmech::nwvec * (ecmech::ndim * ecmech::ndim) ]; // (WVEC,DIMS,DIMS)
         d_rot_mat_wveccp_latop(dW_dxtal_rmat, // xtal_rmat,
                                w_vec_sm);
         //
         vecsMAB<nwvec, nwvec, ndim*ndim>(dWapp_dxi, dW_dxtal_rmat, dxtal_rmat_dxi);
      }
   }

   /**
    * @brief Derivative of sample-frame vector w.r.t. rotation matrix components.
    * 
    * Computes the Jacobian of the transformation from lattice frame to sample frame
    * with respect to the 3×3 rotation matrix C. Given a lattice-frame deviatoric
    * vector, this computes how the corresponding sample-frame vector changes as
    * the rotation matrix C varies.
    * 
    * **Transformation**: vec_sm = Q · vec_lat where Q is the 5×5 rotation matrix
    * derived from the 3×3 rotation matrix C.
    * 
    * **Mathematical operation**:
    * ```
    * dvdc[i,k,l] = ∂(vec_sm[i]) / ∂(C[k,l])
    * ```
    * for i ∈ {0,1,2,3,4} (deviatoric components)
    *     k,l ∈ {0,1,2} (rotation matrix indices)
    * 
    * **Usage context**: Used in implicit Jacobian assembly to compute sensitivities
    * of sample-frame quantities with respect to orientation changes. Common in
    * stress update algorithms where lattice stress is rotated to sample frame.
    * 
    * @param[out] dvdc_raw Jacobian tensor (5×3×3 = 45 components)
    *                      - dvdc[i,k,l] at index i*9 + k*3 + l
    * @param[in] c Rotation matrix C (3×3 = 9 components, row-major)
    *              - Crystal-to-sample rotation matrix
    * @param[in] vec_lat Lattice-frame deviatoric vector (5 components)
    *                    - The vector being transformed to sample frame
    * 
    * @note This is "smop" = sample operation (computing sample-frame result)
    * @note Pairs with d_rot_mat_vecd_latop (lattice operation)
    * @note Uses util/mc_vars_set.h and util/vadl_vars_set.h for component extraction
    * 
    * @see d_rot_mat_vecd_latop for the inverse transformation (sample to lattice)
    * @see get_rot_mat_vecd for computing the 5×5 rotation matrix itself
    */
   __ecmech_hdev__
   inline void
   d_rot_mat_vecd_smop(double* const dvdc_raw, // (TVEC,DIMS,DIMS)
                       const double* const c, // (DIMS, DIMS)
                       const double* const vec_lat // (TVEC)
                       )
   {
      //
      // derivative of 5x5 rotation operation with respect to components of c
      //
      // {vec_sm} = [Q] {vec_lat}
      // dvdc is d({vec_sm})/d{C}
      //

#include "util/mc_vars_set.h"
#include "util/vadl_vars_set.h"

      // include "d_Asm_dC.f90"
      RAJA::View<double, RAJA::Layout<3> > dvdc(dvdc_raw, ecmech::ntvec, ecmech::ndim, ecmech::ndim);
      dvdc(0, 0, 0) = c11 * va1 - c11 * sqr3 * va2 * onethird + c12 * va3 + c13 * va4;
      dvdc(1, 0, 0) = sqr3 * (-3 * c11 * va1 + sqr3 * c11 * va2 - 3 * c12 * va3 - 3 * c13 * va4) * oneninth;
      dvdc(2, 0, 0) = c21 * va1 - c21 * sqr3 * va2 * onethird + c22 * va3 + c23 * va4;
      dvdc(3, 0, 0) = c31 * va1 - c31 * sqr3 * va2 * onethird + c32 * va3 + c33 * va4;
      dvdc(4, 0, 0) = zero;
      dvdc(0, 1, 0) = -c21 * va1 + c21 * sqr3 * va2 * onethird - c22 * va3 - c23 * va4;
      dvdc(1, 1, 0) = -sqr3 * (3 * c21 * va1 - c21 * sqr3 * va2 + 3 * c22 * va3 + 3 * c23 * va4) * oneninth;
      dvdc(2, 1, 0) = c11 * va1 - c11 * sqr3 * va2 * onethird + c12 * va3 + c13 * va4;
      dvdc(3, 1, 0) = zero;
      dvdc(4, 1, 0) = c31 * va1 - c31 * sqr3 * va2 * onethird + c32 * va3 + c33 * va4;
      dvdc(0, 2, 0) = zero;
      dvdc(1, 2, 0) = -twothird * sqr3i * (-3 * c31 * va1 + c31 * sqr3 * va2 - 3 * c32 * va3 - 3 * c33 * va4);
      dvdc(2, 2, 0) = zero;
      dvdc(3, 2, 0) = c11 * va1 - c11 * sqr3 * va2 * onethird + c12 * va3 + c13 * va4;
      dvdc(4, 2, 0) = c21 * va1 - c21 * sqr3 * va2 * onethird + c22 * va3 + c23 * va4;
      dvdc(0, 0, 1) = c11 * va3 - c12 * va1 - c12 * sqr3 * va2 * onethird + c13 * va5;
      dvdc(1, 0, 1) = sqr3 * (-3 * c11 * va3 + 3 * c12 * va1 + c12 * sqr3 * va2 - 3 * c13 * va5) * oneninth;
      dvdc(2, 0, 1) = c21 * va3 - c22 * va1 - c22 * sqr3 * va2 * onethird + c23 * va5;
      dvdc(3, 0, 1) = c31 * va3 - c32 * va1 - c32 * sqr3 * va2 * onethird + c33 * va5;
      dvdc(4, 0, 1) = zero;
      dvdc(0, 1, 1) = -c21 * va3 + c22 * va1 + c22 * sqr3 * va2 * onethird - c23 * va5;
      dvdc(1, 1, 1) = sqr3 * (-3 * c21 * va3 + 3 * c22 * va1 + sqr3 * va2 * c22 - 3 * c23 * va5) * oneninth;
      dvdc(2, 1, 1) = c11 * va3 - c12 * va1 - c12 * sqr3 * va2 * onethird + c13 * va5;
      dvdc(3, 1, 1) = zero;
      dvdc(4, 1, 1) = c31 * va3 - c32 * va1 - c32 * sqr3 * va2 * onethird + c33 * va5;
      dvdc(0, 2, 1) = zero;
      dvdc(1, 2, 1) = -twothird * sqr3i * (-3 * c31 * va3 + 3 * c32 * va1 + sqr3 * va2 * c32 - 3 * c33 * va5);
      dvdc(2, 2, 1) = zero;
      dvdc(3, 2, 1) = c11 * va3 - c12 * va1 - c12 * sqr3 * va2 * onethird + c13 * va5;
      dvdc(4, 2, 1) = c21 * va3 - c22 * va1 - c22 * sqr3 * va2 * onethird + c23 * va5;
      dvdc(0, 0, 2) = c11 * va4 + c12 * va5 + twothird * c13 * sqr3 * va2;
      dvdc(1, 0, 2) = -sqr3 * c11 * va4 * onethird - sqr3 * c12 * va5 * onethird - twothird * c13 * va2;
      dvdc(2, 0, 2) = c21 * va4 + c22 * va5 + twothird * c23 * sqr3 * va2;
      dvdc(3, 0, 2) = c31 * va4 + c32 * va5 + twothird * c33 * sqr3 * va2;
      dvdc(4, 0, 2) = zero;
      dvdc(0, 1, 2) = -c21 * va4 - c22 * va5 - twothird * c23 * sqr3 * va2;
      dvdc(1, 1, 2) = -sqr3 * c21 * va4 * onethird - sqr3 * c22 * va5 * onethird - twothird * c23 * va2;
      dvdc(2, 1, 2) = c11 * va4 + c12 * va5 + twothird * c13 * sqr3 * va2;
      dvdc(3, 1, 2) = zero;
      dvdc(4, 1, 2) = c31 * va4 + c32 * va5 + twothird * c33 * sqr3 * va2;
      dvdc(0, 2, 2) = zero;
      dvdc(1, 2, 2) = twothird * sqr3 * c31 * va4 + twothird * sqr3 * c32 * va5 + fourthirds * c33 * va2;
      dvdc(2, 2, 2) = zero;
      dvdc(3, 2, 2) = c11 * va4 + c12 * va5 + twothird * c13 * sqr3 * va2;
      dvdc(4, 2, 2) = c21 * va4 + c22 * va5 + twothird * c23 * sqr3 * va2;
   } // d_rot_mat_vecd_smop

   /** @} */ // end of rotation_derivatives group

   /**
    * @defgroup tangent_stiffness Tangent Stiffness Utilities
    * @brief Conversion and manipulation of material tangent matrices.
    * 
    * @{
    */

   /**
    * @brief Pre-multiply 6×n matrix by extended 5×5 rotation matrix.
    * 
    * Performs rotation of a 6×n matrix representing stress/strain quantities
    * using a 5×5 deviatoric rotation matrix. The 6×6 extended rotation is:
    * ```
    * [qr5x5,  0  ]
    * [  0  ,  1  ]
    * ```
    * where the last row/column (pressure/volumetric) is unchanged.
    * 
    * **Operations**:
    * - `l_T = false`: M_out = [qr5x5 ⊕ 1] · M_in  (standard rotation)
    * - `l_T = true`:  M_out = [qr5x5 ⊕ 1]^T · M_in  (transposed rotation)
    * 
    * **Physical interpretation**:
    * Rotates deviatoric stress/strain components (rows 0-4) while leaving
    * spherical/pressure component (row 5) invariant. Used for rotating tangent
    * stiffness matrices or stress tensors between crystallographic frames.
    * 
    * **Matrix dimensions**:
    * - M_in, M_out: 6×n matrices (row-major storage)
    * - qr5x5: 5×5 rotation matrix for deviatoric components
    * - Row 5 (index iSvecS) is copied unchanged from M_in to M_out
    * 
    * @tparam n Number of columns in input/output matrices
    * @tparam l_T Transpose flag
    *         - false: Apply rotation directly
    *         - true: Apply transposed rotation
    * @param[out] M_out Rotated matrix (6×n, row-major)
    * @param[in] M_in Input matrix (6×n, row-major)
    * @param[in] qr5x5 Rotation matrix for deviatoric components (5×5, row-major)
    * 
    * @note M_in and M_out MUST be different arrays (not safe for in-place operation)
    * @note Typically used with n=6 for rotating tangent stiffness matrices
    * @note The ⊕ symbol denotes block-diagonal extension: [A ⊕ 1] = diag(A, 1)
    * 
    * @see get_rot_mat_vecd for computing qr5x5 from rotation matrix
    * @see mtan_conv_sd_svec for converting tangent representations
    */
   template<int n, bool l_T>
   __ecmech_hdev__
   inline void
   qr6x6_pre_mul(double* const M_out, // 6xn
                 const double* const M_in, // 6xn
                 const double* const qr5x5 // 5x5
                 )
   {
      for (int iM = 0; iM<ecmech::ntvec; ++iM) { // only up to ntvec on purpose !
         for (int jM = 0; jM<n; ++jM) {
            int ijM = ECMECH_NM_INDX(iM, jM, ecmech::nsvec, n);
            M_out[ijM] = 0.0;
            for (int pTvec = 0; pTvec<ecmech::ntvec; ++pTvec) {
               if (l_T) {
                  M_out[ijM] += qr5x5[ECMECH_NN_INDX(pTvec, iM, ecmech::ntvec)] * M_in[ECMECH_NM_INDX(pTvec, jM, ecmech::nsvec, n)];
               }
               else {
                  M_out[ijM] += qr5x5[ECMECH_NN_INDX(iM, pTvec, ecmech::ntvec)] * M_in[ECMECH_NM_INDX(pTvec, jM, ecmech::nsvec, n)];
               }
            }
         }
      }

      for (int jM = 0; jM<n; ++jM) {
         int ijM = ECMECH_NM_INDX(iSvecS, jM, ecmech::nsvec, n);
         M_out[ijM] = M_in[ijM];
      }
   } // qr6x6_pre_mul


   /**
    * @brief Convert tangent stiffness from deviatoric-pressure to Voigt form.
    * 
    * Transforms the material tangent stiffness matrix from the internal
    * deviatoric-pressure (vecds) representation to standard Voigt (svec) form.
    * The transformation accounts for normalization differences and optionally
    * applies shear strain convention scaling.
    * 
    * **Transformation**: mtanSD = T · mtanSD_vecds · T^{-1}
    * 
    * where T is the transformation matrix between vecds and svec representations.
    * 
    * **Template parameter behavior**:
    * - `l_ddsdde_gamma = true`: Tangent is dσ/dγ (shear strain γ), applies √2/2 scaling to shear components
    * - `l_ddsdde_gamma = false`: Tangent is dσ/dε (engineering strain ε), applies √2 scaling to shear components
    * 
    * **Physical interpretation**:
    * The factor of 1/2 difference accounts for the relationship between engineering shear strain
    * (γ = 2ε₁₂) and tensorial shear strain (ε₁₂). This is critical for correct finite element
    * implementations where different strain conventions may be used.
    * 
    * **Representation details**:
    * - vecds: [dev₀, dev₁, dev₂, dev₃, dev₄, spherical] (6 components)
    * - svec: [σ₁₁, σ₂₂, σ₃₃, σ₂₃, σ₁₃, σ₁₂] (6 components, Voigt notation)
    * 
    * @tparam l_ddsdde_gamma Strain convention flag
    *         - true: Tangent w.r.t. shear strain γ (extra 1/2 factor for shear)
    *         - false: Tangent w.r.t. engineering strain ε (standard √2 factor)
    * @param[out] mtanSD_raw Material tangent in Voigt form (6×6 = 36 components, row-major)
    *                        - Component [i,j] at index i*6 + j
    * @param[in] mtanSD_vecds_raw Material tangent in deviatoric-spherical form (6×6 = 36 components, row-major)
    *                             - Input in vecds representation
    * 
    * @note Preserves stress-strain relationship: dσ = C : dε
    * @note The template parameter must be known at compile time
    * @note Used in finite element assembly to match host code strain conventions
    * 
    * @see qr6x6_pre_mul for rotating tangent stiffness matrices
    */
   template<bool l_ddsdde_gamma>
   __ecmech_hdev__
   inline
   void
   mtan_conv_sd_svec(double* const mtanSD_raw,
                     const double* const mtanSD_vecds_raw) {
      double C_raw[ecmech::nsvec2];
      double t1_vec[ecmech::nsvec], t2_vec[ecmech::nsvec], t3_vec[ecmech::nsvec];

      RAJA::View<double, RAJA::Layout<2> > mtanSD(mtanSD_raw, ecmech::nsvec, ecmech::nsvec);
      RAJA::View<double const, RAJA::Layout<2> > mtanSD_vecds(mtanSD_vecds_raw, ecmech::nsvec, ecmech::nsvec);
      RAJA::View<double, RAJA::Layout<2> > C(C_raw, ecmech::nsvec, ecmech::nsvec);

      // mtanSD = T . mtanSD_vecds . T^{-1}

      // C = T . mtanSD_vecds
      // C(i,:) = T(i,k) . mtanSD_vecds(k,:) -- sum over k
      //
      for (int jSvec = 0; jSvec<ecmech::nsvec; ++jSvec) {
         t3_vec[jSvec] = sqr3i * mtanSD_vecds(iSvecS, jSvec);
      }

      for (int jSvec = 0; jSvec<ecmech::nsvec; ++jSvec) {
         t1_vec[jSvec] = sqr2i * mtanSD_vecds(0, jSvec);
      }

      for (int jSvec = 0; jSvec<ecmech::nsvec; ++jSvec) {
         t2_vec[jSvec] = sqr6i * mtanSD_vecds(1, jSvec);
      }

      //
      for (int jSvec = 0; jSvec<ecmech::nsvec; ++jSvec) {
         C(0, jSvec) = t1_vec[jSvec] - t2_vec[jSvec] + t3_vec[jSvec];
      }

      for (int jSvec = 0; jSvec<ecmech::nsvec; ++jSvec) {
         C(1, jSvec) = -t1_vec[jSvec] - t2_vec[jSvec] + t3_vec[jSvec];
      }

      for (int jSvec = 0; jSvec<ecmech::nsvec; ++jSvec) {
         C(2, jSvec) = sqr2b3 * mtanSD_vecds(1, jSvec) + t3_vec[jSvec];
      }

      for (int jSvec = 0; jSvec<ecmech::nsvec; ++jSvec) {
         C(3, jSvec) = sqr2i * mtanSD_vecds(4, jSvec);
      }

      for (int jSvec = 0; jSvec<ecmech::nsvec; ++jSvec) {
         C(4, jSvec) = sqr2i * mtanSD_vecds(3, jSvec);
      }

      for (int jSvec = 0; jSvec<ecmech::nsvec; ++jSvec) {
         C(5, jSvec) = sqr2i * mtanSD_vecds(2, jSvec);
      }

      // mtanSD = C . T^{-1}
      // mtanSD(:,j) = C(:,k) . [T^{-1}](k,j) -- sum over k
      //
      for (int jSvec = 0; jSvec<ecmech::nsvec; ++jSvec) {
         t3_vec[jSvec] = C(jSvec, iSvecS) * sqr3i;
      }

      for (int jSvec = 0; jSvec<ecmech::nsvec; ++jSvec) {
         t1_vec[jSvec] = C(jSvec, 0) * sqr2i;
      }

      for (int jSvec = 0; jSvec<ecmech::nsvec; ++jSvec) {
         t2_vec[jSvec] = C(jSvec, 1) * sqr6i;
      }

      //
      for (int jSvec = 0; jSvec<ecmech::nsvec; ++jSvec) {
         mtanSD(jSvec, 0) = t1_vec[jSvec] - t2_vec[jSvec] + t3_vec[jSvec];
      }

      for (int jSvec = 0; jSvec<ecmech::nsvec; ++jSvec) {
         mtanSD(jSvec, 1) = -t1_vec[jSvec] - t2_vec[jSvec] + t3_vec[jSvec];
      }

      for (int jSvec = 0; jSvec<ecmech::nsvec; ++jSvec) {
         mtanSD(jSvec, 2) = sqr2b3 * C(jSvec, 1) + t3_vec[jSvec];
      }

      if (l_ddsdde_gamma) {
         // extra factor of 1/2 for shear deformation transformation
         for (int jSvec = 0; jSvec<ecmech::nsvec; ++jSvec) {
            mtanSD(jSvec, 3) = C(jSvec, 4) * sqr2i;
         }

         for (int jSvec = 0; jSvec<ecmech::nsvec; ++jSvec) {
            mtanSD(jSvec, 4) = C(jSvec, 3) * sqr2i;
         }

         for (int jSvec = 0; jSvec<ecmech::nsvec; ++jSvec) {
            mtanSD(jSvec, 5) = C(jSvec, 2) * sqr2i;
         }
      }
      else {
         for (int jSvec = 0; jSvec<ecmech::nsvec; ++jSvec) {
            mtanSD(jSvec, 3) = C(jSvec, 4) * sqr2;
         }

         for (int jSvec = 0; jSvec<ecmech::nsvec; ++jSvec) {
            mtanSD(jSvec, 4) = C(jSvec, 3) * sqr2;
         }

         for (int jSvec = 0; jSvec<ecmech::nsvec; ++jSvec) {
            mtanSD(jSvec, 5) = C(jSvec, 2) * sqr2;
         }
      }
   } // mtan_conv_sd_svec

   /** @} */ // end of tangent_stiffness group

   /**
    * @defgroup miller_indices Miller Index Conversions
    * @brief Convert crystallographic Miller indices to Cartesian coordinates.
    * 
    * For hexagonal crystals, Miller-Bravais indices use 4-index notation:
    * - Directions: [uvtw] where u+v+t = 0
    * - Planes: (hkil) where h+k+i = 0
    * 
    * These functions convert to orthogonal Cartesian coordinates using the
    * c/a ratio of the hexagonal lattice.
    * 
    * @{
    */

   /**
    * @brief Convert Miller direction [uvtw] to Cartesian vector.
    * 
    * Transforms a 4-index Miller direction to Cartesian coordinates for
    * hexagonal crystal systems.
    * 
    * **Coordinate system**:
    * - x-axis: Along a₁ basal vector
    * - y-axis: 30° from a₂
    * - z-axis: Along c-axis
    * 
    * @param[out] orthog Cartesian direction vector (length 3)
    * @param[in] miller Miller indices [u, v, t, w] (length 4, u+v+t=0)
    * @param[in] cOverA Lattice parameter ratio c/a
    * 
    * @note Output is normalized to unit length
    * @note For cubic: cOverA parameter ignored, use standard conversion
    */
   __ecmech_hdev__
   inline
   void
   m_to_o_dir(double* const dir_o, // ecmech::ndim
              const double* const dir_m, // ecmech::nMiller
              double cOverA
              ) {
      // note: this does not assume SUM(dir_m[:]) = 0
      dir_o[0] = dir_m[0] - onehalf * (dir_m[1] + dir_m[2]);
      dir_o[1] = sqr3 * onehalf * (dir_m[1] - dir_m[2]);
      dir_o[2] = dir_m[3] * cOverA;
   } // m_to_o_dir

   /**
    * @brief Convert Miller indices to orthogonal vectors for slip system definition.
    * 
    * Transforms a slip plane (hkil) and slip direction [uvtw] from Miller-Bravais
    * 4-index notation to Cartesian unit vectors for hexagonal crystal systems.
    * This is the single-system version used to define individual slip systems.
    * 
    * **Miller-Bravais indexing** (HCP crystals):
    * - Planes: (hkil) where h+k+i=0, i=-(h+k)
    * - Directions: [uvtw] where u+v+t=0, t=-(u+v)
    * 
    * **Coordinate system** (orthogonal basis):
    * - x-axis: Along a₁ basal vector
    * - y-axis: 30° from a₂ basal vector  
    * - z-axis: Along c-axis (perpendicular to basal plane)
    * 
    * **Algorithm**:
    * 1. Convert direction [uvtw] → vecs using m_to_o_dir, then normalize
    * 2. Convert plane (hkil) → vecm using cross products of plane points, then normalize
    * 3. Align vecm sign with ann[3] (c-component) for unidirectional modes like twinning
    * 
    * **Special cases**:
    * - Basal plane (0001): vecm = [0, 0, 1] directly
    * - Non-basal: Compute via cross product of two vectors in the plane
    * 
    * @param[in] ann Slip plane Miller indices (hkil) (length 4)
    *                - ann[0]=h, ann[1]=k, ann[2]=i, ann[3]=l
    *                - Must satisfy: h+k+i=0
    * @param[in] abb Slip direction Miller indices [uvtw] (length 4)
    *                - abb[0]=u, abb[1]=v, abb[2]=t, abb[3]=w
    *                - Must satisfy: u+v+t=0
    * @param[out] vecm Slip plane normal in Cartesian coordinates (length 3, unit vector)
    *                  - Perpendicular to slip plane
    *                  - Normalized: ‖vecm‖ = 1
    * @param[out] vecs Slip direction in Cartesian coordinates (length 3, unit vector)
    *                  - Tangent to slip plane
    *                  - Normalized: ‖vecs‖ = 1
    * @param[in] cOverA Lattice parameter ratio c/a
    *                   - Determines relative scaling of c-axis
    * 
    * @note Ensures vecm ⊥ vecs (orthogonality verified in DEBUG mode)
    * @note Both outputs are unit vectors
    * @note Sign of vecm is aligned with ann[3] for twinning compatibility
    * @note Fails with ECMECH_FAIL if sum constraints violated (DEBUG builds)
    * 
    * @see m_to_o_dir for Miller direction to Cartesian conversion
    * @see vecCrossProd for cross product computation
    * @see SlipGeomHCPaBRYcaY1 for usage in defining HCP slip systems
    */
   __ecmech_hdev__
   inline
   void
   miller_to_orthog_sngl(const double* const ann,
                         const double* const abb,
                         double* const vecm,
                         double* const vecs,
                         double cOverA
                         ) {
#ifndef NO_CHECKS
      // these are ecmech::nMiller, but sum over only first three entries on purpose
      if (fabs(vecsssum<ecmech::ndim>(ann)) > idp_eps_sqrt) {
         ECMECH_FAIL(__func__, "bad ann");
      }
      if (fabs(vecsssum<ecmech::ndim>(abb)) > idp_eps_sqrt) {
         ECMECH_FAIL(__func__, "bad abb");
      }
#endif

      // first do direction
      //
      m_to_o_dir(vecs, abb, cOverA);
      //
      vecsVNormalize<ecmech::ndim>(vecs);

      // now do plane;
      // fix me : this algorithm is clunky and not particularly efficient
      //
      if (vecsssumabs<ecmech::ndim>(ann) < idp_eps_sqrt) { // sum over only first three entries on purpose
         // basal plane
         vecm[0] = 0.; vecm[1] = 0.; vecm[2] = 1.;
      }
      else {
         double m_a[ecmech::nMiller] = { 0. };
         double m_b[ecmech::nMiller] = { 0. };
         //
         if (fabs(ann[0]) < idp_eps_sqrt) {
            // use second two axes for basal plane points
            m_a[1] = one / ann[1];
            m_b[2] = one / ann[2];
         }
         else if (fabs(ann[1]) < idp_eps_sqrt) {
            // use first and third axes for basal plane points
            m_a[0] = one / ann[0];
            m_b[2] = one / ann[2];
         }
         else {
            // use first and second axes for basal plane points
            m_a[0] = one / ann[0];
            m_b[1] = one / ann[1];
         }
         //
         double pnt_a[ecmech::ndim];
         m_to_o_dir(pnt_a, m_a, cOverA);
         double pnt_b[ecmech::ndim];
         m_to_o_dir(pnt_b, m_b, cOverA);
         double vec_basal[ecmech::ndim];
         for (int iN = 0; iN<ecmech::ndim; ++iN) {
            vec_basal[iN] = pnt_a[iN] - pnt_b[iN];
         }

         //
         double vec_nb[ecmech::ndim];
         if (fabs(ann[3]) < idp_eps_sqrt) {
            // normal is in basal plane
            vec_nb[0] = 0.; vec_nb[1] = 0.; vec_nb[2] = 1.;
            vecCrossProd(vecm, vec_basal, vec_nb);
         }
         else {
            // normal is not in basal plane
            double m_c[ecmech::nMiller] = { 0. };
            m_c[3] = one / ann[3];
            double pnt_c[ecmech::ndim];
            m_to_o_dir(pnt_c, m_c, cOverA);
            for (int iN = 0; iN<ecmech::ndim; ++iN) {
               vec_nb[iN] = pnt_c[iN] - pnt_b[iN];
            }

            vecCrossProd(vecm, vec_basal, vec_nb);
         }
      }
      //
      vecsVNormalize<ecmech::ndim>(vecm);
      //
      // align vecm with ann; this can be important for unidirectional modes like twinning
      if (vecm[2] * ann[3] < 0.) {
         for (int iN = 0; iN<ecmech::ndim; ++iN) {
            vecm[iN] = -vecm[iN];
         }
      }

#ifndef NO_CHECKS
      if (fabs(vecsyadotb<ecmech::ndim>(vecm, vecs)) > idp_eps_sqrt) {
         ECMECH_FAIL(__func__, "internal error");
      }
#endif
   } // miller_to_orthog_sngl

   /** @} */ // end of miller_indices group

   /**
    * @defgroup debug_utilities Debug and Printing Utilities
    * @brief Formatted output for vectors and matrices (host-only, debug builds).
    * 
    * These functions are only available in debug builds (__ecmech_host_only__)
    * and print formatted representations to std::cout.
    * 
    * @{
    */

#if defined(ECMECH_DEBUG) && defined(__ecmech_host_only__)

   /**
    * @brief Print vector to output stream with formatted precision.
    * 
    * Outputs vector components in row format with fixed precision and width.
    * Compile-time template version for known vector size.
    * 
    * **Format**: Each component printed as: width=21, precision=14 decimal places
    * 
    * @tparam n Vector length (compile-time constant)
    * @param[in] y Vector to print (length n)
    * @param[in,out] oss Output stream (default: std::cout)
    * 
    * @note Only available in debug builds (ECMECH_DEBUG && __ecmech_host_only__)
    * @note Outputs newline after all components
    */
   template<int n>
   inline void
   printVec(const double* const y, std::ostream & oss = std::cout) {
      for (int iX = 0; iX<n; ++iX) {
         oss << std::setw(21) << std::setprecision(14) << y[iX] << " ";
      }

      oss << std::endl;
   }

   /**
    * @brief Print vector to output stream with formatted precision (runtime size).
    * 
    * Outputs vector components in row format with fixed precision and width.
    * Runtime version for dynamic vector sizes.
    * 
    * **Format**: Each component printed as: width=21, precision=14 decimal places
    * 
    * @param[in] y Vector to print
    * @param[in] n Vector length (runtime parameter)
    * @param[in,out] oss Output stream
    * 
    * @note Only available in debug builds (ECMECH_DEBUG && __ecmech_host_only__)
    * @note No default stream parameter (must specify)
    * @note Outputs newline after all components
    */
   inline void
   printVec(const double* const y, int n, std::ostream & oss) {
      for (int iX = 0; iX<n; ++iX) {
         oss << std::setw(21) << std::setprecision(14) << y[iX] << " ";
      }

      oss << std::endl;
   }

   /**
    * @brief Print square matrix to output stream with formatted precision.
    * 
    * Outputs matrix in row-major format with each element formatted to
    * fixed width and precision. Each row printed on separate line.
    * Compile-time template version for square matrices.
    * 
    * **Format**: Each element printed as: width=21, precision=14 decimal places
    * 
    * @tparam n Matrix dimension (n×n, compile-time constant)
    * @param[in] A Matrix to print (n×n, row-major storage)
    * @param[in,out] oss Output stream (default: std::cout)
    * 
    * @note Only available in debug builds (ECMECH_DEBUG && __ecmech_host_only__)
    * @note Uses ECMECH_NN_INDX for row-major indexing
    * @note Outputs extra newline after matrix
    */
   template<int n>
   inline void
   printMat(const double* const A, std::ostream & oss = std::cout) {
      for (int iX = 0; iX<n; ++iX) {
         for (int jX = 0; jX<n; ++jX) {
            oss << std::setw(21) << std::setprecision(14) << A[ECMECH_NN_INDX(iX, jX, n)] << " ";
         }
         oss << std::endl;
      }
      oss << std::endl;
   }

   /**
    * @brief Print non-square matrix to output stream with formatted precision.
    * 
    * Outputs matrix in row-major format with each element formatted to
    * fixed width and precision. Each row printed on separate line.
    * Compile-time template version for rectangular matrices.
    * 
    * **Format**: Each element printed as: width=21, precision=14 decimal places
    * 
    * @tparam n Number of rows (compile-time constant)
    * @tparam m Number of columns (compile-time constant)
    * @param[in] A Matrix to print (n×m, row-major storage)
    * @param[in,out] oss Output stream (default: std::cout)
    * 
    * @note Only available in debug builds (ECMECH_DEBUG && __ecmech_host_only__)
    * @note Uses ECMECH_NM_INDX for row-major indexing
    * @note Outputs extra newline after matrix
    */
   template<int n, int m>
   inline void
   printMat(const double* const A, std::ostream & oss = std::cout) {
      for (int iX = 0; iX<n; ++iX) {
         for (int jX = 0; jX<m; ++jX) {
            oss << std::setw(21) << std::setprecision(14) << A[ECMECH_NM_INDX(iX, jX, n, m)] << " ";
         }

         oss << std::endl;
      }
      oss << std::endl;
   }

#endif
   /** @} */ // end of debug_utilities group
} // namespace ecmech

#endif // ECMECH_UTIL_H
