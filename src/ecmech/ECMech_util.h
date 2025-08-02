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

//
// maybe replace this macro with RAJA::View machinery at some point
//
// do row-major storage to be consistent with SNLS solver, and with RAJA
// row-major storage
#define ECMECH_NN_INDX(p, q, nDim) (p) * (nDim) + (q)
#define ECMECH_NM_INDX(p, q, pDim, qDim) (p) * (qDim) + (q)
//
// column-major storage would look like :
// define ECMECH_NN_INDX(p,q,nDim) (p)+(q)*(nDim)
// define ECMECH_NM_INDX(p,q,pDim,qDim) (p)+(q)*(pDim)

namespace ecmech {
   template<int n>
   __ecmech_hdev__
   inline void vecsVapb(double* const v,
                        const double* const a,
                        const double* const b) {
      for (int i = 0; i<n; ++i) {
         v[i] = a[i] + b[i];
      }
   }

   template<int n>
   __ecmech_hdev__
   inline void vecsVAdiagB(double* const v,
                           const double* const a,
                           const double* const b) {
      for (int i = 0; i<n; ++i) {
         v[i] = a[i] * b[i];
      }
   }

   // standard vector inner product
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

   template<int n>
   __ecmech_hdev__
   inline double vecsssumabs(const double* const a) {
      double s = 0.0;
      for (int i = 0; i<n; ++i) {
         s += fabs(a[i]);
      }

      return s;
   }

   template<int n>
   __ecmech_hdev__
   inline double vecsssum(const double* const a) {
      double s = 0.0;
      for (int i = 0; i<n; ++i) {
         s += a[i];
      }

      return s;
   }

   __ecmech_hdev__
   inline double vecsssumabs_n(const double* const a, int n) {
      double s = 0.0;
      for (int i = 0; i<n; ++i) {
         s += fabs(a[i]);
      }

      return s;
   }

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

   template<int n>
   __ecmech_hdev__
   inline void vecsVxa(double* const v,
                       double x,
                       const double* const a) {
      for (int i = 0; i<n; ++i) {
         v[i] = x * a[i];
      }
   }

   /*
    * scale vector in place
    */
   template<int n>
   __ecmech_hdev__
   inline void vecsVsa(double* const a,
                       double s) {
      for (int i = 0; i<n; ++i) {
         a[i] *= s;
      }
   }

   template<int n>
   __ecmech_hdev__
   inline void vecsVNormalize(double* const v){
      const double norm = vecNorm<n>(v);
      const double s = (fabs(norm) > idp_eps) ? 1.0 / norm : 1.0 / idp_eps;
      vecsVsa<n>(v, s);
   }
   // -fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all -fsanitize=float-divide-by-zero -fsanitize=float-cast-overflow -fno-omit-frame-pointer -fno-optimize-sibling-calls 

   /**
    * @brief matrix transposed times vector for square matrix
    *
    * NOTE : unlike SUBROUTINE matt_x_vec_5, output is first argument
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
    * @brief vector transposed times non-square matrix
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
    * @brief non-square matrix times vector
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
    * @brief square matrix times vector
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
    * @brief P = A . B^T where A and B are non-square and both n-by-q, so that P is n-by-n
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
    * @brief P = A . B^T where A (n-by-q) and B (m-by-q) are non-square, so that P is n-by-m
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
    * @brief P = A . B where A (n-by-q) and B (q-by-m) are non-square, so that P is n-by-m
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
    * @brief outer product P_ij = a_i b_j where a and b are n-vectors so that P is n-by-n
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

   /**
      ! w_j = \frac{1}{2} \varepsilon_{ijk} W_{ik}
    */
   // SUBROUTINE skew_to_veccp(veccp, W)
   __ecmech_hdev__
   inline void skewToVeccp(double* const veccp, // (WVEC)
                           const double* const W // (DIMS,DIMS)
                           )
   {
      veccp[0] = W[ECMECH_NN_INDX(2, 1, ecmech::ndim)];
      veccp[1] = W[ECMECH_NN_INDX(0, 2, ecmech::ndim)];
      veccp[2] = W[ECMECH_NN_INDX(1, 0, ecmech::ndim)];
   }

   __ecmech_hdev__
   inline double trace3(const double* const A // (DIMS,DIMS)
                        ) {
      double trace =
         A[ECMECH_NN_INDX(0, 0, ecmech::ndim)] +
         A[ECMECH_NN_INDX(1, 1, ecmech::ndim)] +
         A[ECMECH_NN_INDX(2, 2, ecmech::ndim)];

      return trace;
   }

   // trace_to_vecds_s(vecds_s, dkk)
   __ecmech_hdev__
   inline double traceToVecdsS(double dkk) {
      double vecds_s = sqr3i * dkk;
      return vecds_s;
   }

   //
   __ecmech_hdev__
   inline double vecd_Deff(const double* const vecd) {
      double retval = vecNorm<ecmech::ntvec>(vecd);
      retval = sqr2b3 * retval;
      return retval;
   }

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

   /*
    * see setup_vel_grad_svec_kk and svec_kk_to_vecds in Fortran
    *
    * svec_kk[6] is not accessed, so it need not be there
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
       ! symmetric (non-deviatoric) matrix to vecds representation
    */
   // SUBROUTINE symm_to_vecds(vecds, A)
   __ecmech_hdev__
   inline void symmToVecds(double* const vecds, // (SVEC)
                           const double* const A // (DIMS,DIMS)
                           )
   {
      symmToVecd(vecds, A);
      double Akk = trace3(A);
      vecds[iSvecS] = traceToVecdsS(Akk);
   }

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

   // SUBROUTINE vecds_to_symm(A, vecds)
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

   // SUBROUTINE svec_p_to_svec(a_svec, a_svec_p)
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

   __ecmech_hdev__
   inline void inv_to_quat(double* const quat,
                           const double* const inv) {
      double a = inv[0] * 0.5;
      quat[0] = cos(a);
      a = sin(a);
      vecsVxa<nwvec>(&(quat[1]), a, &(inv[1]));
   }

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

   // emap is simply the condensed 3 component version of an angle-axis vector
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
    * @brief calculate quaternion product q = a . b
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
    * @brief calculate quaternion inverse q^{-1}
    * q^{-1} = 1 / |q| [q0, -q1, -q2, -q3]
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
    * @brief calculate quaternion product q' = q^{-1}_1 . q_2
    */
   __ecmech_hdev__
   inline void quat_rel_rotation(double* const qprime,
                                 const double* const q1,
                                 const double* const q2) {
      double q1_inv[4] = {1.0, 0.0, 0.0, 0.0};
      quat_inverse(q1_inv, q1);
      quat_prod(qprime, q1_inv, q2);
   }

   __ecmech_hdev__
   inline void get_c_quat(double* const c_quat,
                          const double* const dr_quat,
                          const double* const cn_quat) {
      // Compute : c = c_n * dr
      quat_prod(c_quat, cn_quat, dr_quat);
   }

   /**
    * @brief after \cite{kri-etal-94a}
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
    * @brief rot_mat_vecd same as old rot_mat_symm in mat_tensor_mod BUT argument order here is switched to have a uniform convention that the output goes first

       !     construct 5x5 crystal rotation matrix:
       !                                 T
       !             [A_sm]=[C][A_lat][C] <=>  {A_sm} = [Q]{A_lat}
       !     with: {A}={()/sqr2,sqr3b2*(),sqr2*(),sqr2*(),sqr2*()}

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

   /**
    * matrix M such that
    *        v = M . a
    * where v = vecds(V), b = vecds(B), a = vecds(A)
    * V = A . B - B . A
    * and A, B are symmetric matrices (not necessarily deviatoric)
    *
    * as V is linear in A, this is the same as the operation itself --
    * operation of B on A
    *
    * if want _dB, just call with cmv6a and take negative of result;
    * this will give opeartion of A on B
    */
   __ecmech_hdev__
   inline void M35_d_AAoB_dA(double* const M35, // nwvec * ntvec
                             const double* const cmv6b // nsvec or ntvec -- cmv6b[iSvecS] not accessed
                             ) {
#include "util/vb_d_vars_set.h"
      // include "M36_d_AAoB_dA.f90"

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
      ! derivative of quaternion parameters with respect to exponential map
      ! parameters
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
      ! like d_rot_mat_symm_latop

       ! derivative of 5x5 rotation operation with respect to components of c
       !
       ! {vec_lat} = [Q]^T {vec_sm}
       ! dvdc is d({vec_lat})/d{C}
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
      ! like d_rot_mat_skew_latop

       ! derivative of 3x3 rotation operation with respect to components of C
       !
       ! {vec_lat} = [Q]^T {vec_sm}
       ! dvdc is d({vec_lat})/d{C}
       !
       ! dvdc(i,k,l) = \pfrac{(Qw_ji W_j)}{C_kl}
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
       ! modeled after eval_d_dxi in evptl_util_mod;
       ! for fully implicit only;
       ! everything in reference constituent frame
       !
       ! dDapp derivatives are all through lattice rotations, so just TVEC rows instead of SVEC -- trace of applied D does not change with rotation
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

   /**
    * pre-multiply by [qr5x5, 0; 0, 1] or its transpose;
    *
    * NOTE : NOT set up so that M_in and M_out may be the same
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

   /**
    * @brief covert Miller indices into a cartesian vector
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
    * convert slip systems (planes and directions) from miller indices to
    * vectors in an orthogonal coordinate system
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

#if defined(ECMECH_DEBUG) && defined(__ecmech_host_only__)

   template<int n>
   inline void
   printVec(const double* const y, std::ostream & oss = std::cout) {
      for (int iX = 0; iX<n; ++iX) {
         oss << std::setw(21) << std::setprecision(14) << y[iX] << " ";
      }

      oss << std::endl;
   }

   inline void
   printVec(const double* const y, int n, std::ostream & oss) {
      for (int iX = 0; iX<n; ++iX) {
         oss << std::setw(21) << std::setprecision(14) << y[iX] << " ";
      }

      oss << std::endl;
   }
   
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
} // namespace ecmech

#endif // ECMECH_UTIL_H
