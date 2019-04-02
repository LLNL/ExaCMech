// -*-c++-*-

#ifndef ECMECH_UTIL_H
#define ECMECH_UTIL_H

#include <cmath>

#ifdef DEBUG
#ifdef __cuda_host_only__
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#endif
#endif

#include "ECMech_const.h"
#include "RAJA/RAJA.hpp"

//
// maybe replace this macro with RAJA::View machinery at some point
//
// do row-major storage to be consistent with SNLS solver, and with RAJA
// row-major storage
#define ECMECH_NN_INDX(p,q,nDim) (p)*(nDim)+(q)
#define ECMECH_NM_INDX(p,q,pDim,qDim) (p)*(qDim)+(q)
//
// column-major storage would look like :
// define ECMECH_NN_INDX(p,q,nDim) (p)+(q)*(nDim)
// define ECMECH_NM_INDX(p,q,pDim,qDim) (p)+(q)*(pDim)

namespace ecmech {

template< int n >
inline void vecsVapb( real8* const v,
                      const real8* const a,
                      const real8* const b) {
   for (int i=0; i<n; ++i) { v[i] = a[i] + b[i]; }
}

template< int n >
inline void vecsVAdiagB( real8* const v,
                         const real8* const a,
                         const real8* const b) {
   for (int i=0; i<n; ++i) { v[i] = a[i] * b[i]; }
}

// standard vector inner product
template< int n >
inline real8 vecsyadotb( const real8* const a,
                         const real8* const b) {
   real8 y = 0.0 ;
   for (int i=0; i<n; ++i) { y+= a[i] * b[i]; }
   return y ;
}

template< int n >
inline real8 vecsssumabs( const real8* const a) {
   real8 s = 0.0 ;
   for (int i=0; i<n; ++i) { s+= fabs(a[i]); }
   return s ;
}

inline real8 vecsssumabs_n( const real8* const a, int n) {
   real8 s = 0.0 ;
   for (int i=0; i<n; ++i) { s+= fabs(a[i]); }
   return s ;
}

template< int n >
inline real8 vecNorm( const real8* const v ){
   real8 retval = 0.0 ;
   for (int i=0; i<n; ++i) { retval += v[i]*v[i]; }
   retval = sqrt(retval) ;
   return retval ;
}

template< int n >
inline void vecsVxa( real8* const v,
                     real8 x,
                     const real8* const a) {
   for (int i=0; i<n; ++i) { v[i] = x * a[i]; }
}

/*
 * scale vector in place
 */
template< int n >
inline void vecsVsa( real8* const a,
                     real8 s) {
   for (int i=0; i<n; ++i) { a[i] = s * a[i]; }
}

template< int n >
inline void vecsVNormalize( real8* const v ){
   real8 s = vecNorm<n>(v) ;
   vecsVsa<n>(v, s) ;
}

/**
 * @brief matrix transposed times vector for square matrix
 *
 * NOTE : unlike SUBROUTINE matt_x_vec_5, output is first argument
 */
template< int n >
inline void vecsVMTa( real8* const p,
                      const real8* const M,
                      const real8* const v) {
   for (int i=0; i<n; ++i) {
      p[i] = 0.0 ;
      for (int j=0; j<n; ++j) {
         p[i] += M[ECMECH_NN_INDX(j,i,n)] * v[j] ;
      }
   }
}

/**
 * @brief vector transposed times non-square matrix
 */
template< int n, int q >
inline void vecsVaTM( real8* const p,
                      const real8* const a,
                      const real8* const M) {
   for (int iQ=0; iQ<q; ++iQ) {
      p[iQ] = 0.0 ;
      for (int iN=0; iN<n; ++iN) {
         p[iQ] += a[iN] * M[ECMECH_NM_INDX(iN,iQ,n,q)] ;
      }
   }
}

/**
 * @brief non-square matrix times vector
 */
template< int n, int q >
inline void vecsVMa( real8* const p,
                      const real8* const M,
                      const real8* const a) {
   for (int iN=0; iN<n; ++iN) {
      p[iN] = 0.0 ;
      for (int iQ=0; iQ<q; ++iQ) {
         p[iN] += M[ECMECH_NM_INDX(iN,iQ,n,q)] * a[iQ] ;
      }
   }
}

/**
 * @brief square matrix times vector
 */
template< int n  >
inline void vecsVMa( real8* const p,
                      const real8* const M,
                      const real8* const a) {
   for (int iN=0; iN<n; ++iN) {
      p[iN] = 0.0 ;
      for (int jN=0; jN<n; ++jN) {
         p[iN] += M[ECMECH_NN_INDX(iN,jN,n)] * a[jN] ;
      }
   }
}

/**
 * @brief P = A . B^T where A and B are non-square and both n-by-q, so that P is n-by-n
 */
template< int n, int q >
inline void vecsMABT( real8* const P,
                      const real8* const A,
                      const real8* const B)
{
   
   for ( int ij=0; ij<n*n; ++ij ){
      P[ij] = 0.0 ;
   }
   
   for ( int jN = 0; jN < n; ++jN ) {
      for ( int iQ=0; iQ < q; ++iQ ) {
         real8 temp = B[ECMECH_NM_INDX(jN,iQ,n,q)] ;
         for ( int iN = 0; iN < n; ++iN ) {
            P[ECMECH_NM_INDX(iN,jN,n,n)] += A[ECMECH_NM_INDX(iN,iQ,n,q)] * temp ;
         }
      }
   }
   
}

/**
 * @brief P = A . B^T where A (n-by-q) and B (m-by-q) are non-square, so that P is n-by-m
 */
template< int n, int m, int q >
inline void vecsMABT( real8* const P,
                      const real8* const A,
                      const real8* const B) {
   
   for ( int ij=0; ij<n*m; ++ij ){
      P[ij] = 0.0 ;
   }
   
   for ( int jM = 0; jM < m; ++jM ) {
      for ( int iQ=0; iQ < q; ++iQ ) {
         real8 temp = B[ECMECH_NM_INDX(jM,iQ,m,q)] ;
         for ( int iN = 0; iN < n; ++iN ) {
            P[ECMECH_NM_INDX(iN,jM,n,m)] += A[ECMECH_NM_INDX(iN,iQ,n,q)] * temp ;
         }
      }
   }
   
}

/**
 * @brief P = A . B where A (n-by-q) and B (q-by-m) are non-square, so that P is n-by-m
 */
template< int n, int m, int q >
inline void vecsMAB( real8* const P,
                      const real8* const A,
                      const real8* const B) {
   
   for ( int ij=0; ij<n*m; ++ij ){
      P[ij] = 0.0 ;
   }
   
   for ( int iN = 0; iN < n; ++iN ) {
      for ( int jM = 0; jM < m; ++jM ) {
         for ( int iQ=0; iQ < q; ++iQ ) {
            P[ECMECH_NM_INDX(iN,jM,n,m)] += A[ECMECH_NM_INDX(iN,iQ,n,q)] * B[ECMECH_NM_INDX(iQ,jM,q,m)] ;
         }
      }
   }
   
}

/**
 * @brief outer product P_ij = a_i b_j where a and b are n-vectors so that P is n-by-n
 */
template< int n >
inline void vecsMaTb( real8* const P,
                      const real8* const a,
                      const real8* const b)
{
   for ( int iN = 0; iN < n; ++iN ) {
      for ( int jN = 0; jN < n; ++jN ) {
         P[ECMECH_NN_INDX(iN,jN,n)] = a[iN] * b[jN] ;
      }
   }
}

template< int n >
inline void vecsMsymm( real8* const P,
                       const real8* const A)
{
   for ( int iN = 0; iN < n; ++iN ) {
      for ( int jN = 0; jN < n; ++jN ) {
         P[ECMECH_NN_INDX(iN,jN,n)] = 0.5*( A[ECMECH_NN_INDX(iN,jN,n)] + A[ECMECH_NN_INDX(jN,iN,n)] );
      }
   }
}

template< int n >
inline void vecsMskew( real8* const Q,
                       const real8* const A)
{
   for ( int iN = 0; iN < n; ++iN ) {
      for ( int jN = 0; jN < n; ++jN ) {
         Q[ECMECH_NN_INDX(iN,jN,n)] = 0.5*( A[ECMECH_NN_INDX(iN,jN,n)] - A[ECMECH_NN_INDX(jN,iN,n)] );
      }
   }
}

/**
   ! w_j = \frac{1}{2} \varepsilon_{ijk} W_{ik}
 */
// SUBROUTINE skew_to_veccp(veccp, W)
inline void skewToVeccp( real8* const veccp, // (WVEC)
                         const real8* const W // (DIMS,DIMS)
                         )
{
    veccp[0] = W[ECMECH_NN_INDX(2,1,ecmech::ndim)] ;
    veccp[1] = W[ECMECH_NN_INDX(0,2,ecmech::ndim)] ;
    veccp[2] = W[ECMECH_NN_INDX(1,0,ecmech::ndim)] ;
}

inline real8 trace3( const real8* const A // (DIMS,DIMS)
                     ) {
   real8 trace =
      A[ECMECH_NN_INDX(0,0,ecmech::ndim)] +
      A[ECMECH_NN_INDX(1,1,ecmech::ndim)] +
      A[ECMECH_NN_INDX(2,2,ecmech::ndim)] ;
                   
   return trace ;
}

// trace_to_vecds_s(vecds_s, dkk)
inline real8 traceToVecdsS(real8 dkk) {
   real8 vecds_s = sqr3i * dkk ;
   return vecds_s ;
}

// 
inline real8 vecd_Deff(const real8* const vecd) {
   real8 retval = vecNorm<ecmech::ntvec>(vecd) ;
   retval = sqr2b3 * retval ;
   return retval ;
}

inline void symmToVecd( real8* const vecd, // (TVEC)
                        const real8* const A // (DIMS,DIMS)
                        )
{
    vecd[0] = sqr2i * (A[ECMECH_NN_INDX(0,0,ecmech::ndim)] - A[ECMECH_NN_INDX(1,1,ecmech::ndim)]) ;
    vecd[1] = sqr6i * (two * A[ECMECH_NN_INDX(2,2,ecmech::ndim)] - A[ECMECH_NN_INDX(0,0,ecmech::ndim)] - A[ECMECH_NN_INDX(1,1,ecmech::ndim)]) ; // = sqr6i * (3*A33 - Akk) = sqr3b2 * (A33 - Akk/3) = sqr3b2 * Adev33
    vecd[2] = sqr2 * A[ECMECH_NN_INDX(1,0,ecmech::ndim)] ;
    vecd[3] = sqr2 * A[ECMECH_NN_INDX(2,0,ecmech::ndim)] ;
    vecd[4] = sqr2 * A[ECMECH_NN_INDX(2,1,ecmech::ndim)] ;
}

/*
 * see setup_vel_grad_svec_kk and svec_kk_to_vecds in Fortran
 *
 * svec_kk[6] is not accessed, so it need not be there
 */
inline void svecToVecd( real8* const vecd, // (TVEC)
                        const real8* const svec_kk // (SVEC[+1])
                        )
{
   vecd[0] = sqr2i * (svec_kk[0] - svec_kk[1]) ;
   vecd[1] = sqr3b2 * svec_kk[2] ;
   vecd[2] = sqr2 * svec_kk[5] ;
   vecd[3] = sqr2 * svec_kk[4] ;
   vecd[4] = sqr2 * svec_kk[3] ;
   // vecds[5] = traceToVecdsS( svec_kk[6] ) ;
}

/**
    ! symmetric (non-deviatoric) matrix to vecds representation
 */
// SUBROUTINE symm_to_vecds(vecds, A)
inline void symmToVecds( real8* const vecds, // (SVEC)
                         const real8* const A // (DIMS,DIMS)
                         )
{
   symmToVecd(vecds, A) ;
   real8 Akk = trace3(A) ;
   vecds[iSvecS] = traceToVecdsS(Akk) ;
}       

inline real8 vecsInnerSvecDev( const real8* const stressSvec,
                               const real8* const dSvec ) {
   real8 retval = 
      stressSvec[0] * dSvec[0] +
      stressSvec[1] * dSvec[1] +
      stressSvec[2] * dSvec[2] +
      two * (
      stressSvec[3] * dSvec[3] +
      stressSvec[4] * dSvec[4] +
      stressSvec[5] * dSvec[5] ) ;
   return retval ;
}

// SUBROUTINE vecds_to_symm(A, vecds)
inline void vecdsToSvecP( real8* const svecp, // (SVEC+1)
                          const real8* const vecds // (SVEC)
                          )
{
   
   svecp[iSvecP] = -sqr3i * vecds[iSvecS] ; // -Akk_by_3
   
   real8 t1 = sqr2i * vecds[0] ;
   real8 t2 = sqr6i * vecds[1] ;
   //
   svecp[0] =    t1 - t2 ;        // 11'
   svecp[1] =   -t1 - t2 ;        // 22'
   svecp[2] = sqr2b3 * vecds[1] ; // 33'
   svecp[3] = sqr2i  * vecds[4] ; // 23
   svecp[4] = sqr2i  * vecds[3] ; // 31
   svecp[5] = sqr2i  * vecds[2] ; // 12

}       

inline void matToPQ( real8* const P_vecd,  // ntvec
                     real8* const Q_veccp, // nwvec
                     const real8* const T // ndim*ndim
                     )
{
   
   // CALL mat_to_symm_3(crys%p_ref(:,:,is), crys%t_ref(:,:,is))
   // CALL symm_to_vecds(P_ref_svec, crys%p_ref(:,:,is))
   // crys%P_ref_vec(:, is) = P_ref_svec(1:TVEC)
   //
   real8 P[ ecmech::ndim*ecmech::ndim ] ;
   vecsMsymm< ndim >(P,T) ;
   symmToVecd(P_vecd, P) ;
   
   // CALL mat_to_skew_3(crys%q_ref(:,:,is), crys%t_ref(:,:,is))
   // CALL skew_to_veccp(crys%q_ref_vec(:,is), crys%q_ref(:,:,is))
   real8 Q[ ecmech::ndim*ecmech::ndim ] ;
   vecsMskew< ndim >(Q,T) ;
   skewToVeccp(Q_veccp, Q) ;
   
}
   
inline void inv_to_quat(real8* const quat,
                        const real8* const inv) {
   real8 a = inv[0] * 0.5 ;
   quat[0] = cos(a) ;
   a = sin(a) ;
   vecsVxa<nwvec>(&(quat[1]), a, &(inv[1])) ;
}

inline void emap_to_quat(real8* const quat,
                         const real8* const emap) {
   real8 inv[invdim] = { 0.0, 1.0, 0.0, 0.0 };
   inv[0] = vecNorm<emapdim>(emap) ;
   if ( inv[0] > idp_tiny_sqrt ) {
      real8 invInv = 1.0 / inv[0] ;
      vecsVxa<emapdim>(&(inv[1]), invInv, emap) ;
   } // else, emap is effectively zero, so axis does not matter
   inv_to_quat(quat, inv);
}

/**
 * @brief calculate quaternion product q = a . b
 */
inline void quat_prod(real8* const q,
                       const real8* const a,
                       const real8* const b) {

   q[0] = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3] ;

   q[1] = a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2] ;
   q[2] = a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1] ;
   q[3] = a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0] ;
   
}

inline void get_c_quat(real8* const c_quat,
                       const real8* const dr_quat,
                       const real8* const cn_quat) {
   //  Compute : c = c_n * dr
   quat_prod(c_quat, cn_quat, dr_quat) ;
}

/**
 * @breif after \cite{kri-etal-94a}
 */
inline void quat_to_tensor(real8* const c, // ndim * ndim
                           const real8* const quat // qdim
                           ) {

   real8 x0sq = quat[0]*quat[0] ;
   real8 x1sq = quat[1]*quat[1] ;
   real8 x2sq = quat[2]*quat[2] ;
   real8 x3sq = quat[3]*quat[3] ;

   real8 x0x1 = quat[0]*quat[1] ;
   real8 x0x2 = quat[0]*quat[2] ;
   real8 x0x3 = quat[0]*quat[3] ;

   real8 x1x2 = quat[1]*quat[2] ;
   real8 x1x3 = quat[1]*quat[3] ;

   real8 x2x3 = quat[2]*quat[3] ; 

   c[ECMECH_NN_INDX(0,0,ndim)] = x0sq+x1sq-x2sq-x3sq ;
   c[ECMECH_NN_INDX(0,1,ndim)] = two*(x1x2-x0x3)     ;
   c[ECMECH_NN_INDX(0,2,ndim)] = two*(x1x3+x0x2)     ;
   c[ECMECH_NN_INDX(1,0,ndim)] = two*(x1x2+x0x3)     ;
   c[ECMECH_NN_INDX(1,1,ndim)] = x0sq-x1sq+x2sq-x3sq ;
   c[ECMECH_NN_INDX(1,2,ndim)] = two*(x2x3-x0x1)     ;
   c[ECMECH_NN_INDX(2,0,ndim)] = two*(x1x3-x0x2)     ;
   c[ECMECH_NN_INDX(2,1,ndim)] = two*(x2x3+x0x1)     ;
   c[ECMECH_NN_INDX(2,2,ndim)] = x0sq-x1sq-x2sq+x3sq ;

}

/**
 * @brief rot_mat_vecd same as old rot_mat_symm in mat_tensor_mod BUT argument order here is switched to have a uniform convention that the output goes first
 
    !     construct 5x5 crystal rotation matrix: 
    !                                 T
    !             [A_sm]=[C][A_lat][C] <=>  {A_sm} = [Q]{A_lat}
    !     with: {A}={()/sqr2,sqr3b2*(),sqr2*(),sqr2*(),sqr2*()}
    
 */
inline void get_rot_mat_vecd( real8* const qr5x5_raw, // ntvec * ntvec
                              const real8* const c // ndim * ndim
                              ) {

// include "mc_vars.f90"
// include "set_mc.f90"
#include "mc_vars_set.h"   


//    IF ((UBOUND(c,DIM=1) /= DIMS) .OR. (UBOUND(c,DIM=2) /= DIMS)) &
//         & CALL consider_ierr(1,location,CIERR_DIMS_p,IERR_FATAL_p)
//    IF ((UBOUND(qr5x5,DIM=1) /= TVEC) .OR. (UBOUND(qr5x5,DIM=2) /= TVEC)) &
//         & CALL consider_ierr(1,location,CIERR_DIMS_p,IERR_FATAL_p)

   RAJA::View< real8, RAJA::Layout<2> > qr5x5(qr5x5_raw, ecmech::ntvec, ecmech::ntvec) ;

//     ! if do not want to assume (c31**2+c32**2+c33**2=1)
//     qr5x5(1, 1)  =  c33 * c33 - onehalf * (c31 * c31 + c32 * c32)
   
    qr5x5(0, 0)  =  onehalf * (c11 * c11 - c12 * c12 - c21 * c21 + c22 * c22) ;
    qr5x5(0, 1)  =  sqr3 * onehalf * (c13 * c13 - c23 * c23) ;
    qr5x5(0, 2)  =  c11 * c12 - c21 * c22 ;
    qr5x5(0, 3)  =  c11 * c13 - c21 * c23 ;
    qr5x5(0, 4)  =  c12 * c13 - c22 * c23 ;
    qr5x5(1, 0)  =  sqr3 * onehalf * (c31 * c31 - c32 * c32) ;
    qr5x5(1, 1)  =  thrhalf * c33 * c33 - onehalf ;
    qr5x5(1, 2)  =  sqr3 * c31 * c32 ;
    qr5x5(1, 3)  =  sqr3 * c31 * c33 ;
    qr5x5(1, 4)  =  sqr3 * c32 * c33 ;
    qr5x5(2, 0)  =  c11 * c21 - c12 * c22 ;
    qr5x5(2, 1)  =  sqr3 * c13 * c23 ;
    qr5x5(2, 2)  =  c11 * c22 + c12 * c21 ;
    qr5x5(2, 3)  =  c11 * c23 + c13 * c21 ;
    qr5x5(2, 4)  =  c12 * c23 + c13 * c22 ;
    qr5x5(3, 0)  =  c11 * c31 - c12 * c32 ;
    qr5x5(3, 1)  =  sqr3 * c13 * c33 ;
    qr5x5(3, 2)  =  c11 * c32 + c12 * c31 ;
    qr5x5(3, 3)  =  c11 * c33 + c13 * c31 ;
    qr5x5(3, 4)  =  c12 * c33 + c13 * c32 ;
    qr5x5(4, 0)  =  c21 * c31 - c22 * c32 ;
    qr5x5(4, 1)  =  sqr3 * c23 * c33 ;
    qr5x5(4, 2)  =  c21 * c32 + c22 * c31 ;
    qr5x5(4, 3)  =  c21 * c33 + c23 * c31 ;
    qr5x5(4, 4)  =  c22 * c33 + c23 * c32 ;

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
inline void M35_d_AAoB_dA( real8* const M35, // nwvec * ntvec
                       const real8* const cmv6b // nsvec or ntvec -- cmv6b[iSvecS] not accessed
                       ) {

#include "vb_d_vars_set.h"
// include "M36_d_AAoB_dA.f90"
   
   M35[ECMECH_NM_INDX(0,0,nwvec,ntvec)] = vb5*onehalf ;
   M35[ECMECH_NM_INDX(1,0,nwvec,ntvec)] = vb4*onehalf ;
   M35[ECMECH_NM_INDX(2,0,nwvec,ntvec)] = -vb3 ;
   M35[ECMECH_NM_INDX(0,1,nwvec,ntvec)] = vb5*halfsqr3 ;
   M35[ECMECH_NM_INDX(1,1,nwvec,ntvec)] = -vb4*halfsqr3 ;
   M35[ECMECH_NM_INDX(2,1,nwvec,ntvec)] = zero ;
   M35[ECMECH_NM_INDX(0,2,nwvec,ntvec)] = -vb4*onehalf ;
   M35[ECMECH_NM_INDX(1,2,nwvec,ntvec)] = vb5*onehalf ;
   M35[ECMECH_NM_INDX(2,2,nwvec,ntvec)] = vb1 ;
   M35[ECMECH_NM_INDX(0,3,nwvec,ntvec)] = vb3*onehalf ;
   M35[ECMECH_NM_INDX(1,3,nwvec,ntvec)] = halfsqr3*vb2-vb1*onehalf ;
   M35[ECMECH_NM_INDX(2,3,nwvec,ntvec)] = -vb5*onehalf ;
   M35[ECMECH_NM_INDX(0,4,nwvec,ntvec)] = -vb1*onehalf-halfsqr3*vb2 ;
   M35[ECMECH_NM_INDX(1,4,nwvec,ntvec)] = -vb3*onehalf ;
   M35[ECMECH_NM_INDX(2,4,nwvec,ntvec)] = vb4*onehalf ;
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
inline void dquat_demap_T( real8* const dqdeT_raw, // (EMAPDIM_p,QDIM_p)
                           const real8* const emap // (EMAPDIM_p)
                           )
{

   const real8 theta_sm_a = 1e-9 ;
   const real8 oo48 = 1.0/48.0 ;


   real8 theta = vecNorm< emapdim >( emap ) ;

   real8 theta_inv, sthhbyth, halfsthh, na, nb, nc ;
   if ( abs(theta) < theta_sm_a ) {
      sthhbyth = onehalf - theta*theta * oo48 ; // truncated Taylor seriers; probably safe to just use onehalf and be done with it
      halfsthh = theta * oneqrtr ; // truncated Taylor seriers
      if ( abs(theta) < idp_tiny_sqrt ) {
         // n is arbitrary, as theta is effectively zero
         na = one; nb = zero; nc = zero;
      }
      else {
         theta_inv = one/theta ;
         na = emap[0]*theta_inv; nb = emap[1]*theta_inv; nc = emap[2]*theta_inv;
      }
   }
   else {
      halfsthh = sin(theta*onehalf) ;
      sthhbyth = halfsthh/theta ;
      halfsthh = halfsthh * onehalf ;
      theta_inv = one/theta ;
      na = emap[0]*theta_inv; nb = emap[1]*theta_inv; nc = emap[2]*theta_inv;
   }
   //
   real8 halfcthh = cos(theta*onehalf)*onehalf ;
   //
   // now have: halfsthh, sthhbyth, halfcthh, theta, na, nb, nc

   RAJA::View< real8, RAJA::Layout<2> > dqdeT(dqdeT_raw, ecmech::emapdim, ecmech::qdim) ;
   
   dqdeT(0,0) = -halfsthh * na ;
   dqdeT(1,0) = -halfsthh * nb ;
   dqdeT(2,0) = -halfsthh * nc ;

   real8 temp = na*na ;
   dqdeT(0,1) = halfcthh * temp + sthhbyth * (one - temp) ;
   //
   temp = nb*nb ;
   dqdeT(1,2) = halfcthh * temp + sthhbyth * (one - temp) ;
   //
   temp = nc*nc ;
   dqdeT(2,3) = halfcthh * temp + sthhbyth * (one - temp) ;

   temp = halfcthh - sthhbyth ;
   //
   real8 tempb = temp * na*nb ;
   dqdeT(1,1) = tempb ;
   dqdeT(0,2) = tempb ;
   //
   tempb = temp * na*nc ;
   dqdeT(2,1) = tempb ;
   dqdeT(0,3) = tempb ;
   //
   tempb = temp * nb*nc ;
   dqdeT(2,2) = tempb ;
   dqdeT(1,3) = tempb ;

}
   
inline void d_quat_to_tensor(real8* const dcdq_raw, // (DIMS,DIMS,QDIM_p)
                             const real8* const quat // (QDIM_p)
                             )
{
    real8 tqa = two * quat[0] ;
    real8 tqb = two * quat[1] ;
    real8 tqc = two * quat[2] ;
    real8 tqd = two * quat[3] ;

    RAJA::View< real8, RAJA::Layout<3> > dcdq(dcdq_raw, ecmech::ndim, ecmech::ndim, ecmech::qdim) ;
    
    // c(1,1) = x1sq+x2sq-x3sq-x4sq
    dcdq(0,0,0) =  tqa ;
    dcdq(0,0,1) =  tqb ;
    dcdq(0,0,2) = -tqc ;
    dcdq(0,0,3) = -tqd ;

    // c(0,1) = two*(x1x2-x0x3)
    dcdq(0,1,0) = -tqd ;
    dcdq(0,1,1) =  tqc ;
    dcdq(0,1,2) =  tqb ;
    dcdq(0,1,3) = -tqa ;

    // c(0,2) = two*(x1x3+x0x2)
    dcdq(0,2,0) =  tqc ;
    dcdq(0,2,1) =  tqd ;
    dcdq(0,2,2) =  tqa ;
    dcdq(0,2,3) =  tqb ;

    // c(1,0) = two*(x1x2+x0x3)
    dcdq(1,0,0) =  tqd ;
    dcdq(1,0,1) =  tqc ;
    dcdq(1,0,2) =  tqb ;
    dcdq(1,0,3) =  tqa ;

    // c(1,1) = x0sq-x1sq+x2sq-x3sq
    dcdq(1,1,0) =  tqa ;
    dcdq(1,1,1) = -tqb ;
    dcdq(1,1,2) =  tqc ;
    dcdq(1,1,3) = -tqd ;

    // c(1,2) = two*(x2x3-x0x1)
    dcdq(1,2,0) = -tqb ;
    dcdq(1,2,1) = -tqa ;
    dcdq(1,2,2) =  tqd ;
    dcdq(1,2,3) =  tqc ;

    // c(2,0) = two*(x1x3-x0x2)
    dcdq(2,0,0) = -tqc ;
    dcdq(2,0,1) =  tqd ;
    dcdq(2,0,2) = -tqa ;
    dcdq(2,0,3) =  tqb ;

    // c(2,1) = two*(x2x3+x0x1)
    dcdq(2,1,0) =  tqb ;
    dcdq(2,1,1) =  tqa ;
    dcdq(2,1,2) =  tqd ;
    dcdq(2,1,3) =  tqc ;

    // c(2,2) = x0sq-x1sq-x2sq+x3sq
    dcdq(2,2,0) =  tqa ;
    dcdq(2,2,1) = -tqb ;
    dcdq(2,2,2) = -tqc ;
    dcdq(2,2,3) =  tqd ;

}

/**
   ! like d_rot_mat_symm_latop

    ! derivative of 5x5 rotation operation with respect to components of c
    !
    ! {vec_lat} = [Q]^T {vec_sm}
    ! dvdc is d({vec_lat})/d{C}
*/
inline void d_rot_mat_vecd_latop( real8* const dvdc_raw, // (TVEC,DIMS,DIMS)
                                  const real8* const c, // (DIMS, DIMS)
                                  const real8* const vec_sm // (TVEC)
                                  )
{

#include "mc_vars_set.h"
#include "vad_vars_set.h"   

// include "d_Alat_dC.f90"
   RAJA::View< real8, RAJA::Layout<3> > dvdc(dvdc_raw, ecmech::ntvec, ecmech::ndim, ecmech::ndim) ;
   dvdc(0,0,0)   = c11*va1-c11*sqr3*va2*onethird+c21*va3+c31*va4 ;
   dvdc(1,0,0)   = -sqr3*(3*c11*va1-c11*sqr3*va2+3*c21*va3+3*c31*va4)*oneninth ;
   dvdc(2,0,0)   = c12*va1-c12*sqr3*va2*onethird+c22*va3+c32*va4 ;
   dvdc(3,0,0)   = c13*va1-c13*sqr3*va2*onethird+c23*va3+c33*va4 ;
   dvdc(4,0,0)   = zero ;
   dvdc(0,1,0)   = c11*va3-c21*va1-c21*sqr3*va2*onethird+c31*va5 ;
   dvdc(1,1,0)   = sqr3*(-3*c11*va3+3*c21*va1+c21*sqr3*va2-3*c31*va5)*oneninth ;
   dvdc(2,1,0)   = c12*va3-c22*va1-c22*sqr3*va2*onethird+c32*va5 ;
   dvdc(3,1,0)   = c13*va3-c23*va1-c23*sqr3*va2*onethird+c33*va5 ;
   dvdc(4,1,0)   = zero ;
   dvdc(0,2,0)   = c11*va4+c21*va5+twothird*c31*sqr3*va2 ;
   dvdc(1,2,0)   = -sqr3*c11*va4*onethird-sqr3*c21*va5*onethird-twothird*c31*va2 ;
   dvdc(2,2,0)   = c12*va4+c22*va5+twothird*c32*sqr3*va2 ;
   dvdc(3,2,0)   = c13*va4+c23*va5+twothird*c33*sqr3*va2 ;
   dvdc(4,2,0)   = zero ;
   dvdc(0,0,1)   = -c12*va1+c12*sqr3*va2*onethird-c22*va3-c32*va4 ;
   dvdc(1,0,1)   = sqr3*(-3*c12*va1+c12*sqr3*va2-3*c22*va3-3*c32*va4)*oneninth ;
   dvdc(2,0,1)   = c11*va1-c11*sqr3*va2*onethird+c21*va3+c31*va4 ;
   dvdc(3,0,1)   = zero ;
   dvdc(4,0,1)   = c13*va1-c13*sqr3*va2*onethird+c23*va3+c33*va4 ;
   dvdc(0,1,1)   = -c12*va3+c22*va1+c22*sqr3*va2*onethird-c32*va5 ;
   dvdc(1,1,1)   = sqr3*(-3*c12*va3+3*c22*va1+c22*sqr3*va2-3*c32*va5)*oneninth ;
   dvdc(2,1,1)   = c11*va3-c21*va1-c21*sqr3*va2*onethird+c31*va5 ;
   dvdc(3,1,1)   = zero ;
   dvdc(4,1,1)   = c13*va3-c23*va1-c23*sqr3*va2*onethird+c33*va5 ;
   dvdc(0,2,1)   = -c12*va4-c22*va5-twothird*c32*sqr3*va2 ;
   dvdc(1,2,1)   = -sqr3*c12*va4*onethird-sqr3*c22*va5*onethird-twothird*c32*va2 ;
   dvdc(2,2,1)   = c11*va4+c21*va5+twothird*c31*sqr3*va2 ;
   dvdc(3,2,1)   = zero ;
   dvdc(4,2,1)   = c13*va4+c23*va5+twothird*c33*sqr3*va2 ;
   dvdc(0,0,2)   = zero ;
   dvdc(1,0,2)   = -twothird*sqr3i*(-3*c13*va1+sqr3*c13*va2-3*c23*va3-3*c33*va4) ;
   dvdc(2,0,2)   = zero ;
   dvdc(3,0,2)   = c11*va1-c11*sqr3*va2*onethird+c21*va3+c31*va4 ;
   dvdc(4,0,2)   = c12*va1-c12*sqr3*va2*onethird+c22*va3+c32*va4 ;
   dvdc(0,1,2)   = zero ;
   dvdc(1,1,2)   = -twothird*sqr3i*(-3*va3*c13+3*va1*c23+sqr3*c23*va2-3*va5*c33) ;
   dvdc(2,1,2)   = zero ;
   dvdc(3,1,2)   = c11*va3-c21*va1-c21*sqr3*va2*onethird+c31*va5 ;
   dvdc(4,1,2)   = c12*va3-c22*va1-c22*sqr3*va2*onethird+c32*va5 ;
   dvdc(0,2,2)   = zero ;
   dvdc(1,2,2)   = twothird*sqr3*c13*va4+twothird*sqr3*c23*va5+fourthirds*c33*va2 ;
   dvdc(2,2,2)   = zero ;
   dvdc(3,2,2)   = c11*va4+c21*va5+twothird*c31*sqr3*va2 ;
   dvdc(4,2,2)   = c12*va4+c22*va5+twothird*c32*sqr3*va2 ;
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
inline void d_rot_mat_wveccp_latop( real8* const dvdc_raw, // (WVEC,DIMS,DIMS)
                                    // const real8* const c, // (DIMS, DIMS) // not used
                                    const real8* const cmv3w // (WVEC) // vec_sm(WVEC)
                                    )
{

#include "vw_vars_set.h"

// include "d_Wlat_dC.f90"
   RAJA::View< real8, RAJA::Layout<3> > dvdc(dvdc_raw, ecmech::nwvec, ecmech::ndim, ecmech::ndim) ;
   
   dvdc(0,0,0)   = vw1 ;
   dvdc(1,0,0)   = zero ;
   dvdc(2,0,0)   = zero ;
   dvdc(0,1,0)   = vw2 ;
   dvdc(1,1,0)   = zero ;
   dvdc(2,1,0)   = zero ;
   dvdc(0,2,0)   = vw3 ;
   dvdc(1,2,0)   = zero ;
   dvdc(2,2,0)   = zero ;
   dvdc(0,0,1)   = zero ;
   dvdc(1,0,1)   = vw1 ;
   dvdc(2,0,1)   = zero ;
   dvdc(0,1,1)   = zero ;
   dvdc(1,1,1)   = vw2 ;
   dvdc(2,1,1)   = zero ;
   dvdc(0,2,1)   = zero ;
   dvdc(1,2,1)   = vw3 ;
   dvdc(2,2,1)   = zero ;
   dvdc(0,0,2)   = zero ;
   dvdc(1,0,2)   = zero ;
   dvdc(2,0,2)   = vw1 ;
   dvdc(0,1,2)   = zero ;
   dvdc(1,1,2)   = zero ;
   dvdc(2,1,2)   = vw2 ;
   dvdc(0,2,2)   = zero ;
   dvdc(1,2,2)   = zero ;
   dvdc(2,2,2)   = vw3 ;

} // d_rot_mat_wveccp_latop

/**
    ! modeled after eval_d_dxi in evptl_util_mod;
    ! for fully implicit only;
    ! everything in reference constituent frame
    ! 
    ! dDapp derivatives are all through lattice rotations, so just TVEC rows instead of SVEC -- trace of applied D does not change with rotation
 */
inline void eval_d_dxi_impl_quat( real8* const dC_quat_dxi_T, // (WVEC,QDIM_p)
                                  // real8* const dC_matx_dxi, // (DIMS,DIMS,WVEC)
                                  real8* const dDapp_dxi,   // dDapp_dxi(TVEC, WVEC)
                                  real8* const dWapp_dxi,   // dWapp_dxi(WVEC, WVEC)
                                  const real8* const d_vecd_sm, // (TVEC), or (SVEC) is fine too
                                  const real8* const w_vec_sm, // (WVEC)
                                  const real8* const xi, // (WVEC)
                                  const real8* const Cn_quat, // (QDIM_p)
                                  const real8* const C_matx, // (DIMS,DIMS)
                                  const real8* const C_quat // (QDIM_p)
                                  // const real8* const A_quat // (QDIM_p) // not used
                                  ) {

   // working with quats, so do not call eval_d_cA_dxi(dc_dxi, dA_dxi, xi, c_n)
   //
   {
      real8 dA_quat_dxi_T[ ecmech::ndim * ecmech::qdim ] ; // (QDIM_p,DIMS)^T
      dquat_demap_T(dA_quat_dxi_T, xi) ;

      // can get away with these three calls as quat_prod is bilinear in the input arguments
      //
      quat_prod( &(dC_quat_dxi_T[ecmech::qdim*0]), Cn_quat, &(dA_quat_dxi_T[ecmech::qdim*0]) ) ;
      quat_prod( &(dC_quat_dxi_T[ecmech::qdim*1]), Cn_quat, &(dA_quat_dxi_T[ecmech::qdim*1]) ) ;
      quat_prod( &(dC_quat_dxi_T[ecmech::qdim*2]), Cn_quat, &(dA_quat_dxi_T[ecmech::qdim*2]) ) ;
   }
   // now have dC_quat_dxi

   real8 dC_matx_dxi[ (ecmech::ndim * ecmech::ndim) * ecmech::nwvec ] ; // (DIMS,DIMS,WVEC)
   {
      real8 dCmatx_dq[ (ecmech::ndim * ecmech::ndim) * ecmech::qdim ] ; // (DIMS,DIMS,QDIM_p)
      // get dC_matx_dxi
      d_quat_to_tensor(dCmatx_dq, C_quat) ;
      vecsMABT< ndim*ndim, nwvec, qdim >( dC_matx_dxi, dCmatx_dq, dC_quat_dxi_T ) ; // vecsMABT because _T on dC_quat_dxi_T
   }
    
   {
      real8 dD_dC_matx[ ecmech::ntvec * (ecmech::ndim * ecmech::ndim) ] ;
      d_rot_mat_vecd_latop(dD_dC_matx, C_matx, d_vecd_sm) ;
      //
      vecsMAB< ntvec, nwvec, ndim*ndim >( dDapp_dxi, dD_dC_matx, dC_matx_dxi ) ;
      // dDapp_dxi(SVEC,:) = zero
   }

   {
      real8 dW_dC_matx[ ecmech::nwvec * (ecmech::ndim * ecmech::ndim) ] ; // (WVEC,DIMS,DIMS)
      d_rot_mat_wveccp_latop(dW_dC_matx, // C_matx,
                             w_vec_sm) ;
      //
      vecsMAB< nwvec, nwvec, ndim*ndim >( dWapp_dxi, dW_dC_matx, dC_matx_dxi ) ;
   }

}

inline void
d_rot_mat_vecd_smop( real8* const dvdc_raw, // (TVEC,DIMS,DIMS)
                     const real8* const c, // (DIMS, DIMS)
                     const real8* const vec_lat // (TVEC)
                     ) 
{
    //
    // derivative of 5x5 rotation operation with respect to components of c
    //
    // {vec_sm} = [Q] {vec_lat}
    // dvdc is d({vec_sm})/d{C}
    //

#include "mc_vars_set.h"
#include "vadl_vars_set.h"

   // include "d_Asm_dC.f90"
   RAJA::View< real8, RAJA::Layout<3> > dvdc(dvdc_raw, ecmech::ntvec, ecmech::ndim, ecmech::ndim) ;
   dvdc(0,0,0) = c11*va1-c11*sqr3*va2*onethird+c12*va3+c13*va4 ;
   dvdc(1,0,0) = sqr3*(-3*c11*va1+sqr3*c11*va2-3*c12*va3-3*c13*va4)*oneninth ;
   dvdc(2,0,0) = c21*va1-c21*sqr3*va2*onethird+c22*va3+c23*va4 ;
   dvdc(3,0,0) = c31*va1-c31*sqr3*va2*onethird+c32*va3+c33*va4 ;
   dvdc(4,0,0) = zero ;
   dvdc(0,1,0) = -c21*va1+c21*sqr3*va2*onethird-c22*va3-c23*va4 ;
   dvdc(1,1,0) = -sqr3*(3*c21*va1-c21*sqr3*va2+3*c22*va3+3*c23*va4)*oneninth ;
   dvdc(2,1,0) = c11*va1-c11*sqr3*va2*onethird+c12*va3+c13*va4 ;
   dvdc(3,1,0) = zero ;
   dvdc(4,1,0) = c31*va1-c31*sqr3*va2*onethird+c32*va3+c33*va4 ;
   dvdc(0,2,0) = zero ;
   dvdc(1,2,0) = -twothird*sqr3i*(-3*c31*va1+c31*sqr3*va2-3*c32*va3-3*c33*va4) ;
   dvdc(2,2,0) = zero ;
   dvdc(3,2,0) = c11*va1-c11*sqr3*va2*onethird+c12*va3+c13*va4 ;
   dvdc(4,2,0) = c21*va1-c21*sqr3*va2*onethird+c22*va3+c23*va4 ;
   dvdc(0,0,1) = c11*va3-c12*va1-c12*sqr3*va2*onethird+c13*va5 ;
   dvdc(1,0,1) = sqr3*(-3*c11*va3+3*c12*va1+c12*sqr3*va2-3*c13*va5)*oneninth ;
   dvdc(2,0,1) = c21*va3-c22*va1-c22*sqr3*va2*onethird+c23*va5 ;
   dvdc(3,0,1) = c31*va3-c32*va1-c32*sqr3*va2*onethird+c33*va5 ;
   dvdc(4,0,1) = zero ;
   dvdc(0,1,1) = -c21*va3+c22*va1+c22*sqr3*va2*onethird-c23*va5 ;
   dvdc(1,1,1) = sqr3*(-3*c21*va3+3*c22*va1+sqr3*va2*c22-3*c23*va5)*oneninth ;
   dvdc(2,1,1) = c11*va3-c12*va1-c12*sqr3*va2*onethird+c13*va5 ;
   dvdc(3,1,1) = zero ;
   dvdc(4,1,1) = c31*va3-c32*va1-c32*sqr3*va2*onethird+c33*va5 ;
   dvdc(0,2,1) = zero ;
   dvdc(1,2,1) = -twothird*sqr3i*(-3*c31*va3+3*c32*va1+sqr3*va2*c32-3*c33*va5) ;
   dvdc(2,2,1) = zero ;
   dvdc(3,2,1) = c11*va3-c12*va1-c12*sqr3*va2*onethird+c13*va5 ;
   dvdc(4,2,1) = c21*va3-c22*va1-c22*sqr3*va2*onethird+c23*va5 ;
   dvdc(0,0,2) = c11*va4+c12*va5+twothird*c13*sqr3*va2 ;
   dvdc(1,0,2) = -sqr3*c11*va4*onethird-sqr3*c12*va5*onethird-twothird*c13*va2 ;
   dvdc(2,0,2) = c21*va4+c22*va5+twothird*c23*sqr3*va2 ;
   dvdc(3,0,2) = c31*va4+c32*va5+twothird*c33*sqr3*va2 ;
   dvdc(4,0,2) = zero ;
   dvdc(0,1,2) = -c21*va4-c22*va5-twothird*c23*sqr3*va2 ;
   dvdc(1,1,2) = -sqr3*c21*va4*onethird-sqr3*c22*va5*onethird-twothird*c23*va2 ;
   dvdc(2,1,2) = c11*va4+c12*va5+twothird*c13*sqr3*va2 ;
   dvdc(3,1,2) = zero ;
   dvdc(4,1,2) = c31*va4+c32*va5+twothird*c33*sqr3*va2 ;
   dvdc(0,2,2) = zero ;
   dvdc(1,2,2) = twothird*sqr3*c31*va4+twothird*sqr3*c32*va5+fourthirds*c33*va2 ;
   dvdc(2,2,2) = zero ;
   dvdc(3,2,2) = c11*va4+c12*va5+twothird*c13*sqr3*va2 ;
   dvdc(4,2,2) = c21*va4+c22*va5+twothird*c23*sqr3*va2 ;

} // d_rot_mat_vecd_smop

/**
 * pre-multiply by [qr5x5, 0; 0, 1] or its transpose;
 *
 * NOTE : NOT set up so that M_in and M_out may be the same 
 */
template< int n, bool l_T >
inline void
qr6x6_pre_mul( real8* const M_out, // 6xn
               const real8* const M_in, // 6xn
               const real8* const qr5x5 // 5x5
               )
{

   for (int iM=0; iM<ecmech::ntvec; ++iM ) { // only up to ntvec on purpose !
      for (int jM=0; jM<n; ++jM) {
         int ijM = ECMECH_NM_INDX(iM,jM,ecmech::nsvec,n) ;
         M_out[ijM] = 0.0 ;
         for (int pTvec=0; pTvec<ecmech::ntvec; ++pTvec) { 
            if (l_T) {
               M_out[ijM] += qr5x5[ECMECH_NN_INDX(pTvec,iM,ecmech::ntvec)] * M_in[ECMECH_NM_INDX(pTvec,jM,ecmech::nsvec,n)] ;
            }
            else {
               M_out[ijM] += qr5x5[ECMECH_NN_INDX(iM,pTvec,ecmech::ntvec)] * M_in[ECMECH_NM_INDX(pTvec,jM,ecmech::nsvec,n)] ;
            }
         }
      }
   }
   for (int jM=0; jM<n; ++jM) {
      int ijM = ECMECH_NM_INDX(iSvecS,jM,ecmech::nsvec,n) ;
      M_out[ijM] = M_in[ijM] ;
   }
   
} // qr6x6_pre_mul

template< bool l_ddsdde_gamma >
inline void
mtan_conv_sd_svec(real8* const mtanSD_raw,
                  const real8* const mtanSD_vecds_raw) {
   real8 C_raw[ecmech::nsvec2] ;
   real8 t1_vec[ecmech::nsvec], t2_vec[ecmech::nsvec], t3_vec[ecmech::nsvec] ;

   RAJA::View< real8, RAJA::Layout<2> > mtanSD(mtanSD_raw, ecmech::nsvec, ecmech::nsvec) ;   
   RAJA::View< real8 const, RAJA::Layout<2> > mtanSD_vecds(mtanSD_vecds_raw, ecmech::nsvec, ecmech::nsvec) ;   
   RAJA::View< real8, RAJA::Layout<2> > C(C_raw, ecmech::nsvec, ecmech::nsvec) ;   

   // mtanSD = T . mtanSD_vecds . T^{-1}

   // C = T . mtanSD_vecds
   // C(i,:) = T(i,k) . mtanSD_vecds(k,:) -- sum over k
   //
   for ( int jSvec=0; jSvec<ecmech::nsvec; ++jSvec ) { t3_vec[jSvec] = sqr3i * mtanSD_vecds(iSvecS,jSvec) ; }
   for ( int jSvec=0; jSvec<ecmech::nsvec; ++jSvec ) { t1_vec[jSvec] = sqr2i * mtanSD_vecds(     0,jSvec) ; }
   for ( int jSvec=0; jSvec<ecmech::nsvec; ++jSvec ) { t2_vec[jSvec] = sqr6i * mtanSD_vecds(     1,jSvec) ; }
   //
   for ( int jSvec=0; jSvec<ecmech::nsvec; ++jSvec ) { C(0,jSvec) =    t1_vec[jSvec] - t2_vec[jSvec] + t3_vec[jSvec] ; }
   for ( int jSvec=0; jSvec<ecmech::nsvec; ++jSvec ) { C(1,jSvec) =   -t1_vec[jSvec] - t2_vec[jSvec] + t3_vec[jSvec] ; }
   for ( int jSvec=0; jSvec<ecmech::nsvec; ++jSvec ) { C(2,jSvec) =   sqr2b3 * mtanSD_vecds(1,jSvec) + t3_vec[jSvec] ; }
   for ( int jSvec=0; jSvec<ecmech::nsvec; ++jSvec ) { C(3,jSvec) =    sqr2i * mtanSD_vecds(4,jSvec) ; }
   for ( int jSvec=0; jSvec<ecmech::nsvec; ++jSvec ) { C(4,jSvec) =    sqr2i * mtanSD_vecds(3,jSvec) ; }
   for ( int jSvec=0; jSvec<ecmech::nsvec; ++jSvec ) { C(5,jSvec) =    sqr2i * mtanSD_vecds(2,jSvec) ; }

   // mtanSD = C . T^{-1}
   // mtanSD(:,j) = C(:,k) . [T^{-1}](k,j) -- sum over k
   //
   for ( int jSvec=0; jSvec<ecmech::nsvec; ++jSvec ) { t3_vec[jSvec] = C(jSvec,iSvecS) * sqr3i ; }
   for ( int jSvec=0; jSvec<ecmech::nsvec; ++jSvec ) { t1_vec[jSvec] = C(jSvec,     0) * sqr2i ; }
   for ( int jSvec=0; jSvec<ecmech::nsvec; ++jSvec ) { t2_vec[jSvec] = C(jSvec,     1) * sqr6i ; }
   //
   for ( int jSvec=0; jSvec<ecmech::nsvec; ++jSvec ) { mtanSD(jSvec,0) =  t1_vec[jSvec] - t2_vec[jSvec] + t3_vec[jSvec] ; }
   for ( int jSvec=0; jSvec<ecmech::nsvec; ++jSvec ) { mtanSD(jSvec,1) = -t1_vec[jSvec] - t2_vec[jSvec] + t3_vec[jSvec] ; }
   for ( int jSvec=0; jSvec<ecmech::nsvec; ++jSvec ) { mtanSD(jSvec,2) =           sqr2b3 * C(jSvec,1)  + t3_vec[jSvec] ; }
   if ( l_ddsdde_gamma ) {
      // extra factor of 1/2 for shear deformation transformation
      for ( int jSvec=0; jSvec<ecmech::nsvec; ++jSvec ) { mtanSD(jSvec,3) =  C(jSvec,4) * sqr2i ; }
      for ( int jSvec=0; jSvec<ecmech::nsvec; ++jSvec ) { mtanSD(jSvec,4) =  C(jSvec,3) * sqr2i ; }
      for ( int jSvec=0; jSvec<ecmech::nsvec; ++jSvec ) { mtanSD(jSvec,5) =  C(jSvec,2) * sqr2i ; }
   } else {
      for ( int jSvec=0; jSvec<ecmech::nsvec; ++jSvec ) { mtanSD(jSvec,3) =  C(jSvec,4) * sqr2 ; }
      for ( int jSvec=0; jSvec<ecmech::nsvec; ++jSvec ) { mtanSD(jSvec,4) =  C(jSvec,3) * sqr2 ; }
      for ( int jSvec=0; jSvec<ecmech::nsvec; ++jSvec ) { mtanSD(jSvec,5) =  C(jSvec,2) * sqr2 ; }
   }
   
} // mtan_conv_sd_svec

#ifdef __cuda_host_only__
template< int n >
inline void
printVec(const real8* const y, std::ostream & oss ) {
   oss << std::setprecision(14) ;
   for ( int iX=0; iX<n; ++iX) {
      oss << y[iX] << " " ;
   }
   oss << std::endl ;
}

template< int n >
inline void
printMat(const real8* const A, std::ostream & oss ) {
   oss << std::setprecision(14) ;
   for ( int iX=0; iX<n; ++iX) {
      for ( int jX=0; jX<n; ++jX) {
         oss << A[ECMECH_NN_INDX(iX,jX,n)] << " " ;
      }
      oss << std::endl ;
   } 
}
#endif

} // namespace ecmech

#endif  // ECMECH_UTIL_H
