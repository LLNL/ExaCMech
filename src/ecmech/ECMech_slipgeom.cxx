// -*-c++-*-

#include "ECMech_slipgeom.h"
#include "ECMech_port.h"

namespace ecmech {

__ecmech_hdev__
static void
fillFromMS( real8* const P, // ntvec * nslip
            real8* const Q, // nwvec * nslip
            const real8* const mVecs, // nslip * ndim
            const real8* const sVecs, // nslip * ndim
            int nslip )
{
   for ( int iSlip=0; iSlip<nslip; ++iSlip) {

      // CALL vec_x_vect_mn(crys%vecs(:,is),crys%vecm(:,is),crys%t_ref(:,:,is),DIMS,DIMS)
      real8 T_ref[ ecmech::ndim*ecmech::ndim ] ;
      vecsMaTb< ndim >( T_ref, &(sVecs[iSlip*ecmech::ndim]), &(mVecs[iSlip*ecmech::ndim]) ) ;

      real8 P_vecd[ ecmech::ntvec ] ;
      real8 Q_veccp[ ecmech::nwvec ] ;
      matToPQ( P_vecd, Q_veccp, T_ref ) ;

      for ( int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec ) {
         P[ECMECH_NM_INDX(iTvec,iSlip,ecmech::ntvec,nslip)] = P_vecd[iTvec] ;
      }
      for ( int iWvec = 0; iWvec < ecmech::nwvec; ++iWvec ) {
         Q[ECMECH_NM_INDX(iWvec,iSlip,ecmech::nwvec,nslip)] = Q_vecd[iWvec] ;
      }

   }
}   
   
// constructor
__ecmech_hdev__
SlipGeomFCCA::SlipGeomFCCA()
{

   //   m = (/ sqr3i, sqr3i, sqr3i /)
   //   s = (/ zero, sqr2i, -sqr2i /)
   //
   // do not yet bother with making slip systems from symmetry group -- just write them out
   const real8 mVecs[ nslip * ecmech::ndim ] = {
      sqr3i ,  sqr3i ,  sqr3i ,
      sqr3i ,  sqr3i ,  sqr3i ,
      sqr3i ,  sqr3i ,  sqr3i ,
     -sqr3i ,  sqr3i ,  sqr3i ,
     -sqr3i ,  sqr3i ,  sqr3i ,
     -sqr3i ,  sqr3i ,  sqr3i ,
     -sqr3i , -sqr3i ,  sqr3i ,
     -sqr3i , -sqr3i ,  sqr3i ,
     -sqr3i , -sqr3i ,  sqr3i ,
      sqr3i , -sqr3i ,  sqr3i ,
      sqr3i , -sqr3i ,  sqr3i ,
      sqr3i , -sqr3i ,  sqr3i ,
   } ;
   const real8 sVecs[ nslip * ecmech::ndim ] = {
      zero  ,  sqr2i ,  -sqr2i ,
     -sqr2i ,  zero  ,   sqr2i ,
      sqr2i , -sqr2i ,   zero  ,
     -sqr2i ,  zero  ,  -sqr2i ,
      zero  , -sqr2i ,   sqr2i ,
      sqr2i ,  sqr2i ,   zero  ,
      zero  , -sqr2i ,  -sqr2i ,
      sqr2i ,  zero  ,   sqr2i ,
     -sqr2i ,  sqr2i ,   zero  ,
      sqr2i ,  zero  ,  -sqr2i ,
      zero  ,  sqr2i ,   sqr2i ,
     -sqr2i , -sqr2i ,   zero  ,
   } ;

   fillFromMS( this->_P_ref_vec, this->_Q_ref_vec,
               mVecs, sVecs, this->nslip ) ;
      
}

} // namespace ecmech

#endif  // ECMECH_SLIPGEOM_H
