// -*-c++-*-

#include "ECMech_slipgeom.h"
#include "ECMech_util.h"

namespace ecmech {
   __ecmech_hdev__
   static void
   fillFromMS(double* const P, // ntvec * nslip
              double* const Q, // nwvec * nslip
              const double* const mVecs, // nslip * ndim
              const double* const sVecs, // nslip * ndim
              int nslip)
   {
      for (int iSlip = 0; iSlip<nslip; ++iSlip) {
         // CALL vec_x_vect_mn(crys%vecs(:,is),crys%vecm(:,is),crys%t_ref(:,:,is),DIMS,DIMS)
         double T_ref[ ecmech::ndim*ecmech::ndim ];
         vecsMaTb<ndim>(T_ref, &(sVecs[iSlip * ecmech::ndim]), &(mVecs[iSlip * ecmech::ndim]) );

         double P_vecd[ ecmech::ntvec ];
         double Q_veccp[ ecmech::nwvec ];
         matToPQ(P_vecd, Q_veccp, T_ref);

         for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
            P[ECMECH_NM_INDX(iTvec, iSlip, ecmech::ntvec, nslip)] = P_vecd[iTvec];
         }

         for (int iWvec = 0; iWvec < ecmech::nwvec; ++iWvec) {
            Q[ECMECH_NM_INDX(iWvec, iSlip, ecmech::nwvec, nslip)] = Q_veccp[iWvec];
         }

         //
         // in some approaches, it is useful to form the outer product of P_vecd with itself, for tangent stiffness contributions
      }
   }

   // constructor
   __ecmech_hdev__
   SlipGeomFCC::SlipGeomFCC()
   {
      // m = (/ sqr3i, sqr3i, sqr3i /)
      // s = (/ zero, sqr2i, -sqr2i /)
      //
      // do not yet bother with making slip systems from symmetry group -- just write them out
      const double mVecs[ nslip * ecmech::ndim ] = {
         sqr3i, sqr3i, sqr3i,
         sqr3i, sqr3i, sqr3i,
         sqr3i, sqr3i, sqr3i,
         -sqr3i, sqr3i, sqr3i,
         -sqr3i, sqr3i, sqr3i,
         -sqr3i, sqr3i, sqr3i,
         -sqr3i, -sqr3i, sqr3i,
         -sqr3i, -sqr3i, sqr3i,
         -sqr3i, -sqr3i, sqr3i,
         sqr3i, -sqr3i, sqr3i,
         sqr3i, -sqr3i, sqr3i,
         sqr3i, -sqr3i, sqr3i,
      };
      const double sVecs[ nslip * ecmech::ndim ] = {
         zero, sqr2i, -sqr2i,
         -sqr2i, zero, sqr2i,
         sqr2i, -sqr2i, zero,
         -sqr2i, zero, -sqr2i,
         zero, -sqr2i, sqr2i,
         sqr2i, sqr2i, zero,
         zero, -sqr2i, -sqr2i,
         sqr2i, zero, sqr2i,
         -sqr2i, sqr2i, zero,
         sqr2i, zero, -sqr2i,
         zero, sqr2i, sqr2i,
         -sqr2i, -sqr2i, zero,
      };

      fillFromMS(this->_P_ref_vec, this->_Q_ref_vec,
                 mVecs, sVecs, this->nslip);
   };
} // namespace ecmech
