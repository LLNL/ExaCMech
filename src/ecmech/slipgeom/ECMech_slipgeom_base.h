#pragma once

#include "ECMech_core.h"
#include "ECMech_util.h"

namespace ecmech {
   // A version using constexpr could probably be created for the below
   __ecmech_hdev__
   static void
   fillFromMS(double* const P, // ntvec * nslip
              double* const Q, // nwvec * nslip
              const double* const mVecs, // nslip * ndim
              const double* const sVecs, // nslip * ndim
              const int nslip)
   {
      for (int iSlip = 0; iSlip<nslip; ++iSlip) {
         const double* mVec = &(mVecs[iSlip * ecmech::ndim]);
         const double* sVec = &(sVecs[iSlip * ecmech::ndim]);
#ifndef NO_CHECKS
         if (fabs(vecsyadotb<ecmech::ndim>(mVec, sVec)) > idp_eps_sqrt) {
            ECMECH_FAIL(__func__, "internal error");
         }
#endif

         // CALL vec_x_vect_mn(crys%vecs(:,is),crys%vecm(:,is),crys%t_ref(:,:,is),DIMS,DIMS)
         double T_ref[ ecmech::ndim * ecmech::ndim ];
         vecsMaTb<ndim>(T_ref, sVec, mVec);

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
   
   template<size_t num_slip, size_t num_extra = 0>
   class SlipGeom {
      public:
         static constexpr int nslip = num_slip;
         static constexpr int nSlipExtra = num_extra;

         __ecmech_hdev__
         virtual ~SlipGeom(){}
          
         __ecmech_hdev__ inline virtual const double* getP() const { return m_P_ref_vec; }
         __ecmech_hdev__ inline virtual const double* getQ() const { return m_Q_ref_vec; }
         __ecmech_hdev__ inline const double* getM() const { return m_m_ref_vec; }
         __ecmech_hdev__ inline const double* getS() const { return m_s_ref_vec; }

         __ecmech_hdev__ inline virtual void getPQ(double* /* chia */, 
                                                   double* P_vec, 
                                                   double* Q_vec, 
                                                   const double* const /* SvecP = nullptr */) const 
         {
             for (int iTvec = 0; iTvec < ecmech::ntvec * nslip; ++iTvec) {
                 P_vec[iTvec] = m_P_ref_vec[iTvec];
             }
             for (int iWvec = 0; iWvec < ecmech::nwvec * nslip; ++iWvec) {
                 Q_vec[iWvec] = m_Q_ref_vec[iWvec];
             }
         }
         
         __ecmech_hdev__ inline virtual void evalRSS(double* taua, 
                                                     const double* const kirchoff, 
                                                     const double* P_vec) const
         {
             // resolve stress onto slip systems
             vecsVaTM<ecmech::ntvec, nslip>(taua, kirchoff, P_vec);
         }
       
      protected:
         double m_m_ref_vec[ ecmech::ndim * nslip ];
         double m_s_ref_vec[ ecmech::ndim * nslip ];
         double m_P_ref_vec[ ecmech::ntvec * nslip ];
         double m_Q_ref_vec[ ecmech::nwvec * nslip ];
   };
}