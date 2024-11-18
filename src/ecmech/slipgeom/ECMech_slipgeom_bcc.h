#pragma once

#include "ECMech_slipgeom_base.h"

namespace ecmech {

   /**
    * BCC with 12, 24, or 48 slip systems
    *
    */
   template<int nSlipTmplt>
   class SlipGeomBCC : public SlipGeom<nSlipTmplt>
   {
      private:
         static constexpr int nslipAddBase = 12;
         static constexpr int nslipAddPGa = 12;
         static constexpr int nslipAddPGb = 24;

      public:
         static const bool dynamic = false;
         static constexpr int nslip = nSlipTmplt;
         static constexpr int nParams = 0;

         static constexpr int nslipBase = nslipAddBase;
         static constexpr int nslipPGa = nslipAddBase + nslipAddPGa;
         static constexpr int nslipPGb = nslipAddBase + nslipAddPGa + nslipAddPGb;

         // constructor and destructor
         __ecmech_hdev__
         SlipGeomBCC() {
            assert(nslip == nslipBase || nslip == nslipPGa || nslip == nslipPGb);
         }
         __ecmech_hdev__
         ~SlipGeomBCC() {}

         __ecmech_hdev__
         SlipGeomBCC(const double* const params) {
            setParams(params);
         }

         __ecmech_host__
         void setParams(const std::vector<double> & params)
         {
            setParams(params.data());
         }

         __ecmech_hdev__
         void setParams(const double* const)
         {
            double mVecs[nslipPGb * ecmech::ndim] = {};
            double sVecs[nslipPGb * ecmech::ndim] = {};

            auto add_vec_data = [=] (const double* const array_src, double* const array_dst, const size_t length) {
               for(size_t ivec = 0; ivec < length; ivec++)
               {
                  array_dst[ivec] = array_src[ivec];
               }
            };

            {
               // m = (/ zero, sqr2i, -sqr2i /)
               // s = (/ sqr3i, sqr3i, sqr3i /)
               const int nslipThese = nslipAddBase;
               //
               // do not yet bother with making slip systems from symmetry group -- just write them out
               const double P3 = sqr3i, M3 = -sqr3i;
               const double P2 = sqr2i, M2 = -sqr2i;
               const double Z = zero;
               //#Slip direction CUB111
               const double sVecsThese[ nslip * ecmech::ndim ] = {
                     P3, P3, P3,
                     P3, P3, P3,
                     P3, P3, P3,
                     P3, P3, M3,
                     P3, P3, M3,
                     P3, P3, M3,
                     P3, M3, P3,
                     P3, M3, P3,
                     P3, M3, P3,
                     P3, M3, M3,
                     P3, M3, M3,
                     P3, M3, M3};
               //#Slip plane normal CUB110
               const double mVecsThese[ nslip * ecmech::ndim ] = {
                     Z,  P2, M2,
                     P2, Z,  M2,
                     P2, M2, Z,
                     Z,  P2, P2,
                     P2, Z,  P2,
                     P2, M2, Z,
                     Z,  P2, P2,
                     P2, Z,  M2,
                     P2, P2, Z,
                     Z,  P2, M2,
                     P2, Z,  P2,
                     P2, P2, Z};

               add_vec_data(mVecsThese, mVecs, nslipThese * ecmech::ndim);
               add_vec_data(sVecsThese, sVecs, nslipThese * ecmech::ndim);
            }

            if (nslip >= nslipPGa) {
               const double twSqr6i = 2.0 * sqr6i;

               // 12 {112}<111> slip systems
               const int nslipThese = nslipAddPGa;

               const double mVecsThese[ nslipThese * ecmech::ndim ] = {
                  -twSqr6i, sqr6i, sqr6i,
                  sqr6i, -twSqr6i, sqr6i,
                  sqr6i, sqr6i, -twSqr6i,
                  -sqr6i, -twSqr6i, sqr6i,
                  twSqr6i, sqr6i, sqr6i,
                  -sqr6i, sqr6i, -twSqr6i,
                  twSqr6i, -sqr6i, sqr6i,
                  -sqr6i, twSqr6i, sqr6i,
                  -sqr6i, -sqr6i, -twSqr6i,
                  sqr6i, twSqr6i, sqr6i,
                  -twSqr6i, -sqr6i, sqr6i,
                  sqr6i, -sqr6i, -twSqr6i,
               };
               const double sVecsThese[ nslipThese * ecmech::ndim ] = {
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
               add_vec_data(mVecsThese, &mVecs[nslipAddBase], nslipThese * ecmech::ndim);
               add_vec_data(sVecsThese, &sVecs[nslipAddBase], nslipThese * ecmech::ndim);
            }

            if (nslip >= nslipPGb) {
               const double mPg2a = 1.0 / sqrt(14.0);
               const double mPg2b = 2.0 / sqrt(14.0);
               const double mPg2c = 3.0 / sqrt(14.0);

               // 24 {123}<111> slip systems
               const int nslipThese = nslipAddPGb;

               const double mVecsThese[ nslipThese * ecmech::ndim ] = {
                  mPg2c, -mPg2a, -mPg2b,
                  -mPg2b, mPg2c, -mPg2a,
                  -mPg2a, -mPg2b, mPg2c,
                  mPg2a, mPg2c, -mPg2b,
                  -mPg2c, -mPg2b, -mPg2a,
                  mPg2b, -mPg2a, mPg2c,
                  -mPg2c, mPg2a, -mPg2b,
                  mPg2b, -mPg2c, -mPg2a,
                  mPg2a, mPg2b, mPg2c,
                  -mPg2a, -mPg2c, -mPg2b,
                  mPg2c, mPg2b, -mPg2a,
                  -mPg2b, mPg2a, mPg2c,
                  -mPg2a, mPg2c, mPg2b,
                  mPg2c, -mPg2b, mPg2a,
                  -mPg2b, -mPg2a, -mPg2c,
                  -mPg2c, -mPg2a, mPg2b,
                  mPg2b, mPg2c, mPg2a,
                  mPg2a, -mPg2b, -mPg2c,
                  mPg2a, -mPg2c, mPg2b,
                  -mPg2c, mPg2b, mPg2a,
                  mPg2b, mPg2a, -mPg2c,
                  mPg2c, mPg2a, mPg2b,
                  -mPg2b, -mPg2c, mPg2a,
                  -mPg2a, mPg2b, -mPg2c,
               };
               const double sVecsThese[ nslipThese * ecmech::ndim ] = {
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
                  sqr3i, sqr3i, -sqr3i,
                  sqr3i, sqr3i, -sqr3i,
                  sqr3i, sqr3i, -sqr3i,
                  -sqr3i, sqr3i, -sqr3i,
                  -sqr3i, sqr3i, -sqr3i,
                  -sqr3i, sqr3i, -sqr3i,
                  -sqr3i, -sqr3i, -sqr3i,
                  -sqr3i, -sqr3i, -sqr3i,
                  -sqr3i, -sqr3i, -sqr3i,
                  sqr3i, -sqr3i, -sqr3i,
                  sqr3i, -sqr3i, -sqr3i,
                  sqr3i, -sqr3i, -sqr3i,
               };
               add_vec_data(mVecsThese, &mVecs[nslipAddPGa], nslipThese * ecmech::ndim);
               add_vec_data(sVecsThese, &sVecs[nslipAddPGa], nslipThese * ecmech::ndim);
            }

            fillFromMS(this->m_P_ref_vec, this->m_Q_ref_vec, mVecs, sVecs, this->nslip);

            for (int i = 0; i < nslip * ecmech::ndim; i++) {
               this->m_s_ref_vec[i] = sVecs[i];
               this->m_m_ref_vec[i] = mVecs[i];
            }
         }

         __ecmech_host__
         void getParams(std::vector<double> & /* params */
                        ) const {
            // do not clear params in case adding to an existing set
         }

   }; // SlipGeomBCC

   class SlipGeomBCCPencil : public SlipGeom<4, 4>
   {
      public:

         static const bool dynamic = true;
         static constexpr int nParams = 0;

         // constructor and destructor
         SlipGeomBCCPencil() = default;
         __ecmech_hdev__
         ~SlipGeomBCCPencil() {}

         __ecmech_hdev__
         SlipGeomBCCPencil(const double* const params) {
            setParams(params);
         }

         __ecmech_host__
         void setParams(const std::vector<double> & params)
         {
            setParams(params.data());
         }

         __ecmech_hdev__
         void setParams(const double* const)
         {
            // s = (/ sqr3i, sqr3i, sqr3i /)
            //
            const double P3 = sqr3i, M3 = -sqr3i;
            const double P2 = sqr2i, M2 = -sqr2i;
            const double Z = zero;
         
            const double sVecs[ nslip * ecmech::ndim ] = {
               P3, P3, P3,
               M3, P3, P3,
               P3, M3, P3,
               P3, P3, M3};
               
            const double mVecs[ nslip * ecmech::ndim ] = {
               Z, M2, P2,
               Z, M2, P2,
               Z, P2, P2,
               Z, P2, P2};

            fillFromMS(this->m_P_ref_vec, this->m_Q_ref_vec,
                       mVecs, sVecs, this->nslip);

            for (int i = 0; i < nslip * ecmech::ndim; i++) {
               m_s_ref_vec[i] = sVecs[i];
               m_m_ref_vec[i] = mVecs[i];
            }
         }

         __ecmech_host__
         void getParams(std::vector<double> & /* params */
                        ) const {
            // do not clear params in case adding to an existing set
         }

         __ecmech_hdev__ inline void getPQ(double* chia, 
                                           double* P_vec, 
                                           double* Q_vec, 
                                           const double* const SvecP) const override final
         {
             double eps = 1e-10;
             double mVecs[nslip * ecmech::ndim];
             
             double S[ecmech::ndim * ecmech::ndim];
             // Svec: 11' 22' 33' 23 31 12 p
             S[ECMECH_NN_INDX(0, 0, 3)] = SvecP[0] + ecmech::onethird * SvecP[6];
             S[ECMECH_NN_INDX(1, 1, 3)] = SvecP[1] + ecmech::onethird * SvecP[6];
             S[ECMECH_NN_INDX(2, 2, 3)] = SvecP[2] + ecmech::onethird * SvecP[6];
             S[ECMECH_NN_INDX(1, 2, 3)] = S[ECMECH_NN_INDX(2, 1, 3)] = SvecP[3];
             S[ECMECH_NN_INDX(2, 0, 3)] = S[ECMECH_NN_INDX(0, 2, 3)] = SvecP[4];
             S[ECMECH_NN_INDX(0, 1, 3)] = S[ECMECH_NN_INDX(1, 0, 3)] = SvecP[5];
             
             for (int iSlip = 0; iSlip < nslip; ++iSlip) {
                 const double* sVec = &m_s_ref_vec[iSlip * ecmech::ndim];
                 
                 // PK force direction
                 double fpk[ecmech::ndim] = {0.0};
                 double Sb[ecmech::ndim];
                 vecsVMa<ecmech::ndim>(Sb, S, sVec);
                 if (vecNorm<ecmech::ndim>(Sb) > eps) {
                     vecCrossProd(fpk, Sb, sVec);
                 }
                 
                 // Normal direction
                 double* mVec = &mVecs[iSlip * ecmech::ndim];
                 if (vecNorm<ecmech::ndim>(fpk) > eps) {
                     vecCrossProd(mVec, sVec, fpk);
                     vecsVNormalize<ecmech::ndim>(mVec);
                 } else {
                     for (int i = 0; i < ecmech::ndim; i++)
                        mVec[i] = m_m_ref_vec[iSlip * ecmech::ndim + i];
                 }
                 
                 // MRSSP angle
                 double n0Vec[ecmech::ndim] = { //n0 = 1/sqrt(2)*(2*b[0],-b[1],-b[2])
                      2.0*sqr2i*sVec[0],
                     -1.0*sqr2i*sVec[1],
                     -1.0*sqr2i*sVec[2]
                 };
                 double t0Vec[ecmech::ndim];
                 vecCrossProd(t0Vec, n0Vec, sVec);
                 double fx = vecsyadotb<ecmech::ndim>(fpk, t0Vec);
                 double fy = vecsyadotb<ecmech::ndim>(fpk, n0Vec);
                 double chi = atan2(fy, fx)-M_PI/6.0;
                 // Fold into T/AT primary region (-30:30)
                 if (chi > 1.0*M_PI/6.0 && chi <= 3.0*M_PI/6.0) {
                     chi = M_PI/3.0-chi;
                 } else if (chi > 3.0*M_PI/6.0 && chi <= 5.0*M_PI/6.0) {
                     chi -= 2.0*M_PI/3.0;
                 } else if (chi >= -7.0*M_PI/6.0 && chi < -5.0*M_PI/6.0) {
                     chi = -M_PI-chi;
                 } else if (chi >= -5.0*M_PI/6.0 && chi < -3.0*M_PI/6.0) {
                     chi += 2.0*M_PI/3.0;
                 } else if (chi >= -3.0*M_PI/6.0 && chi < -1.0*M_PI/6.0) {
                     chi = -M_PI/3.0-chi;
                 }
                 chia[iSlip] = chi;
                 m_chia[iSlip] = chia[iSlip];
             }
             
             fillFromMS(P_vec, Q_vec, mVecs, m_s_ref_vec, nslip);
             fillFromMS(m_P_vec, m_Q_vec, mVecs, m_s_ref_vec, nslip);
         }
         __ecmech_hdev__ inline virtual const double* getP() const override final { return m_P_vec; }
         __ecmech_hdev__ inline virtual const double* getQ() const override final { return m_Q_vec; }
         __ecmech_hdev__
         inline
         void getExtras(double* const chia) const
         {
            for (size_t islip = 0; islip < nslip; islip++) {
               chia[islip] = m_chia[islip];
            }
         }
      private:
         mutable double m_P_vec[ ecmech::ntvec * nslip ];
         mutable double m_Q_vec[ ecmech::nwvec * nslip ];
         mutable double m_chia[ nslip ];

   }; // SlipGeomBCCPencil
   
   // Note the slip systems within this model have slight variations of the slip planes compared
   // to the slip planes in the normal SlipGeomBCC<12> model...
   // Due to this some models could end up with slightly different answers when comparing the two cases.
   class SlipGeomBCCNonSchmid : public SlipGeom<12>
   {
      public:
         static const bool dynamic = true;
         static constexpr int nParams = 3;
         static constexpr size_t nSlipExtra = 0;

         // constructor and destructor
         SlipGeomBCCNonSchmid() = default;
         __ecmech_hdev__
         ~SlipGeomBCCNonSchmid() {}

         __ecmech_hdev__
         SlipGeomBCCNonSchmid(const double* const params) {
            setParams(params);
         }

         __ecmech_host__
         void setParams(const std::vector<double> & params)
         {
            setParams(params.data());
         }

         __ecmech_hdev__
         void setParams(const double* const params)
         {
            const double* parsIt = params;

            m_omegas[0] = *parsIt; ++parsIt;
            m_omegas[1] = *parsIt; ++parsIt;
            m_omegas[2] = *parsIt; ++parsIt;
            const double sum_omegas = fabs(m_omegas[0]) + fabs(m_omegas[1]) + fabs(m_omegas[2]);
            m_isotropic = sum_omegas < ecmech::idp_eps_sqrt;

            const double P3 = sqr3i, M3 = -sqr3i;
            const double P2 = sqr2i, M2 = -sqr2i;
            const double Z = zero;
            
            const double sVecs[ nslip * ecmech::ndim ] = {
                 P3, P3, P3,
                 P3, P3, P3,
                 P3, P3, P3,
                 P3, P3, M3,
                 P3, P3, M3,
                 P3, P3, M3,
                 P3, M3, P3,
                 P3, M3, P3,
                 P3, M3, P3,
                 P3, M3, M3,
                 P3, M3, M3,
                 P3, M3, M3};
           
            // This list of planes has been generated to be 
            // consistent with T/AT directions
            const double mVecs[ nslip * ecmech::ndim ] = {
                 P2, M2,  Z,
                 Z,  P2, M2,
                 M2,  Z, P2,
                 P2,  Z, P2,
                 M2, P2,  Z,
                 Z,  M2, M2,
                 P2,  Z, M2,
                 M2, M2,  Z,
                 Z,  P2, P2,
                 P2, P2,  Z,
                 Z,  M2, P2,
                 M2,  Z, M2};
            
            fillFromMS(this->m_P_ref_vec, this->m_Q_ref_vec,
                       mVecs, sVecs, this->nslip);

            for (int i = 0; i < nslip * ecmech::ndim; i++) {
               m_s_ref_vec[i] = sVecs[i];
               m_m_ref_vec[i] = mVecs[i];
            }

#if defined(ECMECH_DEBUG)
            int iParam = parsIt - params;
            if (iParam != nParams) {
               ECMECH_FAIL(__func__, "iParam != nParams");
            }
#endif
         }

         __ecmech_host__
         void getParams(std::vector<double> & params
                        ) const {
#ifdef ECMECH_DEBUG
            // do not clear params in case adding to an existing set
            int paramsStart = params.size();
#endif
            params.push_back(m_omegas[0]);
            params.push_back(m_omegas[1]);
            params.push_back(m_omegas[2]);

#ifdef ECMECH_DEBUG
            assert((params.size() - paramsStart) == nParams);
#endif
         }
         
         __ecmech_hdev__ inline void NSprojection(double* taua,
                                                  double* P_vec, 
                                                  double* Q_vec, 
                                                  const double* const kirchoff,
                                                  bool fill_PQ) const
         {
             // Resolve stress onto slip systems
             // Compute RSS considering both senses of the slip direction
             // and select the most favorable direction
             for (int iSlip = 0; iSlip < nslip; ++iSlip) {
                 double P_s[2 * ecmech::ntvec];
                 double Q_s[2 * ecmech::nwvec];
                 double tau_s[2] = { 0.0 };
                 
                 for (int iS = 0; iS < 2; ++iS) {
                     
                     double P_tmp[ecmech::ntvec];
                     double Q_tmp[ecmech::nwvec];
                     
                     const double* mVec = &m_m_ref_vec[iSlip * ecmech::ndim];
                     double sVec[ecmech::ndim];
                     for (int i = 0; i < ecmech::ndim; i++)
                         sVec[i] = (1 - 2*iS) * m_s_ref_vec[iSlip * ecmech::ndim + i];
                         
                     // Schmid
                     fillFromMS(P_tmp, Q_tmp, mVec, sVec, 1);
                     for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                         P_s[iS * ecmech::ntvec + iTvec] = P_tmp[iTvec];
                         tau_s[iS] += kirchoff[iTvec] * P_tmp[iTvec];
                     }
                     for (int iWvec = 0; iWvec < ecmech::nwvec; ++iWvec) {
                         Q_s[iS * ecmech::nwvec + iWvec] = Q_tmp[iWvec];
                     }
                     
                     // Non-Schmid
                     double smVec[ecmech::ndim];
                     double mpVec[ecmech::ndim];
                     vecCrossProd(smVec, sVec, mVec);
                     for (int i = 0; i < ecmech::ndim; i++)
                         mpVec[i] = 0.5*mVec[i] + 0.8660254037844386*smVec[i];
                         
                     // Omega 1
                     fillFromMS(P_tmp, Q_tmp, mpVec, sVec, 1);
                     for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                         tau_s[iS] += m_omegas[0] * kirchoff[iTvec] * P_tmp[iTvec];
                     }
                     
                     // Omega 2
                     fillFromMS(P_tmp, Q_tmp, mVec, smVec, 1);
                     for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                         tau_s[iS] -= m_omegas[1] * kirchoff[iTvec] * P_tmp[iTvec];
                     }
                     
                     // Omega 3
                     double mpsVec[ecmech::ndim];
                     vecCrossProd(mpsVec, mpVec, sVec);
                     fillFromMS(P_tmp, Q_tmp, mpVec, mpsVec, 1);
                     for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                         tau_s[iS] += m_omegas[2] * kirchoff[iTvec] * P_tmp[iTvec];
                     }
                 }
                 
                 // Keep highest value
                 int iS = (int)(tau_s[1] > tau_s[0]);
                 
                 taua[iSlip] = tau_s[iS];
                 if (taua[iSlip] < 0.0) taua[iSlip] = 0.0;
                 
                 if (fill_PQ) {
                     for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                         P_vec[ECMECH_NM_INDX(iTvec, iSlip, ecmech::ntvec, nslip)] = P_s[iS * ecmech::ntvec + iTvec];
                         m_P_vec[ECMECH_NM_INDX(iTvec, iSlip, ecmech::ntvec, nslip)] = P_s[iS * ecmech::ntvec + iTvec];
                     }
                     for (int iWvec = 0; iWvec < ecmech::nwvec; ++iWvec) {
                         Q_vec[ECMECH_NM_INDX(iWvec, iSlip, ecmech::nwvec, nslip)] = Q_s[iS * ecmech::nwvec + iWvec];
                         m_Q_vec[ECMECH_NM_INDX(iWvec, iSlip, ecmech::nwvec, nslip)] = Q_s[iS * ecmech::nwvec + iWvec];
                     }
                 }
             }
         }
         
         __ecmech_hdev__ inline void getPQ(double* /*chia*/,
                                           double* P_vec, 
                                           double* Q_vec, 
                                           const double* const SvecP) const override final
         {
            if (m_isotropic) {
                  for (int iTvec = 0; iTvec < ecmech::ntvec * nslip; ++iTvec) {
                        P_vec[iTvec] = m_P_ref_vec[iTvec];
                        m_P_vec[iTvec] = m_P_ref_vec[iTvec];
                  }
                  for (int iWvec = 0; iWvec < ecmech::nwvec * nslip; ++iWvec) {
                        Q_vec[iWvec] = m_Q_ref_vec[iWvec];
                        m_Q_vec[iWvec] = m_Q_ref_vec[iWvec];
                  }
            } else {
               // we need to reverse the stress first...
               double kirchoff[ecmech::nsvec];
               kirchoff[iSvecS] = -sqr3 * SvecP[iSvecP];
               kirchoff[0] = sqr2i * SvecP[0] - sqr2i * SvecP[1];
               kirchoff[1] = - sqr3b2 * SvecP[0] - sqr3b2 * SvecP[1];
               kirchoff[4] = sqr2 * SvecP[3]; // 23
               kirchoff[3] = sqr2 * SvecP[4]; // 31
               kirchoff[2] = sqr2 * SvecP[5]; // 12
               
               double taua[nslip];
               NSprojection(taua, P_vec, Q_vec, kirchoff, true);
            }
         }
         
         __ecmech_hdev__ inline void evalRSS(double* taua, 
                                             const double* const kirchoff, 
                                             const double* /*P_vec*/) const override final
         {
            if (m_isotropic) {
               SlipGeom::evalRSS(taua, kirchoff, m_P_ref_vec);
            } else {
               NSprojection(taua, NULL, NULL, kirchoff, false);
            }
         }

         __ecmech_hdev__ inline virtual const double* getP() const override final { return m_P_vec; }
         __ecmech_hdev__ inline virtual const double* getQ() const override final { return m_Q_vec; }
         __ecmech_hdev__
         inline
         void getExtras(double* const) const {}
      private:
         mutable double m_P_vec[ ecmech::ntvec * nslip ];
         mutable double m_Q_vec[ ecmech::nwvec * nslip ];
         double m_omegas[3];
         bool m_isotropic = false;

   }; // SlipGeomBCCNonSchmid

}