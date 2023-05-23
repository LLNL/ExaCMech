// -*-c++-*-

#ifndef ECMECH_SLIPGEOM_H
#define ECMECH_SLIPGEOM_H

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
   
   template<int Nslip>
   class SlipGeom {
      public:
         static const int nslip = Nslip;
          
         __ecmech_hdev__ inline const double* getP() const { return _P_ref_vec; };
         __ecmech_hdev__ inline const double* getQ() const { return _Q_ref_vec; };
         __ecmech_hdev__ inline const double* getM() const { return _m_ref_vec; };
         __ecmech_hdev__ inline const double* getS() const { return _s_ref_vec; };
         
         __ecmech_hdev__ inline virtual void getPQ(double* chia, 
                                                   double* _P_vec, 
                                                   double* _Q_vec, 
                                                   const double* const /* SvecP = nullptr */) const 
         {
             for (int iSlip = 0; iSlip < nslip; iSlip++) {
                 chia[iSlip] = ecmech::zero;
             }
             for (int iTvec = 0; iTvec < ecmech::ntvec * nslip; ++iTvec) {
                 _P_vec[iTvec] = _P_ref_vec[iTvec];
             }
             for (int iWvec = 0; iWvec < ecmech::nwvec * nslip; ++iWvec) {
                 _Q_vec[iWvec] = _Q_ref_vec[iWvec];
             }
         };
         
         __ecmech_hdev__ inline virtual void evalRSS(double* taua, 
                                                     const double* const T_vecds, 
                                                     const double* P_vec) const
         {
             // resolve stress onto slip systems
             vecsVaTM<ecmech::ntvec, nslip>(taua, T_vecds, P_vec);
         }
       
      protected:
         double _m_ref_vec[ ecmech::ndim * nslip ];
         double _s_ref_vec[ ecmech::ndim * nslip ];
         double _P_ref_vec[ ecmech::ntvec * nslip ];
         double _Q_ref_vec[ ecmech::nwvec * nslip ];
   };
   
   
   class SlipGeomFCC : public SlipGeom<12>
   {
      public:
         static const bool dynamic = false;
         static const int nParams = 0;

         // constructor and destructor
         __ecmech_hdev__  SlipGeomFCC() {};
         __ecmech_hdev__ ~SlipGeomFCC() {};

         __ecmech_host__
         void setParams(const std::vector<double> & /* params */
                        )
         {

            // m = (/ sqr3i, sqr3i, sqr3i /)
            // s = (/ zero, sqr2i, -sqr2i /)
            //
            // do not yet bother with making slip systems from symmetry group -- just write them out
         const double P3 = sqr3i, M3 = -sqr3i;
         const double P2 = sqr2i, M2 = -sqr2i;
         const double Z = zero;
         //#Slip plane normal CUB111
         const double mVecs[ nslip * ecmech::ndim ] = {
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
            //#Slip direction CUB110
         const double sVecs[ nslip * ecmech::ndim ] = {
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

            fillFromMS(this->_P_ref_vec, this->_Q_ref_vec,
                       mVecs, sVecs, this->nslip);

            for (int i = 0; i < nslip * ecmech::ndim; i++) {
               _s_ref_vec[i] = sVecs[i];
               _m_ref_vec[i] = mVecs[i];
            }

            // assert((parsIt - params.begin()) == nParams);
         };

         __ecmech_host__
         void getParams(std::vector<double> & /* params */
                        ) const {
            // do not clear params in case adding to an existing set
         }
         
   }; // SlipGeomFCC

   /**
    * BCC with 12, 24, or 48 slip systems
    *
    */
   template<int nSlipTmplt>
   class SlipGeomBCC : public SlipGeom<nSlipTmplt>
   {
      private:
         static const int _nslipAddBase = 12;
         static const int _nslipAddPGa = 12;
         static const int _nslipAddPGb = 24;

      public:
         static const bool dynamic = false;
         static const int nslip = nSlipTmplt;
         static const int nParams = 0;

         static const int nslipBase = _nslipAddBase;
         static const int nslipPGa = _nslipAddBase + _nslipAddPGa;
         static const int nslipPGb = _nslipAddBase + _nslipAddPGa + _nslipAddPGb;

         // constructor and destructor
         __ecmech_hdev__  SlipGeomBCC() {
            assert(nslip == nslipBase || nslip == nslipPGa || nslip == nslipPGb);
         };
         __ecmech_hdev__ ~SlipGeomBCC() {};

         __ecmech_host__
         void setParams(const std::vector<double> & /* params */
                        )
         {

            std::vector<double> mVecs;
            std::vector<double> sVecs;

            {
               // m = (/ zero, sqr2i, -sqr2i /)
               // s = (/ sqr3i, sqr3i, sqr3i /)
               const int nslipThese = _nslipAddBase;
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

               mVecs.insert(mVecs.end(), &(mVecsThese[0]), &(mVecsThese[nslipThese * ecmech::ndim]));
               sVecs.insert(sVecs.end(), &(sVecsThese[0]), &(sVecsThese[nslipThese * ecmech::ndim]));
            }

            if (nslip >= nslipPGa) {
               const double twSqr6i = 2.0 * sqr6i;

               // 12 {112}<111> slip systems
               const int nslipThese = _nslipAddPGa;

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
               mVecs.insert(mVecs.end(), &(mVecsThese[0]), &(mVecsThese[nslipThese * ecmech::ndim]));
               sVecs.insert(sVecs.end(), &(sVecsThese[0]), &(sVecsThese[nslipThese * ecmech::ndim]));
            }

            if (nslip >= nslipPGb) {
               const double mPg2a = 1.0 / sqrt(14.0);
               const double mPg2b = 2.0 / sqrt(14.0);
               const double mPg2c = 3.0 / sqrt(14.0);

               // 24 {123}<111> slip systems
               const int nslipThese = _nslipAddPGb;

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
               mVecs.insert(mVecs.end(), &(mVecsThese[0]), &(mVecsThese[nslipThese * ecmech::ndim]));
               sVecs.insert(sVecs.end(), &(sVecsThese[0]), &(sVecsThese[nslipThese * ecmech::ndim]));
            }

            fillFromMS(this->_P_ref_vec, this->_Q_ref_vec,
                       &(mVecs[0]), &(sVecs[0]), this->nslip);

            for (int i = 0; i < nslip * ecmech::ndim; i++) {
               this->_s_ref_vec[i] = sVecs.at(i);
               this->_m_ref_vec[i] = mVecs.at(i);
            }

            // assert((parsIt - params.begin()) == nParams);
         };

         __ecmech_host__
         void getParams(std::vector<double> & /* params */
                        ) const {
            // do not clear params in case adding to an existing set
         }

   }; // SlipGeomBCC

   /**
    * HCP with <a> slip on basal, prisamtic, and pyramidal families and type-1 <c+a> pyramidal slip
    *
    * the name aBRYcaY1 traces back to EVP_HCP_a_BRY_ca_Y1 (integer code 32) in the old Fortran coding
    *
    * fix me : the coding below is a hack just to get things going ;
    * it is not the best way of doing things, and modifications should be made with great care
    */
   class SlipGeomHCPaBRYcaY1 : public SlipGeom<3 + 3 + 6 + 12>
   {
      public:
         static const bool dynamic = false;
         // 3  slip systems in basal <a> family
         // 3  slip systems in prismatic <a> family
         // 6  slip systems in pyramidal <a> family
         // 12  slip systems in pyramidal 1 <c+a> family
         //static const int nslip = 3 + 3 + 6 + 12;
         static const int nParams = 1;

         // constructor and destructor
         __ecmech_hdev__  SlipGeomHCPaBRYcaY1() {};
         __ecmech_hdev__ ~SlipGeomHCPaBRYcaY1() {};

         __ecmech_host__
         void setParams(const std::vector<double> & params
                        )
         {
            std::vector<double>::const_iterator parsIt = params.begin();

            _cOverA = *parsIt; ++parsIt;

            // pyramidal 10-11 1-210 depends on c/a
            //
            double m_ya[ecmech::ndim], s_ya[ecmech::ndim];
            {
               double an[ecmech::nMiller] = { one, zero, -one, one }; // plane
               double ab[ecmech::nMiller] = { one, -two, one, zero }; // direction
               //
               miller_to_orthog_sngl(an, ab,
                                     m_ya, s_ya,
                                     _cOverA);
            }
            double m_ya_pp = sqrt(1.0 - m_ya[2] * m_ya[2]);

            // pyramidal 10-11 -1-123 depends on c/a
            //
            double m_y1ca[ecmech::ndim], s_y1ca[ecmech::ndim];
            {
               double an[ecmech::nMiller] = { one, zero, -one, one }; // plane
               double ab[ecmech::nMiller] = { -one, -one, two, three }; // direction
               //
               miller_to_orthog_sngl(an, ab,
                                     m_y1ca, s_y1ca,
                                     _cOverA);
            }
            double m_y1ca_pp = sqrt(1.0 - m_y1ca[2] * m_y1ca[2]);
            double s_y1ca_pp = sqrt(1.0 - s_y1ca[2] * s_y1ca[2]);

            const double mVecs[ nslip * ecmech::ndim ] = {
               zero, zero, one,
               zero, zero, one,
               zero, zero, one,

               -halfsqr3, onehalf, zero,
               -halfsqr3, -onehalf, zero,
               zero, -one, zero,

               m_ya[0], m_ya[1], m_ya[2],
               m_ya[0], -m_ya[1], -m_ya[2],
               m_ya[0], m_ya[1], -m_ya[2],
               -m_ya[0], m_ya[1], -m_ya[2],
               zero, m_ya_pp, -m_ya[2],
               zero, -m_ya_pp, -m_ya[2],

               m_y1ca[0], m_y1ca[1], m_y1ca[2],
               m_y1ca[0], -m_y1ca[1], -m_y1ca[2],
               m_y1ca[0], m_y1ca[1], -m_y1ca[2],
               zero, m_y1ca_pp, -m_y1ca[2],
               -m_y1ca[0], m_y1ca[1], -m_y1ca[2],
               -m_y1ca[0], -m_y1ca[1], -m_y1ca[2],
               zero, -m_y1ca_pp, -m_y1ca[2],
               zero, m_y1ca_pp, m_y1ca[2],
               -m_y1ca[0], m_y1ca[1], m_y1ca[2],
               -m_y1ca[0], -m_y1ca[1], m_y1ca[2],
               m_y1ca[0], -m_y1ca[1], m_y1ca[2],
               zero, -m_y1ca_pp, m_y1ca[2]
            };
            const double sVecs[ nslip * ecmech::ndim ] = {
               onehalf, halfsqr3, zero,
               onehalf, -halfsqr3, zero,
               one, zero, zero,

               onehalf, halfsqr3, zero,
               onehalf, -halfsqr3, zero,
               one, zero, zero,

               s_ya[0], s_ya[1], zero,
               s_ya[0], -s_ya[1], zero,
               -s_ya[0], -s_ya[1], zero,
               -s_ya[0], s_ya[1], zero,
               -one, zero, zero,
               one, zero, zero,

               s_y1ca[0], s_y1ca[1], s_y1ca[2],
               s_y1ca[0], -s_y1ca[1], -s_y1ca[2],
               -s_y1ca_pp, zero, -s_y1ca[2],
               s_y1ca[0], s_y1ca[1], -s_y1ca[2],
               -s_y1ca[0], s_y1ca[1], -s_y1ca[2],
               s_y1ca_pp, zero, -s_y1ca[2],
               -s_y1ca[0], -s_y1ca[1], -s_y1ca[2],
               -s_y1ca[0], s_y1ca[1], s_y1ca[2],
               s_y1ca_pp, zero, s_y1ca[2],
               -s_y1ca[0], -s_y1ca[1], s_y1ca[2],
               -s_y1ca_pp, zero, s_y1ca[2],
               s_y1ca[0], -s_y1ca[1], s_y1ca[2]
            };

            fillFromMS(this->_P_ref_vec, this->_Q_ref_vec,
                       mVecs, sVecs, this->nslip);

            for (int i = 0; i < nslip * ecmech::ndim; i++) {
               _s_ref_vec[i] = sVecs[i];
               _m_ref_vec[i] = mVecs[i];
            }

            assert((parsIt - params.begin()) == nParams);
         };

         __ecmech_host__
         void getParams(std::vector<double> & params
                        ) const {
#ifdef ECMECH_DEBUG
            // do not clear params in case adding to an existing set
            int paramsStart = params.size();
#endif

            params.push_back(_cOverA);

#ifdef ECMECH_DEBUG
            assert((params.size() - paramsStart) == nParams);
#endif
         }

      private:
         double _cOverA;
    
   }; // SlipGeomHCPaBRYcaY1
   
   
   class SlipGeomBCCPencil : public SlipGeom<4>
   {
      public:

         static const bool dynamic = true;
         static const int nParams = 0;

         // constructor and destructor
         __ecmech_hdev__  SlipGeomBCCPencil() {};
         __ecmech_hdev__ ~SlipGeomBCCPencil() {};

         __ecmech_host__
         void setParams(const std::vector<double> & /* params */
                        )
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

            fillFromMS(this->_P_ref_vec, this->_Q_ref_vec,
                       mVecs, sVecs, this->nslip);

            for (int i = 0; i < nslip * ecmech::ndim; i++) {
               _s_ref_vec[i] = sVecs[i];
               _m_ref_vec[i] = mVecs[i];
            }

            // assert((parsIt - params.begin()) == nParams);
         };

         __ecmech_host__
         void getParams(std::vector<double> & /* params */
                        ) const {
            // do not clear params in case adding to an existing set
         }

         __ecmech_hdev__ inline void getPQ(double* chia, 
                                           double* _P_vec, 
                                           double* _Q_vec, 
                                           const double* const SvecP) const
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
                 const double* sVec = &_s_ref_vec[iSlip * ecmech::ndim];
                 
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
                        mVec[i] = _m_ref_vec[iSlip * ecmech::ndim + i];
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
             }
             
             fillFromMS(_P_vec, _Q_vec, mVecs, _s_ref_vec, nslip);
         };

   }; // SlipGeomBCCPencil
   
   
   class SlipGeomBCCNonSchmid : public SlipGeom<12>
   {
      public:
         static const bool dynamic = true;
         static const int nParams = 3;

         // constructor and destructor
         __ecmech_hdev__  SlipGeomBCCNonSchmid() {};
         __ecmech_hdev__ ~SlipGeomBCCNonSchmid() {};

         __ecmech_host__
         void setParams(const std::vector<double> & params
                        )
         {
            std::vector<double>::const_iterator parsIt = params.begin();
            
            _omegas[0] = *parsIt; ++parsIt;
            _omegas[1] = *parsIt; ++parsIt;
            _omegas[2] = *parsIt; ++parsIt;
            
             
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
            
            fillFromMS(this->_P_ref_vec, this->_Q_ref_vec,
                       mVecs, sVecs, this->nslip);

            for (int i = 0; i < nslip * ecmech::ndim; i++) {
               _s_ref_vec[i] = sVecs[i];
               _m_ref_vec[i] = mVecs[i];
            }

            assert((parsIt - params.begin()) == nParams);
         };

         __ecmech_host__
         void getParams(std::vector<double> & params
                        ) const {
#ifdef ECMECH_DEBUG
            // do not clear params in case adding to an existing set
            int paramsStart = params.size();
#endif
            params.push_back(_omegas[0]);
            params.push_back(_omegas[1]);
            params.push_back(_omegas[2]);

#ifdef ECMECH_DEBUG
            assert((params.size() - paramsStart) == nParams);
#endif
         }
         
         __ecmech_hdev__ inline void NSprojection(double* taua,
                                                  double* _P_vec, 
                                                  double* _Q_vec, 
                                                  const double* const T_vecds,
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
                     
                     const double* mVec = &_m_ref_vec[iSlip * ecmech::ndim];
                     double sVec[ecmech::ndim];
                     for (int i = 0; i < ecmech::ndim; i++)
                         sVec[i] = (1 - 2*iS) * _s_ref_vec[iSlip * ecmech::ndim + i];
                         
                     // Schmid
                     fillFromMS(P_tmp, Q_tmp, mVec, sVec, 1);
                     for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                         P_s[iS * ecmech::ntvec + iTvec] = P_tmp[iTvec];
                         tau_s[iS] += T_vecds[iTvec] * P_tmp[iTvec];
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
                         tau_s[iS] += _omegas[0] * T_vecds[iTvec] * P_tmp[iTvec];
                     }
                     
                     // Omega 2
                     fillFromMS(P_tmp, Q_tmp, mVec, smVec, 1);
                     for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                         tau_s[iS] -= _omegas[1] * T_vecds[iTvec] * P_tmp[iTvec];
                     }
                     
                     // Omega 3
                     double mpsVec[ecmech::ndim];
                     vecCrossProd(mpsVec, mpVec, sVec);
                     fillFromMS(P_tmp, Q_tmp, mpVec, mpsVec, 1);
                     for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                         tau_s[iS] += _omegas[2] * T_vecds[iTvec] * P_tmp[iTvec];
                     }
                 }
                 
                 // Keep highest value
                 int iS = (int)(tau_s[1] > tau_s[0]);
                 
                 taua[iSlip] = tau_s[iS];
                 if (taua[iSlip] < 0.0) taua[iSlip] = 0.0;
                 
                 if (fill_PQ) {
                     for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                         _P_vec[ECMECH_NM_INDX(iTvec, iSlip, ecmech::ntvec, nslip)] = P_s[iS * ecmech::ntvec + iTvec];
                     }
                     for (int iWvec = 0; iWvec < ecmech::nwvec; ++iWvec) {
                         _Q_vec[ECMECH_NM_INDX(iWvec, iSlip, ecmech::nwvec, nslip)] = Q_s[iS * ecmech::nwvec + iWvec];
                     }
                 }
             }
         }
         
         __ecmech_hdev__ inline void getPQ(double* /*chia*/,
                                           double* _P_vec, 
                                           double* _Q_vec, 
                                           const double* const SvecP) const
         {
             // we need to reverse the stress first...
             double T_vecds[ecmech::nsvec];
             T_vecds[iSvecS] = -sqr3 * SvecP[iSvecP];
             T_vecds[0] = sqr2i * SvecP[0] - sqr2i * SvecP[1];
             T_vecds[1] = - sqr3b2 * SvecP[0] - sqr3b2 * SvecP[1];
             T_vecds[4] = sqr2 * SvecP[3]; // 23
             T_vecds[3] = sqr2 * SvecP[4]; // 31
             T_vecds[2] = sqr2 * SvecP[5]; // 12
             
             double taua[nslip];
             NSprojection(taua, _P_vec, _Q_vec, T_vecds, true);
         }
         
         __ecmech_hdev__ inline void evalRSS(double* taua, 
                                             const double* const T_vecds, 
                                             const double* /*P_vec*/) const
         {
             NSprojection(taua, NULL, NULL, T_vecds, false);
         }
         
     private:
         double _omegas[3];

   }; // SlipGeomBCCNonSchmid
   
   
} // namespace ecmech

#endif // ECMECH_SLIPGEOM_H
