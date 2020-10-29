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
              int nslip)
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

   class SlipGeomFCC
   {
      public:

         static const int nslip = 12;
         static const int nParams = 0;

         // constructor and destructor
         __ecmech_hdev__  SlipGeomFCC() {};
         __ecmech_hdev__ ~SlipGeomFCC() {};

         __ecmech_host__
         void setParams(const std::vector<double> & params
                        )
         {
            std::vector<double>::const_iterator parsIt = params.begin();

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

            int iParam = parsIt - params.begin();
            assert(iParam == nParams);
         };

         __ecmech_host__
         void getParams(std::vector<double> & params
                        ) const {
            // do not clear params in case adding to an existing set
            int paramsStart = params.size();

            // params.push_back(); // no parameters

            int iParam = params.size() - paramsStart;
            assert(iParam == nParams);
         }

         __ecmech_hdev__ inline const double* getP() const { return _P_ref_vec; };
         __ecmech_hdev__ inline const double* getQ() const { return _Q_ref_vec; };

      private:
         double _P_ref_vec[ ecmech::ntvec * nslip ];
         double _Q_ref_vec[ ecmech::nwvec * nslip ];
   }; // SlipGeomFCC

   /**
    * BCC with 12, 24, or 48 slip systems
    *
    */
   template<int nSlipTmplt>
   class SlipGeomBCC
   {
      private:
         static const int _nslipAddBase = 12;
         static const int _nslipAddPGa  = 12;
         static const int _nslipAddPGb  = 24;

      public:

         static const int nslip = nSlipTmplt;
         static const int nParams = 0;

         static const int nslipBase = _nslipAddBase;
         static const int nslipPGa  = _nslipAddBase+_nslipAddPGa;
         static const int nslipPGb  = _nslipAddBase+_nslipAddPGa+_nslipAddPGb;

         // constructor and destructor
         __ecmech_hdev__  SlipGeomBCC() {
            assert(nslip == nslipBase || nslip == nslipPGa || nslip == nslipPGb);
         };
         __ecmech_hdev__ ~SlipGeomBCC() {};

         __ecmech_host__
         void setParams(const std::vector<double> & params
                        )
         {
            std::vector<double>::const_iterator parsIt = params.begin();

            std::vector<double> mVecs;
            std::vector<double> sVecs;

            {
               // m = (/ zero, sqr2i, -sqr2i /)
               // s = (/ sqr3i, sqr3i, sqr3i /)
               const int nslipThese = _nslipAddBase;
               //
               // do not yet bother with making slip systems from symmetry group -- just write them out
               const double mVecsThese[ nslipThese * ecmech::ndim ] = {
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

            if ( nslip >= nslipPGa ) {
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

            if ( nslip >= nslipPGb ) {
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

            int iParam = parsIt - params.begin();
            assert(iParam == nParams);
         };

         __ecmech_host__
         void getParams(std::vector<double> & params
                        ) const {
            // do not clear params in case adding to an existing set
            int paramsStart = params.size();

            // params.push_back(); // no parameters

            int iParam = params.size() - paramsStart;
            assert(iParam == nParams);
         }

         __ecmech_hdev__ inline const double* getP() const { return _P_ref_vec; };
         __ecmech_hdev__ inline const double* getQ() const { return _Q_ref_vec; };

      private:
         double _P_ref_vec[ ecmech::ntvec * nslip ];
         double _Q_ref_vec[ ecmech::nwvec * nslip ];

   }; // SlipGeomBCC

   /**
    * HCP with <a> slip on basal, prisamtic, and pyramidal families and type-1 <c+a> pyramidal slip
    *
    * the name aBRYcaY1 traces back to EVP_HCP_a_BRY_ca_Y1 (integer code 32) in the old Fortran coding
    *
    * fix me : the coding below is a hack just to get things going ;
    * it is not the best way of doing things, and modifications should be made with great care
    */
   class SlipGeomHCPaBRYcaY1
   {
      public:

         // 3  slip systems in basal <a> family
         // 3  slip systems in prismatic <a> family
         // 6  slip systems in pyramidal <a> family
         // 12  slip systems in pyramidal 1 <c+a> family
         static const int nslip = 3 + 3 + 6 + 12;
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

            int iParam = parsIt - params.begin();
            assert(iParam == nParams);
         };

         __ecmech_host__
         void getParams(std::vector<double> & params
                        ) const {
            // do not clear params in case adding to an existing set
            int paramsStart = params.size();

            params.push_back(_cOverA);

            int iParam = params.size() - paramsStart;
            assert(iParam == nParams);
         }

         __ecmech_hdev__ inline const double* getP() const { return _P_ref_vec; };
         __ecmech_hdev__ inline const double* getQ() const { return _Q_ref_vec; };

      private:
         double _cOverA;
         double _P_ref_vec[ ecmech::ntvec * nslip ];
         double _Q_ref_vec[ ecmech::nwvec * nslip ];
   }; // SlipGeomHCPaBRYcaY1
} // namespace ecmech

#endif // ECMECH_SLIPGEOM_H
