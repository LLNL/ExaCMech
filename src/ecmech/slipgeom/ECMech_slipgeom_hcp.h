#pragma once

#include "ECMech_slipgeom_base.h"

namespace ecmech {

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
         //static constexpr int nslip = 3 + 3 + 6 + 12;
         static constexpr int nParams = 1;

         // constructor and destructor
         SlipGeomHCPaBRYcaY1() = default;
         __ecmech_hdev__
         ~SlipGeomHCPaBRYcaY1() {}

         __ecmech_hdev__
         SlipGeomHCPaBRYcaY1(const double* const params) {
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

            m_cOverA = *parsIt; ++parsIt;

            // pyramidal 10-11 1-210 depends on c/a
            //
            double m_ya[ecmech::ndim], s_ya[ecmech::ndim];
            {
               double an[ecmech::nMiller] = { one, zero, -one, one }; // plane
               double ab[ecmech::nMiller] = { one, -two, one, zero }; // direction
               //
               miller_to_orthog_sngl(an, ab,
                                     m_ya, s_ya,
                                     m_cOverA);
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
                                     m_cOverA);
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
            params.push_back(m_cOverA);
#ifdef ECMECH_DEBUG
            assert((params.size() - paramsStart) == nParams);
#endif
         }

      private:
         double m_cOverA;
    
   }; // SlipGeomHCPaBRYcaY1

}