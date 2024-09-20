#pragma once

#include "ECMech_slipgeom_base.h"

namespace ecmech {

   class SlipGeomFCC : public SlipGeom<12>
   {
      public:
         static const bool dynamic = false;
         static constexpr int nParams = 0;

         // constructor and destructor
         SlipGeomFCC() = default;
         __ecmech_hdev__
         ~SlipGeomFCC() {}

         __ecmech_hdev__
         SlipGeomFCC(const double* const params) {
            setParams(params);
         }

         __ecmech_host__
         void setParams(const std::vector<double> & params
                        )
         {
            setParams(params.data());
         }

         __ecmech_hdev__
         void setParams(const double* const)
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
         
   }; // SlipGeomFCC

}