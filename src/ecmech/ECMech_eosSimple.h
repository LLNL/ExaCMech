// -*-c++-*-

#ifndef ECMECH_EOS_SIMPLE_H
#define ECMECH_EOS_SIMPLE_H

#include "ECMech_core.h"
#include "ECMech_util.h"

#include <string>
#include <vector>

namespace ecmech {
   template<bool isothermal>
   class EosModelConst
   {
      public:
         static const int nParams = 5;

         // constructor
         __ecmech_hdev__
         EosModelConst() {};

         // deconstructor
         __ecmech_hdev__
         ~EosModelConst() {};


         __ecmech_host__
         inline
         void setParams(const std::vector<double> & params // const double* const params
                        ) {
            std::vector<double>::const_iterator parsIt = params.begin();

            //////////////////////////////

            m_rho0 = *parsIt; ++parsIt;
            m_bulkMod = *parsIt; ++parsIt;
            m_cvav = *parsIt; ++parsIt;
            m_gamma = *parsIt; ++parsIt;
            m_ec0 = *parsIt; ++parsIt;

            m_dtde = one / m_cvav;
            m_tK0 = -m_ec0 * m_dtde;

            //////////////////////////////

            int iParam = parsIt - params.begin();
            if (iParam != nParams) {
               ECMECH_FAIL(__func__, "iParam != nParams");
            }
         };

         __ecmech_host__
         inline
         void getParams(std::vector<double> & params
                        ) const {
            // do not clear params in case adding to an existing set
            int paramsStart = params.size();

            //////////////////////////////

            params.push_back(m_rho0);
            params.push_back(m_bulkMod);
            params.push_back(m_cvav);
            params.push_back(m_gamma);
            params.push_back(m_ec0);

            //////////////////////////////

            int iParam = params.size() - paramsStart;
            if (iParam != nParams) {
               ECMECH_FAIL(__func__, "iParam != nParams");
            }
         };

         __ecmech_hdev__
         inline void evalPT(double &p,
                            double &tK,
                            double  v,
                            double  e) const {
            double mu = one / v - one;

            if (isothermal) {
               p = m_bulkMod * mu;
               tK = m_tK0;
            }
            else {
               p = m_bulkMod * mu + m_gamma * e;
               tK = m_tK0 + e * m_dtde;
            }
         }

         __ecmech_hdev__
         inline
         void evalPTDiff(double &p,
                         double &tK,
                         double &bulkNew,
                         double &dpde,
                         double &dtde,
                         double  v,
                         double  e) const {
            double eta = one / v;
            double mu = eta - one;

            tK = this->evalT(e);

            if (isothermal) {
               p = m_bulkMod * mu;
               dpde = zero;
               dtde = 1e-8 * m_dtde; // instead of zero, to prevent divide-by-zero elsewhere
            }
            else {
               p = m_bulkMod * mu + m_gamma * e;
               dpde = m_gamma;
               dtde = m_dtde;
            }
            bulkNew = m_bulkMod * eta;
         }

         __ecmech_hdev__
         inline
         void getInfo(double &vMin,
                      double &vMax,
                      double &e0,
                      double &v0) const {
            vMin = 0.1;
            vMax = 10.0;
            e0 = 0.0;
            v0 = 1.0;
         }

         __ecmech_hdev__
         inline
         double getBulkRef() const {
            return m_bulkMod;
         }

         __ecmech_hdev__
         inline
         double getRho0() const {
            return m_rho0;
         }

      private:

         __ecmech_hdev__
         inline double evalT(double  e) const {
            double tK;
            if (isothermal) {
               tK = m_tK0;
            }
            else {
               tK = m_tK0 + e * m_dtde;
            }
            return tK;
         }

      private:

         // parameters
         double m_rho0, m_bulkMod, m_gamma, m_ec0, m_cvav;

         // derived from parameters
         double m_dtde, m_tK0;
   }; // class EosModelConst

   template<class EosModel>
   __ecmech_hdev__
   inline
   void updateSimple(const EosModel& eos,
                     double &press,
                     double &tK,
                     double &eNew,
                     double &bulkNew,
                     double &dpde,
                     double &dpdv,
                     double &dtde,
                     double  vNew,
                     double  volInc,
                     double  eOld,
                     double  pOld)
   {

      eNew = eOld - volInc * pOld;

      eos.evalPTDiff(press, tK, bulkNew, dpde, dtde, vNew, eNew);
      dpdv = -bulkNew / vNew;

      double bulkMin = 1e-5 * eos.getBulkRef();
      bulkNew = fmax(bulkMin, bulkNew + dpde * pOld * vNew);
   }
} // namespace ecmech

#endif // ECMECH_EOS_SIMPLE_H
