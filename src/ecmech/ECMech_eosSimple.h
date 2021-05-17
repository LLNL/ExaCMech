// -*-c++-*-

#ifndef ECMECH_EOS_SIMPLE_H
#define ECMECH_EOS_SIMPLE_H

#include "ECMech_core.h"

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

            _rho0 = *parsIt; ++parsIt;
            _bulkMod = *parsIt; ++parsIt;
            _cvav = *parsIt; ++parsIt;
            _gamma = *parsIt; ++parsIt;
            _ec0 = *parsIt; ++parsIt;

            _dtde = one / _cvav;
            _tK0 = -_ec0 * _dtde;

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

            params.push_back(_rho0);
            params.push_back(_bulkMod);
            params.push_back(_cvav);
            params.push_back(_gamma);
            params.push_back(_ec0);

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
               p = _bulkMod * mu;
               tK = _tK0;
            }
            else {
               p = _bulkMod * mu + _gamma * e;
               tK = _tK0 + e * _dtde;
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

            if (isothermal) {
               p = _bulkMod * mu;
               tK = _tK0;
               dpde = zero;
               dtde = 1e-8 * _dtde; // instead of zero, to prevent divide-by-zero elsewhere
            }
            else {
               p = _bulkMod * mu + _gamma * e;
               tK = _tK0 + e * _dtde;
               dpde = _gamma;
               dtde = _dtde;
            }
            bulkNew = _bulkMod * eta;
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
            return _bulkMod;
         }

         __ecmech_hdev__
         inline
         double getRho0() const {
            return _rho0;
         }

      private:

         // parameters
         double _rho0, _bulkMod, _gamma, _ec0, _cvav;

         // derived from parameters
         double _dtde, _tK0;
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
      // double vOld = volRatio[0] ; // not needed

      eNew = eOld - volInc * pOld;

      eos.evalPTDiff(press, tK, bulkNew, dpde, dtde, vNew, eNew);
      dpdv = -bulkNew / vNew;

      bulkNew = bulkNew + dpde * pOld * vNew;
   }
} // namespace ecmech

#endif // ECMECH_EOS_SIMPLE_H
