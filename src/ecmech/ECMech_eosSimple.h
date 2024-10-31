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
         static constexpr int nParams = 5;

         // constructor
         EosModelConst() = default;
         __ecmech_hdev__
         EosModelConst(const double* const params) {
            setParams(params);
         }

         // deconstructor
         ~EosModelConst() = default;

         __ecmech_host__
         inline
         void setParams(const std::vector<double>& params) {
            setParams(params.data());
         }

         __ecmech_hdev__
         inline
         void setParams(const double* const params) {
            const double* parsIt = params;
            //////////////////////////////
            m_density0 = *parsIt; ++parsIt;
            m_bulk_modulus = *parsIt; ++parsIt;
            m_cvav = *parsIt; ++parsIt;
            m_gamma = *parsIt; ++parsIt;
            m_cold_energy0 = *parsIt; ++parsIt;

            m_dtde = one / m_cvav;
            m_tkelv0 = -m_cold_energy0 * m_dtde;

            //////////////////////////////
#if defined(ECMECH_DEBUG)
            int iParam = parsIt - params;
            if (iParam != nParams) {
               ECMECH_FAIL(__func__, "iParam != nParams");
            }
#endif
         }

         __ecmech_host__
         inline
         void getParams(std::vector<double> & params
                        ) const {
            // do not clear params in case adding to an existing set
            int paramsStart = params.size();

            //////////////////////////////

            params.push_back(m_density0);
            params.push_back(m_bulk_modulus);
            params.push_back(m_cvav);
            params.push_back(m_gamma);
            params.push_back(m_cold_energy0);

            //////////////////////////////

            int iParam = params.size() - paramsStart;
            if (iParam != nParams) {
               ECMECH_FAIL(__func__, "iParam != nParams");
            }
         }

         __ecmech_hdev__
         inline void evalPT(double &pressure,
                            double &tkelv,
                            double  rel_vol,
                            double  energy) const {
            double mu = one / rel_vol - one;

            if (isothermal) {
               pressure = m_bulk_modulus * mu;
               tkelv = m_tkelv0;
            }
            else {
               pressure = m_bulk_modulus * mu + m_gamma * energy;
               tkelv = m_tkelv0 + energy * m_dtde;
            }
         }

         __ecmech_hdev__
         inline
         void evalPTDiff(double &pressure,
                         double &tkelv,
                         double &bulk_modulus_new,
                         double &dpde,
                         double &dtde,
                         double  rel_vol,
                         double  energy) const {
            double eta = one / rel_vol;
            double mu = eta - one;

            tkelv = this->evalT(energy);

            if (isothermal) {
               pressure = m_bulk_modulus * mu;
               dpde = zero;
               dtde = 1e-8 * m_dtde; // instead of zero, to prevent divide-by-zero elsewhere
            }
            else {
               pressure = m_bulk_modulus * mu + m_gamma * energy;
               dpde = m_gamma;
               dtde = m_dtde;
            }
            bulk_modulus_new = m_bulk_modulus * eta;
         }

         __ecmech_hdev__
         inline
         void getInfo(double &rel_vol_min,
                      double &rel_vol_max,
                      double &energy0,
                      double &rel_vol0) const {
            rel_vol_min = 0.1;
            rel_vol_max = 10.0;
            energy0 = 0.0;
            rel_vol0 = 1.0;
         }

         __ecmech_hdev__
         inline
         double getBulkRef() const {
            return m_bulk_modulus;
         }

         __ecmech_hdev__
         inline
         double getRho0() const {
            return m_density0;
         }

      private:

         __ecmech_hdev__
         inline double evalT(double  energy) const {
            double tkelv;
            if (isothermal) {
               tkelv = m_tkelv0;
            }
            else {
               tkelv = m_tkelv0 + energy * m_dtde;
            }
            return tkelv;
         }

      private:

         // parameters
         double m_density0, m_bulk_modulus, m_gamma, m_cold_energy0, m_cvav;

         // derived from parameters
         double m_dtde, m_tkelv0;
   }; // class EosModelConst

   template<class EosModel>
   __ecmech_hdev__
   inline
   void updateSimple(const EosModel& eos,
                     double &press,
                     double &tkelv,
                     double &energy_new,
                     double &bulk_modulus_new,
                     double &dpde,
                     double &dpdv,
                     double &dtde,
                     double  rel_vol_new,
                     double  rel_vol_increment,
                     double  energy_old,
                     double  pressure_old)
   {

      energy_new = energy_old - rel_vol_increment * pressure_old;

      eos.evalPTDiff(press, tkelv, bulk_modulus_new, dpde, dtde, rel_vol_new, energy_new);
      dpdv = -bulk_modulus_new / rel_vol_new;

      double bulk_modulus_min = 1e-5 * eos.getBulkRef();
      bulk_modulus_new = fmax(bulk_modulus_min, bulk_modulus_new + dpde * pressure_old * rel_vol_new);
   }
} // namespace ecmech

#endif // ECMECH_EOS_SIMPLE_H
