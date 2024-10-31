// -*-c++-*-

#ifndef ECMECH_KINETICS_VOCEPL_H
#define ECMECH_KINETICS_VOCEPL_H

#include <cassert>
#include <cmath>

namespace ecmech {
   /**
    * slip and hardening kinetics
    *
    * power-law slip kinetics with Voce hardening law -- meant to be about as simple as it gets
    */
   template<bool nonlinear>
   class KineticsVocePL
   {
      public:
         static constexpr int nH = 1;
         static constexpr int nParams = 3 + 5 + nH + (nonlinear ? 1 : 0);
         static constexpr int nVals = 1;
         static constexpr int nEvolVals = 2;
         // constructor
         __ecmech_hdev__
         KineticsVocePL(int _nslip) : nslip(_nslip) {}

         // constructor
         __ecmech_hdev__
         KineticsVocePL(const double* const params, int _nslip) :
         nslip(_nslip)
         {
            setParams(params);
         }


         // deconstructor
         ~KineticsVocePL() = default;

         __ecmech_host__
         inline void setParams(const std::vector<double> & params)
         {
            setParams(params.data());
         }

         __ecmech_hdev__
         inline
         void setParams(const double* const params) {
            const double* parsIt = params;

            //////////////////////////////
            // power-law stuff

            m_mu = *parsIt; ++parsIt;
            m_xm = *parsIt; ++parsIt;
            m_gam_w = *parsIt; ++parsIt;

            // CALL fill_power_law(pl)
            // xmm  = xm - one ;
            m_xnn = one / m_xm;
            m_xn = m_xnn - one;
            // xMp1 = xnn + one
            //
            // CALL set_t_min_max(pl)
            m_t_min = pow(ecmech::gam_ratio_min, m_xm);
            m_t_max = pow(ecmech::gam_ratio_ovf, m_xm);

            //////////////////////////////
            // Voce hardening stuff

            m_h0 = *parsIt; ++parsIt;
            m_tausi = *parsIt; ++parsIt;
            m_taus0 = *parsIt; ++parsIt;
            if (nonlinear) {
               m_xmprime = *parsIt; ++parsIt;
               m_xmprime1 = m_xmprime - one;
            }
            else {
               m_xmprime = one;
               m_xmprime1 = zero;
            }
            m_xms = *parsIt; ++parsIt;
            m_gamss0 = *parsIt; ++parsIt;

            //////////////////////////////
            // nH

            m_hdn_init = *parsIt; ++parsIt;

            //////////////////////////////

#if defined(ECMECH_DEBUG)
            int iParam = parsIt - params;
            if (iParam != nParams) {
               ECMECH_FAIL(__func__, "iParam != nParams");
            }
#endif
         }

         __ecmech_host__
         inline void getParams(std::vector<double> & params
                               ) const {
#ifdef ECMECH_DEBUG
            // do not clear params in case adding to an existing set
            int paramsStart = params.size();
#endif

            //////////////////////////////
            // power-law stuff

            params.push_back(m_mu);
            params.push_back(m_xm);
            params.push_back(m_gam_w);

            //////////////////////////////
            // Voce hardening stuff

            params.push_back(m_h0);
            params.push_back(m_tausi);
            params.push_back(m_taus0);
            params.push_back(m_xms);
            params.push_back(m_gamss0);

            //////////////////////////////
            // nH

            params.push_back(m_hdn_init);

            //////////////////////////////

#ifdef ECMECH_DEBUG            
            assert((params.size() - paramsStart) == nParams);
#endif
         }

         __ecmech_host__
         void getHistInfo(std::vector<std::string> & names,
                          std::vector<double>       & init,
                          std::vector<bool>        & plot,
                          std::vector<bool>        & state) const {
            names.push_back("h");
            init.push_back(m_hdn_init);
            plot.push_back(true);
            state.push_back(true);
         }

      private:

         const int nslip; // could template on this if there were call to do so

         // static const _nXnDim = nH*nH ; // do not bother

         //////////////////////////////
         // power-law stuff

         // parameters
         double m_mu; // may evetually set for current conditions
         double m_xm;
         double m_gam_w; // pl%adots, adots0

         // derived from parameters
         double m_t_max, m_t_min, m_xn, m_xnn;

         //////////////////////////////
         // Voce hardening stuff

         double m_h0, m_tausi, m_taus0, m_xms, m_gamss0;
         double m_xmprime, m_xmprime1;

         //////////////////////////////

         double m_hdn_init;

      public:

         __ecmech_hdev__
         inline double getFixedRefRate(const double* const // vals, not used
                                       ) const
         {
            return m_gam_w;
         }

         __ecmech_hdev__
         inline
         double
         getVals(double* const vals,
                 double, // p, not currently used
                 double, // tkelv, not currently used
                 const double* const h_state
                 ) const
         {
            vals[0] = h_state[0]; // _gAll
            assert(vals[0] > zero);
            return vals[0];
         }

         __ecmech_hdev__
         inline
         void
         evalGdots(double* const gdot,
                   double* const dgdot_dtau,
                   const double* const tau,
                   const double* const vals
                   ) const
         {
            double gAll = vals[0]; // gss%h(islip) // _gAll
            for (int iSlip = 0; iSlip<this->nslip; ++iSlip) {
               bool l_act;
               this->evalGdot(gdot[iSlip], l_act, dgdot_dtau[iSlip],
                              gAll,
                              tau[iSlip],
                              m_mu // gss%ctrl%mu(islip)
                              );
            }
         }

         /**
          * see kinetics_pl_d
          */
         __ecmech_hdev__
         inline
         void
         evalGdot(
            double & gdot,
            bool  & l_act,
            double & dgdot_dtau, // wrt resolved shear stress
            double   gIn,
            double   tau,
            double // mu not currently used
            ) const
         {
            // zero things so that can more easily just return in inactive
            //// gdot_w = zero; gdot_r = zero; ! not used by l_linear or l_pl
            gdot = zero;
            //
            dgdot_dtau = zero;
            l_act = false;

            double g_i = one / gIn; // assume have checked gIn>0 elsewhere
            double t_frac = tau * g_i; // has sign of tau
            double at = fabs(t_frac);

            if (at > m_t_min) {
               //
               l_act = true;

               if (at > m_t_max) {
                  // ierr = IERR_OVF_p
                  // set gdot big, evpp may need this for recovery
                  gdot = ecmech::gam_ratio_ovffx * m_gam_w;
                  gdot = copysign(gdot, tau);
                  // do not set any of deriviatives (they are, in truth, zero)
               }
               else {
                  double abslog = log(at);
                  double blog = m_xn * abslog;
                  double temp = m_gam_w * exp(blog);

                  gdot = temp * t_frac;

                  dgdot_dtau = m_xnn * gdot / tau;
                  // dgdot_dtau = temp * m_xnn * g_i; // note: always positive, = xnn * gdot/t
                  // dgdot_dh and dgdot_dg are the same thing for the voce model
               }
            }
         } // evalGdot

         __ecmech_hdev__
         inline
         int
         updateH(double* const hs_u,
                 const double* const hs_o,
                 double dt,
                 const double* const gdot,
                 const double* const /*hvals*/,
                 double tkelv,
                 int outputLevel = 0) const
         {
            double hs_u_1;
            int nFEvals = updateH1<KineticsVocePL>(this,
                                                   hs_u_1, hs_o[0], dt, gdot, tkelv,
                                                   outputLevel);
            hs_u[0] = hs_u_1;

            return nFEvals;
         }

         __ecmech_hdev__
         inline
         void
         getEvolVals(double* const evolVals,
                     const double* const gdot
                     ) const
         {
            // recompute effective shear rate here versus using a stored value
            double shrate_eff = vecsssumabs_n(gdot, nslip); // could switch to template if template class on nslip

            double sv_sat = m_taus0;
            if (shrate_eff > ecmech::idp_tiny_sqrt) {
               sv_sat = m_taus0 * pow((shrate_eff / m_gamss0), m_xms);
            }
            evolVals[0] = shrate_eff;
            evolVals[1] = sv_sat;
         }

         __ecmech_hdev__
         inline
         void
         getSdot1(double &sdot,
                  double &dsdot_ds,
                  double h,
                  const double* const evolVals,
                  double /*tkelv*/
                  ) const
         {
            double shrate_eff = evolVals[0];
            double sv_sat = evolVals[1];
            // When the below ternary op is true then sdot and dsdot_ds remain zero.
            double temp2 = (sv_sat <= m_tausi) ? zero : one / (sv_sat - m_tausi);
            // IF (PRESENT(dfdtkelv)) THEN
            // dfdtkelv(1) = zero
            // END IF

            if (nonlinear) {
               double temp1 = pow((sv_sat - h) * temp2, m_xmprime1);
               sdot = m_h0 * temp1 * (sv_sat - h) * temp2 * shrate_eff;
               dsdot_ds = -m_h0 * temp2 * shrate_eff * m_xmprime * temp1;
            }
            else {
               double temp1 = m_h0 * ((sv_sat - h) * temp2);
               sdot = temp1 * shrate_eff;
               // double dfdshr = temp1 + m_h0 * ( (h - m_tausi) / (temp2*temp2)) * m_xms * sv_sat ;
               dsdot_ds = -m_h0 * temp2 * shrate_eff;
            }
         }
   }; // class KineticsVocePL
} // namespace ecmech

#endif // ECMECH_KINETICS_VOCEPL_H
