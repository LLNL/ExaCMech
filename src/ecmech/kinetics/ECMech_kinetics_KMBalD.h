// -*-c++-*-

#ifndef ECMECH_KINETICS_KMBALD_H
#define ECMECH_KINETICS_KMBALD_H

#include <cassert>
#include <cmath>

#include <string>
#include <vector>

namespace ecmech {
   /**
    * slip and hardening kinetics
    * based on a single Kocks-Mecking dislocation density
    * balanced thermally activated MTS-like slip kinetics with phonon drag effects
    *
    * see \cite{hmx}
    *
    * if withGAthermal then
    *  see subroutine kinetics_mtspwr_d in mdef : (l_mts, l_mtsp, l_plwr)
    *    ! like kinetics_mtswr_d, but with pl%tau_a (possible associated
    *    ! with the Peierls barrier) being the thermally activated part and
    *    ! g being athermal
    *       ! use balanced and pegged MTS model;
    *       ! add to it a low rate sensitivity power law model to take over for high stresses;
    *       ! and combine with drag limited kinetics
    * else then see subroutine kinetics_mtswr_d in mdef
    *
    *   ! note: gdot_w, gdot_r are always positive by definition
    *   !
    *   ! tkelv should only be used for derivative calculations
    *
    * templated on p and q being 1 or not;
    * might eventually template on number of slip systems, but do not do so just yet
    */
   template<bool withGAthermal,
            bool pOne, // l_p_1
            bool qOne, // l_q_1
            bool perSS,
            int  nVPer>
   class KineticsKMBalD
   {
      public:
         static constexpr int nH = 1;
         static constexpr int nParams = 8 + 3 * nVPer + 4 + nH;
         static constexpr int nVals = 2 + nVPer + nVPer;
         static constexpr int nEvolVals = 2;
         // constructor
         __ecmech_hdev__
         KineticsKMBalD(int _nslip) : nslip(_nslip) {
            if (perSS) {
               assert(nslip == nVPer);
            }
            else {
               assert(nVPer == 1);
            }
         }
         // deconstructor
         ~KineticsKMBalD() = default;

         // constructor
         __ecmech_hdev__
         KineticsKMBalD(const double* const params, int _nslip) :
         nslip(_nslip)
         {
            if (perSS) {
               assert(nslip == nVPer);
            }
            else {
               assert(nVPer == 1);
            }
            setParams(params);
         }

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

            m_mu_ref = *parsIt; ++parsIt;
            m_tkelv_ref = *parsIt; ++parsIt;
            for (int iVal = 0; iVal<nVPer; ++iVal) {
               m_c_1[iVal] = *parsIt; ++parsIt;
            }

            m_tau_a = *parsIt; ++parsIt;
            m_p = *parsIt; ++parsIt;
            m_q = *parsIt; ++parsIt;
            m_gam_wo = *parsIt; ++parsIt;
            m_gam_ro = *parsIt; ++parsIt;
            m_wrD = *parsIt; ++parsIt;
            for (int iVal = 0; iVal<nVPer; ++iVal) {
               m_go[iVal] = *parsIt; ++parsIt;
            }

            for (int iVal = 0; iVal<nVPer; ++iVal) {
               m_s[iVal] = *parsIt; ++parsIt;
            }

            if (pOne) {
               assert(m_p == one);
            }
            if (qOne) {
               assert(m_q == one);
            }

            // plaw_from_elawRef
            //
            for (int iVal = 0; iVal<nVPer; ++iVal) {
               // pl%xm = getMtsxmEffective(pl, mu_ref, T_ref)
               double xm = one / (two * ((m_c_1[iVal] / m_tkelv_ref) * m_mu_ref * m_p * m_q));
               //
               // CALL fill_power_law(pl)
               // xmm  = xm - one ;
               m_xnn[iVal] = one / xm;
               m_xn[iVal] = m_xnn[iVal] - one;
               // xMp1 = xnn + one
               //
               // CALL set_t_min_max(pl)
               m_t_min[iVal] = pow(ecmech::gam_ratio_min, xm);
               m_t_max[iVal] = pow(ecmech::gam_ratio_ovf, xm);
            }

            //////////////////////////////
            // Kocks-Mecking stuff

            m_k1 = *parsIt; ++parsIt;
            m_k2o = *parsIt; ++parsIt;
            m_ninv = *parsIt; ++parsIt;
            m_gamma_o = *parsIt; ++parsIt;

            //////////////////////////////
            // nH

            m_hdn_init = *parsIt; ++parsIt;

            m_hdn_min = 1e-4 * m_hdn_init;

            //////////////////////////////

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

            //////////////////////////////
            // power-law stuff

            params.push_back(m_mu_ref);
            params.push_back(m_tkelv_ref);
            for (int iVal = 0; iVal<nVPer; ++iVal) {
               params.push_back(m_c_1[iVal]);
            }

            params.push_back(m_tau_a);
            params.push_back(m_p);
            params.push_back(m_q);
            params.push_back(m_gam_wo);
            params.push_back(m_gam_ro);
            params.push_back(m_wrD);
            for (int iVal = 0; iVal<nVPer; ++iVal) {
               params.push_back(m_go[iVal]);
            }

            for (int iVal = 0; iVal<nVPer; ++iVal) {
               params.push_back(m_s[iVal]);
            }

            //////////////////////////////
            // Kocks-Mecking stuff

            params.push_back(m_k1);
            params.push_back(m_k2o);
            params.push_back(m_ninv);
            params.push_back(m_gamma_o);

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
            names.push_back("rho_dd");
            init.push_back(m_hdn_init);
            plot.push_back(true);
            state.push_back(true);
         }

      private:

         const int nslip; // could template on this if there were call to do so

         //////////////////////////////
         // MTS-like stuff

         // parameters
         double m_mu_ref; // may evetually set for current conditions
         double m_tkelv_ref;
         double m_tau_a; // if withGAthermal then is Peierls barrier
         double m_p; // only used if pOne is false
         double m_q; // only used if qOne is false
         double m_gam_ro;
         double m_gam_wo; // adots0
         double m_c_1[nVPer];
         double m_wrD;
         double m_go[nVPer], m_s[nVPer];

         // derived from parameters
         double m_t_max[nVPer], m_t_min[nVPer], m_xn[nVPer], m_xnn[nVPer];

         //////////////////////////////
         // Kocks-Mecking stuff

         double m_k1, m_k2o, m_ninv, m_gamma_o;

         //////////////////////////////

         double m_hdn_init, m_hdn_min;

      public:

         __ecmech_hdev__
         inline
         double
         getFixedRefRate(const double* const vals) const
         {
            return 1.0 / (1.0 / vals[0] + 1.0 / vals[1]);
         }

         /**
          * @brief Akin to hs_to_gss, power_law_tdep_vals, and plaw_from_hs
          *
          * Could eventually bring in additional pressure and temperature dependence through the dependence of _mu on such ;
          * see use of mu_factors in Fortran code
          */
         __ecmech_hdev__
         inline
         double
         getVals(double* const vals, // [nVals]
                 double, // p, not used
                 double tkelv,
                 const double* const h_state
                 ) const
         {
            double const nVPerInv = 1.0 / nVPer;

            // double sqrtDDens = exp(onehalf * h_state[0]) ; // this is for h_state[0] storing the log of the dislocation density
            double sqrtDDens = sqrt(h_state[0]);

            vals[0] = m_gam_wo / sqrtDDens; // _gam_w
            vals[1] = m_gam_ro * sqrtDDens * sqrtDDens; // _gam_r

            double hdnScale = 0.;
            for (int iVal = 0; iVal<nVPer; ++iVal) {
               double hdnI = m_go[iVal] + m_s[iVal] * sqrtDDens; // _gAll
               hdnScale += hdnI;
               vals[2 + iVal] = hdnI;
               vals[2 + nVPer + iVal] = m_c_1[iVal] / tkelv; // _c_t
               if (!withGAthermal) {
                  assert(vals[2 + iVal] > zero);
               }
            }

            hdnScale = hdnScale * nVPerInv;

            if (withGAthermal) {
               assert(m_tau_a > 0);
            }

            return hdnScale;
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
            for (int iSlip = 0; iSlip<this->nslip; ++iSlip) {
               bool l_act;
               this->evalGdot(gdot[iSlip], l_act, dgdot_dtau[iSlip],
                              vals, iSlip,
                              tau[iSlip],
                              m_mu_ref // gss%ctrl%mu(islip)
                              );
            }
         }

         /**
          * like mts_dG, but with output args first
          */
         __ecmech_hdev__
         inline
         void
         get_mts_dG(double &exp_arg,
                    double &mts_dfac,
                    double c_e, double denom_i, double t_frac) const {
            mts_dfac = c_e * denom_i;

            double p_func;
            if (pOne) {
               p_func = t_frac;
            }
            else {
               if (fabs(t_frac) < idp_tiny_sqrt) {
                  // !! p_dfac is either zero or blows up
                  // !IF (pl%p > one) THEN ! no longer allowed
                  // !   mts_dfac = zero
                  // !ELSE
                  // ! blows up, but just set big
                  p_func = zero;
                  mts_dfac = mts_dfac * 1e10;
                  // !END IF
               }
               else {
                  p_func = pow(fabs(t_frac), m_p);
                  p_func = copysign(p_func, t_frac);
                  mts_dfac = mts_dfac *
                             m_p * p_func / t_frac; // always positive
               }
            }

            double q_arg = one - p_func;
            double pq_fac;
            if (q_arg < idp_tiny_sqrt) {
               // peg
               q_arg = zero;
               mts_dfac = zero;
               pq_fac = zero;
            }
            else {
               if (qOne) {
                  pq_fac = q_arg;
               }
               else {
                  double temp = pow(fabs(q_arg), m_q);
                  mts_dfac = mts_dfac *
                             m_q * temp / fabs(q_arg); // always positive
                  pq_fac = copysign(temp, q_arg);
               }
            }

            exp_arg = -c_e * pq_fac;
         }

         /**
          * see subroutine kinetics_mtspwr_d in mdef
          */
         __ecmech_hdev__
         inline
         void
         evalGdot(
            double & gdot,
            bool   & l_act,
            double & dgdot_dtau, // wrt resolved shear stress
            const double* const vals,
            int      iSlip,
            double   tau,
            double   mu
            ) const
         {
            static const double gdot_w_pl_scaling = 10.0;
            static const double one = 1.0, zero = 0.0;

            const double gam_w = vals[0];
            const double gam_r = vals[1];
            const int iVal = perSS ? iSlip : 0;
            const double gIn = vals[2 + iVal];
            const double c_t = vals[2 + nVPer + iVal];
            const double xn = m_xn[iVal];
            const double xnn = m_xnn[iVal];
            const double t_max = m_t_max[iVal];
            const double t_min = m_t_min[iVal];

            // zero things so that can more easily just return if inactive
            gdot = zero;
            //
            dgdot_dtau = zero;
            l_act = false;

            double g_i;
            double gAth;
            if (withGAthermal) {
               gAth = gIn;
               g_i = one / m_tau_a;
            }
            else {
               gAth = m_tau_a;
               if (tau == zero) {
                  return;
               }
               g_i = one / gIn;
            }
            double at_0 = fmax(zero, fabs(tau) - gAth) * g_i;

            // calculate drag limited kinetics
            //
            double gdot_r, dgdot_r;
            {
               double exp_arg = (fabs(tau) - gAth) / m_wrD;
               double temp;
               if (exp_arg < gam_ratio_min) { // ! IF (gdot_r < gam_ratio_min) THEN
                  // note that this should catch tau <= g
                  return;
               }
               else if (exp_arg < idp_eps_sqrt) {
                  // linear expansion is cheaper and more accurate
                  gdot_r = gam_r * exp_arg;
                  temp = one - exp_arg; // still use temp below as approximation to exp(-fabs(tau)/m_wrD)
               }
               else {
                  temp = exp(-exp_arg);
                  gdot_r = gam_r * (one - temp);
               }
               dgdot_r = gam_r * temp / m_wrD;
            }
            //
            if (at_0 > t_max) {
               // have overflow of thermally activated kinetics, purely drag limited

               gdot = gdot_r;

               dgdot_dtau = dgdot_r;
               gdot = copysign(gdot, tau);

               l_act = true;
               return;
            }

            double gdot_w, dgdot_w;
            double dgdot_wg; // only used if !withGAthermal
            //
            // calculate thermally activated kinetics
            {
               double c_e = c_t * mu;
               //
               double t_frac = (fabs(tau) - gAth) * g_i;
               double exp_arg, mts_dfac;
               get_mts_dG(exp_arg, mts_dfac, c_e, g_i, t_frac);
               //
               if (exp_arg < ln_gam_ratio_min) {
                  // effectively zero due to thermally activated kinetics
                  l_act = false;
                  return;
               }
               //
               //
               // !IF (exp_arg > ln_gam_ratio_ovf) THEN
               // !END IF
               // ! do not need to check the above condition because have pegged the MTS part of the kinetics
               //
               gdot_w = gam_w * exp(exp_arg);
               dgdot_w = mts_dfac * gdot_w;
               if (!withGAthermal) {
                  dgdot_wg = dgdot_w * t_frac;
               }
               //
               double t_frac_m = (-fabs(tau) - gAth) * g_i;
               double exp_arg_m, mts_dfac_m;
               get_mts_dG(exp_arg_m, mts_dfac_m, c_e, g_i, t_frac_m);
               //
               if (exp_arg_m > ln_gam_ratio_min) {
                  // non-vanishing contribution from balancing MTS-like kinetics
                  double gdot_w_m = gam_w * exp(exp_arg_m);
                  gdot_w = gdot_w - gdot_w_m;
                  double contrib = mts_dfac_m * gdot_w_m;
                  dgdot_w = dgdot_w - contrib; // sign used to be the other way, but suspect that was a bug
                  if (!withGAthermal) {
                     dgdot_wg = dgdot_wg - contrib * t_frac_m;
                  }
                  if (fabs(gdot_w / gam_w) < gam_ratio_min) {
                     // effectively zero from roundoff
                     l_act = false;
                     return;
                  }
               }
            }

            if (at_0 > t_min) {
               // need power-law part

               double abslog = log(at_0);
               double blog = xn * abslog;
               double temp = (gam_w * gdot_w_pl_scaling) * exp(blog);

               double gdot_w_pl = temp * at_0; // not signed ! copysign(at_0,tau)
               gdot_w = gdot_w + gdot_w_pl;

               double contrib = temp * xnn * g_i;
               dgdot_w = dgdot_w + contrib;
               if (!withGAthermal) {
                  dgdot_wg = dgdot_wg + contrib * at_0;
               }
            }

            l_act = true;
            //
            {
               gdot = one / (one / gdot_w + one / gdot_r);
               double gdrdiv2 = one / (gdot_r * gdot_r);
               double gdwdiv2 = one / (gdot_w * gdot_w);
               dgdot_dtau = (gdot * gdot) * (dgdot_w * gdwdiv2 + dgdot_r * gdrdiv2);
            }

            gdot = copysign(gdot, tau);
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
            // do not yet both with l_overdriven and setting-to-saturation machinery as in Fortran coding

            // update is done on log(h) -- h treated as a nomralized (unitless) dislocation density
            double log_hs_u;
            double log_hs_o = log(fmax(hs_o[0], m_hdn_min));
            int nFEvals = updateH1<KineticsKMBalD>(this,
                                                   log_hs_u, log_hs_o, dt, gdot, tkelv,
                                                   outputLevel);
            hs_u[0] = exp(log_hs_u);

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

            double k2 = m_k2o;
            if (shrate_eff > ecmech::idp_tiny_sqrt) {
               k2 = m_k2o * pow((m_gamma_o / shrate_eff), m_ninv);
            }

            evolVals[0] = shrate_eff;
            evolVals[1] = k2;
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
            double k2 = evolVals[1];
            double temp_hs_a = exp(-onehalf * h);
            double temp1 = m_k1 * temp_hs_a - k2;
            sdot = temp1 * shrate_eff;
            dsdot_ds = (-m_k1 * onehalf * temp_hs_a) * shrate_eff;
            // }
         }
   }; // class KineticsKMBalD
} // namespace ecmech

#endif // ECMECH_KINETICS_KMBALD_H
