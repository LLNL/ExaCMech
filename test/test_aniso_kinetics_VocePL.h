// -*-c++-*-

/**
 * @file test_aniso_kinetics_VocePL.h
 *
 * @brief Test-only companion to `kinetics/ECMech_kinetics_VocePL.h`: `KineticsAnisoVocePL`,
 * a per-slip-system-hardening variant of the linear/nonlinear-Voce power-law kinetics
 * model, used exclusively by `test_aniso_hardening.cxx` to exercise `updateHN`'s
 * vector-hardening solve path (`kinetics/ECMech_kinetics.h`) with a real, if simplified,
 * multi-component hardening law -- the production `KineticsVocePL` only ever has a
 * single (`nH == 1`) hardening state.
 *
 * The physics is otherwise identical to `KineticsVocePL`: the same power-law slip-rate
 * law (`evalGdot`) and the same per-slip-system Voce hardening-rate form (`getSdotN`),
 * just with `nH == Nslip` independent CRSS states (`m_tausi[]`, one per slip system)
 * instead of one shared value. **Simplification**: each slip system's hardening ODE only
 * depends on its own CRSS -- `getSdotN`'s Jacobian `dsdot_ds` is purely diagonal (no
 * cross-slip-system hardening coupling terms), which keeps this test class simple but
 * means it isn't a template for how a "real" anisotropic hardening model's Jacobian
 * should look.
 */

#ifndef TEST_ANISO_KINETICS_VOCEPL_H
#define TEST_ANISO_KINETICS_VOCEPL_H

#include <cassert>
#include <cmath>
#include "ECMech_port.h"

/** @brief Row-major flattening of a `(p, q)` index into an `nDim × nDim` matrix; shared with the production kinetics headers. */
#define ECMECH_NN_INDX(p, q, nDim) (p) * (nDim) + (q)

namespace ecmech {
   /**
    * @brief Linear/nonlinear-Voce power-law kinetics with one independent hardening
    * state per slip system, rather than `KineticsVocePL`'s single shared state.
    *
    * @tparam nonlinear Selects the linear (`false`) or nonlinear (`true`) Voce
    * saturation-stress rate-sensitivity form, exactly as in `KineticsVocePL`.
    * @tparam Nslip Number of slip systems, and (since `nH == Nslip` here) the number of
    * independent hardening states.
    */
   template<bool nonlinear,
            int Nslip>
   class KineticsAnisoVocePL
   {
      public:
         /** @brief One hardening state per slip system (unlike `KineticsVocePL`'s single shared state). */
         static constexpr int nH = Nslip;
         static constexpr int nslip = Nslip;
         /** @brief `mu`, `xm`, `gam_w` (power-law) + `h0`, `tausi[Nslip]`, `taus0`, [`xmprime` if `nonlinear`], `xms`, `gamss0` (Voce) + `hdn_init` -- see `setParams`. */
         static constexpr int nParams = 3 + 5 + nH + (nonlinear ? 1 : 0);
         static constexpr int nVals = nslip;
         static constexpr int nEvolVals = 2;

         // constructor
         __ecmech_hdev__
         KineticsAnisoVocePL(int /* nslip_ */) {};
         // deconstructor
         __ecmech_hdev__
         ~KineticsAnisoVocePL() {}

         /**
          * @brief Set parameters from a flat array; order matches `KineticsVocePL::setParams`
          * except `tausi` is `nslip` values (one initial CRSS per slip system) rather than one.
          * @param params Flat parameter array of length #nParams.
          */
         __ecmech_host__
         inline void setParams(const std::vector<double> & params // const double* const params
                               ) {
            std::vector<double>::const_iterator parsIt = params.begin();

            //////////////////////////////
            // power-law stuff

            m_shear_modulus = *parsIt; ++parsIt;
            m_xm = *parsIt; ++parsIt;
            m_gam_w = *parsIt; ++parsIt;

            // CALL fill_power_law(pl)
            // xmm  = xm - one ;
            m_xnn = one / m_xm;
            m_xn = m_xnn - one;
            // xMp1 = xnn + one
            //
            // CALL setm_t_min_max(pl)
            m_t_min = pow(ecmech::gam_ratio_min, m_xm);
            m_t_max = pow(ecmech::gam_ratio_ovf, m_xm);

            //////////////////////////////
            // Voce hardening stuff

            m_h0 = *parsIt; ++parsIt;
            for (int i = 0; i < nslip; i++) {
               m_tausi[i] = *parsIt; ++parsIt;
            }

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

            assert((parsIt - params.begin()) == nParams);
         }

         /** @brief Inverse of `setParams`: appends this instance's current parameters (in the same order) onto `params`. */
         __ecmech_host__
         inline void getParams(std::vector<double> & params
                               ) const {
#ifdef ECMECH_DEBUG
            // do not clear params in case adding to an existing set
            int paramsStart = params.size();
#endif

            //////////////////////////////
            // power-law stuff

            params.push_back(m_shear_modulus);
            params.push_back(m_xm);
            params.push_back(m_gam_w);

            //////////////////////////////
            // Voce hardening stuff

            params.push_back(m_h0);
            for (int i = 0; i < nslip; i++) {
               params.push_back(m_tausi[i]);
            }

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

         /**
          * @brief Report a single representative initial hardening value, *not* one per
          * slip system.
          * @note Deviates from the usual kinetics-class contract of reporting `nH`
          * history entries: this always pushes exactly one `"h"` entry (`m_hdn_init`)
          * regardless of `nH == nslip`. `test_aniso_hardening.cxx` compensates by
          * broadcasting this single value across all `nslip` initial-state entries
          * itself (`std::fill`) rather than relying on `getHistInfo` to supply them.
          */
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

         // static const _nXnDim = nH*nH ; // do not bother

         //////////////////////////////
         // power-law stuff

         // parameters
         /** @brief Reference shear modulus for the power-law slip rate. */
         double m_shear_modulus; // may evetually set for current conditions
         /** @brief Power-law rate-sensitivity exponent (1/m). */
         double m_xm;
         /** @brief Reference/normalizing shear rate for the power-law slip rate. */
         double m_gam_w; // pl%adots, adots0

         // derived from parameters
         /** @brief Overflow/underflow stress-ratio thresholds and power-law exponent helpers derived from `m_xm`. */
         double m_t_max, m_t_min, m_xn, m_xnn;

         //////////////////////////////
         // Voce hardening stuff

         /** @brief Voce hardening-rate coefficient, saturation-stress reference value, and saturation-stress rate-sensitivity exponent/reference rate (shared across all slip systems). */
         double m_h0, m_taus0, m_xms, m_gamss0;
         /** @brief Per-slip-system initial CRSS -- the one place this class differs structurally from `KineticsVocePL`. */
         double m_tausi[nslip];
         /** @brief Nonlinear-Voce exponent and its (exponent - 1) precompute; fixed at `1`/`0` (a no-op multiplicatively) when `nonlinear == false`. */
         double m_xmprime, m_xmprime1;

         //////////////////////////////

         /** @brief Single shared initial value for every slip system's hardening state (see `getHistInfo`'s `@note`). */
         double m_hdn_init;

      public:

         /** @brief Reference shear rate used to non-dimensionalize slip rates elsewhere in the solver machinery. */
         __ecmech_hdev__
         inline double getFixedRefRate(const double* const // vals, not used
                                       ) const
         {
            return m_gam_w;
         }

         /**
          * @brief Copy the per-slip-system CRSS directly out of the hardening state
          * (no forest-hardening combination step, unlike the production kinetics
          * classes) and return their mean.
          * @param[out] vals Per-slip-system CRSS, copied from `h_state`.
          * @param[in] h_state Current hardening state (one value per slip system).
          * @return Mean CRSS across all slip systems.
          */
         __ecmech_hdev__
         inline
         double
         getVals(double* const vals,
                 double, // p, not currently used
                 double, // tkelv, not currently used
                 const double* const h_state
                 ) const
         {
            double mVals = ecmech::zero;
            for (int iSlip = 0; iSlip < nslip; ++iSlip) {
               vals[iSlip] = h_state[iSlip]; // _gAll
               mVals += vals[iSlip];
               assert(vals[iSlip] > zero);
            }

            mVals /= nslip;

            return mVals;
         }

         /** @brief Evaluate the power-law slip rate (and its derivatives) on every slip system independently, using each system's own CRSS from `vals`. */
         __ecmech_hdev__
         inline
         void
         evalGdots(double* const gdot,
                   double* const dgdot_dtau,
                   double* const dgdot_dg,
                   const double* const tau,
                   const double* const vals
                   ) const
         {
            for (int iSlip = 0; iSlip<this->nslip; ++iSlip) {
               bool l_act;
               double gAll = vals[iSlip];
               this->evalGdot(gdot[iSlip], l_act, dgdot_dtau[iSlip], dgdot_dg[iSlip],
                              gAll,
                              tau[iSlip],
                              m_shear_modulus // gss%ctrl%shear_modulus(islip)
                              );
            }
         }

         /**
          * @brief Single-slip-system power-law slip rate: `γ̇ = gam_w · sign(τ) · |τ/g|^(1/xm)`,
          * inactive (zero rate/derivatives) below the `m_t_min` rate-independent floor and
          * clamped to a large finite value above the `m_t_max` overflow guard -- same
          * form as `KineticsVocePL::evalGdot`.
          * @param[out] gdot Slip rate.
          * @param[out] l_act Whether this slip system is above the rate-independent floor.
          * @param[out] dgdot_dtau Derivative of `gdot` with respect to resolved shear stress.
          * @param[out] dgdot_dg Derivative of `gdot` with respect to slip-system strength.
          * @param[in] gIn Current CRSS for this slip system.
          * @param[in] tau Resolved shear stress on this slip system.
          */
         __ecmech_hdev__
         inline
         void
         evalGdot(
            double & gdot,
            bool  & l_act,
            double & dgdot_dtau, // wrt resolved shear stress
            double & dgdot_dg, // wrt slip system strength
            double   gIn,
            double   tau,
            double // shear_modulus, not currently used
            ) const
         {
            // zero things so that can more easily just return in inactive
            //// gdot_w = zero; gdot_r = zero; ! not used by l_linear or l_pl
            gdot = zero;
            //
            dgdot_dtau = zero;
            dgdot_dg = zero;
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

                  dgdot_dtau = temp * m_xnn * g_i; // note: always positive, = xnn * gdot/t
                  dgdot_dg = -dgdot_dtau * t_frac; // = - gdot * xnn * g_i
               }
            }
         } // evalGdot

         /**
          * @brief Solve the `nslip`-component (vector) hardening-state update via
          * `updateHN` (`kinetics/ECMech_kinetics.h`)'s trust-region dogleg solver.
          * @param[out] hs_u Updated per-slip-system hardening state.
          * @param[in] hs_o Beginning-of-step per-slip-system hardening state.
          * @param[in] dt Time-step size.
          * @param[in] gdot Slip rates driving the hardening update.
          * @return Number of solver function evaluations.
          */
         __ecmech_hdev__
         inline
         int
         updateH(double* const hs_u,
                 const double* const hs_o,
                 double dt,
                 const double* const gdot,
                 const double* const /*hvals*/,
                 double /*tkelv*/,
                 int outputLevel = 0) const
         {
            double hs_u_1[nslip];
            int nFEvals = updateHN<KineticsAnisoVocePL>(this,
                                                        &hs_u_1[0], &hs_o[0], dt, gdot, nullptr, 0.0,
                                                        outputLevel);

            for (int i = 0; i < nslip; i++) {
               hs_u[i] = hs_u_1[i];
            }

            return nFEvals;
         }

         /**
          * @brief Precompute the total effective shear rate and the (rate-dependent)
          * saturation stress shared by every slip system's hardening ODE this step.
          * @param[out] evolVals `[0]`: total effective shear rate `Σ|γ̇ᵢ|`; `[1]`: Voce
          * saturation stress at that rate.
          * @param[in] gdot Per-slip-system slip rates.
          */
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

         /**
          * @brief Per-slip-system Voce hardening rate `ḣᵢ = h0 · (1 - hᵢ/sv_sat)^xmprime ·
          * shrate_eff` (or its linear form when `nonlinear == false`) and its Jacobian.
          * @note The Jacobian `dsdot_ds` is purely diagonal -- slip system `i`'s hardening
          * rate depends only on its own state `h[i]`, not on any other system's (see this
          * file's `@file` doc for why that's a deliberate test-only simplification).
          * @param[out] sdot Per-slip-system hardening rate.
          * @param[out] dsdot_ds `nslip × nslip` Jacobian of `sdot` with respect to `h` (only the diagonal is nonzero).
          * @param[in] h Current per-slip-system hardening state.
          * @param[in] evolVals `[0]`/`[1]` from `getEvolVals` (effective shear rate, saturation stress).
          */
         __ecmech_hdev__
         inline
         void
         getSdotN(double *sdot,
                  double *dsdot_ds,
                  const double* const h,
                  const double* const evolVals,
                  const double* const /*hvals*/,
                  double /*tkelv*/,
                  double* const dsdot_dgdot = nullptr) const
         {
            double shrate_eff = evolVals[0];
            double sv_sat = evolVals[1];

            for (int i = 0; i < nslip * nslip; i++) {
               dsdot_ds[i] = 0.0;
            }

            for (int iSlip = 0; iSlip < nslip; iSlip++) {
               // When the below ternary op is true then sdot and dsdot_ds remain zero.
               double temp2 = (sv_sat <= m_tausi[iSlip]) ? zero : one / (sv_sat - m_tausi[iSlip]);

               // IF (PRESENT(dfdtkelv)) THEN
               // dfdtkelv(1) = zero
               // END IF
               // Just throwing this in here for the tests
               // in reality we would need to set dsdot_ds in another section
               // after checking if it's a nullptr or not
               assert(dsdot_ds != nullptr);

               if (nonlinear) {
                  double temp1 = pow((sv_sat - h[iSlip]) * temp2, m_xmprime1);
                  sdot[iSlip] = m_h0 * temp1 * (sv_sat - h[iSlip]) * temp2 * shrate_eff;
                  dsdot_ds[ECMECH_NN_INDX(iSlip, iSlip, nslip)] = -m_h0 * temp2 * shrate_eff * m_xmprime * temp1;
               }
               else {
                  double temp1 = m_h0 * ((sv_sat - h[iSlip]) * temp2);
                  sdot[iSlip] = temp1 * shrate_eff;
                  // double dfdshr = temp1 + m_h0 * ( (h - m_tausi) / (temp2*temp2)) * m_xms * sv_sat ;
                  dsdot_ds[ECMECH_NN_INDX(iSlip, iSlip, nslip)] = -m_h0 * temp2 * shrate_eff;
               }
            }
         if (dsdot_dgdot) {
            ECMECH_FAIL("test_aniso", "This model does not implement the dsdot_dgdot feature");
         }
         }
   }; // class KineticsVocePL
} // namespace ecmech

#endif // ECMECH_KINETICS_VOCEPL_H
