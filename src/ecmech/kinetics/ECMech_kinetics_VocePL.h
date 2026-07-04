/**
 * @file ECMech_kinetics_VocePL.h
 * @brief Power-law slip kinetics with (optionally nonlinear) Voce isotropic hardening --
 * the simplest kinetics model in ExaCMech.
 *
 * **Slip-rate law** (see KineticsVocePL::evalGdot): a single isotropic power law shared
 * by every slip system,
 *
 *   γ̇ = γ̇_w · (|τ| / g)^(1/m) · sign(τ)
 *
 * where γ̇ = `gdot`, γ̇_w = `m_gam_w` (a reference slip rate), τ = `tau` (resolved shear
 * stress), g = `gIn` (the current CRSS, from the hardening state, shared by all slip
 * systems), and m = `m_xm` (the rate-sensitivity exponent).
 *
 * Below a small stress-ratio threshold the slip rate is treated as exactly zero; above a
 * large threshold it is capped rather than allowed to overflow (see #gam_ratio_min /
 * #gam_ratio_ovf in ECMech_const.h).
 *
 * **Hardening law** (see KineticsVocePL::getSdot1): a Voce-type saturation law,
 *
 *   dh/dγ = h₀ · γ̇_eff · [(τ_s - h) / (τ_s - τ_si)]^m′
 *   τ_s = τ_s0 · (γ̇_eff / γ̇_s0)^x_s
 *
 * where h = the current CRSS (hardening state), h₀ = `m_h0` (initial hardening rate),
 * γ̇_eff = `shrate_eff` (effective, sum-of-absolute-value, shear rate across all slip
 * systems), τ_s = `sv_sat` (rate-dependent saturation stress), τ_si = `m_tausi` (initial,
 * lower-bound, CRSS), and m′ = `m_xmprime` (fixed at 1 when the `nonlinear` template
 * parameter is `false`, recovering the classic linear Voce law; when `true`, this extra
 * parameter generalizes the shape of the hardening curve as it approaches saturation).
 * τ_s0 = `m_taus0`, γ̇_s0 = `m_gamss0`, and x_s = `m_xms` together define how the
 * saturation stress itself scales with the current effective shear rate.
 *
 * @see ECMech_kinetics.h for the kinetics model interface contract this class
 * implements, and KineticsKMBalD/KineticsOrowanD/KineticsBCCMD for the other available
 * kinetics models
 */

// -*-c++-*-

#ifndef ECMECH_KINETICS_VOCEPL_H
#define ECMECH_KINETICS_VOCEPL_H

#include <cassert>
#include <cmath>

namespace ecmech {
   /**
    * @brief Power-law slip kinetics with (optionally nonlinear) Voce isotropic hardening.
    *
    * All slip systems share a single, isotropic critical resolved shear stress (CRSS);
    * hardening follows a Voce-type saturation law. See the file-level documentation in
    * ECMech_kinetics_VocePL.h for the governing equations.
    *
    * @tparam nonlinear If `false`, the Voce hardening-rate exponent is fixed at 1
    * (classic linear Voce law). If `true`, an extra parameter (`m_xmprime`) generalizes
    * the exponent, allowing a non-linear approach to saturation.
    *
    * @ingroup ECMech_kinetics
    * @see ECMech_kinetics.h for the required-traits/required-member-functions interface
    * every kinetics model must satisfy
    */
   template<bool nonlinear>
   class KineticsVocePL
   {
      public:
         /** @brief Number of hardening state variables: 1 (a single isotropic CRSS shared by all slip systems). */
         static constexpr int nH = 1;
         /** @brief Number of parameters: 3 power-law + 5 Voce hardening + nH initial-state + (1 more if `nonlinear`). */
         static constexpr int nParams = 3 + 5 + nH + (nonlinear ? 1 : 0);
         /** @brief Number of kinetic values precomputed by getVals and reused by evalGdots: 1 (the shared CRSS). */
         static constexpr int nVals = 1;
         /** @brief Number of intermediate values precomputed by getEvolVals and reused by getSdot1: 2 (effective shear rate and saturation stress). */
         static constexpr int nEvolVals = 2;
         /**
          * @brief Construct with a given number of slip systems; parameters must be set
          * separately via setParams.
          * @param _nslip Number of slip systems.
          */
         __ecmech_hdev__
         KineticsVocePL(int _nslip) : nslip(_nslip) {}

         /**
          * @brief Construct with a given number of slip systems and immediately set
          * parameters.
          * @param params Parameter array; see setParams for the expected order.
          * @param _nslip Number of slip systems.
          */
         __ecmech_hdev__
         KineticsVocePL(const double* const params, int _nslip) :
         nslip(_nslip)
         {
            setParams(params);
         }


         /** @brief Destructor (default; no owned resources). */
         ~KineticsVocePL() = default;

         /**
          * @brief Set parameters from a `std::vector` (host-side convenience wrapper
          * around the array-based overload).
          * @param params Parameter vector; see the array overload for the expected order.
          */
         __ecmech_host__
         inline void setParams(const std::vector<double> & params)
         {
            setParams(params.data());
         }

         /**
          * @brief Set parameters from a flat array.
          *
          * Expected order: `mu`, `xm`, `gam_w` (power-law), then `h0`, `tausi`, `taus0`,
          * [`xmprime` if `nonlinear`], `xms`, `gamss0` (Voce hardening), then
          * `hdn_init` (initial hardening state). Also derives the power-law exponent
          * helpers (`m_xnn`, `m_xn`) and the overflow/underflow stress-ratio thresholds
          * (`m_t_min`, `m_t_max`).
          * @param params Flat parameter array of length #nParams.
          */
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

         /**
          * @brief Append this model's current parameters to `params`, in the same order
          * setParams expects them.
          * @param[in,out] params Parameter vector to append to (not cleared first).
          */
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

         /**
          * @brief Describe the single hardening history variable ("h", the isotropic
          * CRSS) for history-array bookkeeping.
          * @param[out] names Appended with `"h"`.
          * @param[out] init Appended with #m_hdn_init.
          * @param[out] plot Appended with `true`.
          * @param[out] state Appended with `true`.
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

         /** @brief Number of slip systems handled by this instance. */
         const int nslip; // could template on this if there were call to do so

         // static const _nXnDim = nH*nH ; // do not bother

         //////////////////////////////
         // power-law stuff

         // parameters
         /** @brief Shear modulus [stress units]; not currently used by evalGdot, reserved for future pressure/temperature-dependent kinetics. */
         double m_mu; // may evetually set for current conditions
         /** @brief Power-law rate-sensitivity exponent `m`. */
         double m_xm;
         /** @brief Reference/characteristic slip rate `gam_w` [1/time]; also returned by getFixedRefRate. */
         double m_gam_w; // pl%adots, adots0

         // derived from parameters
         /** @brief Derived overflow/underflow thresholds on the stress ratio `|tau|/g` (see #gam_ratio_min, #gam_ratio_ovf) and power-law exponent helpers `xnn = 1/xm`, `xn = xnn - 1`. */
         double m_t_max, m_t_min, m_xn, m_xnn;

         //////////////////////////////
         // Voce hardening stuff

         /** @brief Voce hardening parameters: initial hardening rate `h0`, initial (lower-bound) CRSS `tausi`, reference saturation stress `taus0` [stress units], saturation-stress rate-sensitivity exponent `xms`, reference shear rate `gamss0` [1/time] at which `taus0` is defined. */
         double m_h0, m_tausi, m_taus0, m_xms, m_gamss0;
         /** @brief Nonlinear Voce shape exponent `xmprime` and its `xmprime - 1` (only meaningful when `nonlinear` is true; fixed at 1/0 otherwise). */
         double m_xmprime, m_xmprime1;

         //////////////////////////////

         /** @brief Initial value of the (single) hardening state variable `h`. */
         double m_hdn_init;

      public:

         /**
          * @brief Reference slip rate used for scaling elsewhere in the solve (e.g. by
          * the evptn elastic-strain/rotation solver).
          * @return #m_gam_w.
          */
         __ecmech_hdev__
         inline double getFixedRefRate(const double* const // vals, not used
                                       ) const
         {
            return m_gam_w;
         }

         /**
          * @brief Precompute the kinetic values used by evalGdots: for this isotropic
          * model, simply the current CRSS.
          * @param[out] vals `vals[0]` set to the current CRSS, `h_state[0]`.
          * @param p Pressure; not currently used by this model.
          * @param tkelv Temperature [Kelvin]; not currently used by this model.
          * @param[in] h_state Current hardening state, `h_state[0]` the CRSS.
          * @return `vals[0]`, the current CRSS.
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
            vals[0] = h_state[0]; // _gAll
            assert(vals[0] > zero);
            return vals[0];
         }

         /**
          * @brief Evaluate the slip rate and its derivative w.r.t. resolved shear stress
          * on every slip system, using the shared isotropic CRSS from getVals.
          * @param[out] gdot Slip rate on each slip system [1/time].
          * @param[out] dgdot_dtau Derivative of slip rate w.r.t. resolved shear stress on
          * each slip system.
          * @param[in] tau Resolved shear stress on each slip system [stress units].
          * @param[in] vals Kinetic values from getVals; `vals[0]` is the shared CRSS.
          */
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
          * @brief Evaluate the power-law slip rate and its derivative on a single slip
          * system:
          *
          *   γ̇ = γ̇_w · (|τ| / g)^(1/m) · sign(τ)
          *   ∂γ̇/∂τ = (1/m) · (γ̇/τ)
          *
          * where γ̇ = `gdot`, γ̇_w = `m_gam_w`, τ = `tau`, g = `gIn`, and m = `m_xm`.
          *
          * The stress ratio |τ/g| (`at`) is compared against thresholds derived from
          * #gam_ratio_min/#gam_ratio_ovf (via `m_t_min`/`m_t_max`): below `m_t_min` the
          * system is treated as inactive (`gdot = 0`); above `m_t_max` the rate is capped
          * at `gam_ratio_ovffx * gam_w` rather than allowed to overflow, with derivatives
          * left at zero.
          * @param[out] gdot Slip rate γ̇ [1/time].
          * @param[out] l_act `true` if the slip system is active (stress ratio above
          * `m_t_min`).
          * @param[out] dgdot_dtau Derivative of `gdot` w.r.t. resolved shear stress.
          * @param gIn Current CRSS (flow strength) for this slip system, g.
          * @param tau Resolved shear stress τ [stress units].
          * @param mu Shear modulus; not currently used.
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

         /**
          * @brief Advance the single CRSS hardening state variable one time step by
          * delegating to the shared scalar-hardness SNLS solve.
          * @param[out] hs_u End-of-step hardening state, `hs_u[0]`.
          * @param[in] hs_o Start-of-step hardening state, `hs_o[0]`.
          * @param dt Time step size.
          * @param[in] gdot Slip rates on all slip systems, used to compute the evolution
          * inputs via getEvolVals.
          * @param hvals Unused by this model.
          * @param tkelv Temperature [Kelvin].
          * @param outputLevel Verbosity passed through to the SNLS solver.
          * @return Function-evaluation count from updateH1, or a negative value if the
          * solve failed to converge.
          * @see updateH1 in ECMech_kinetics.h
          */
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

         /**
          * @brief Precompute the effective shear rate γ̇_eff and rate-dependent
          * saturation stress τ_s used by getSdot1 (see the file-level documentation in
          * ECMech_kinetics_VocePL.h for both formulas).
          * @param[out] evolVals `evolVals[0]` = γ̇_eff (sum of absolute slip rates across
          * all slip systems); `evolVals[1]` = τ_s (`m_taus0` if the effective shear rate
          * is negligible).
          * @param[in] gdot Slip rates on all slip systems [1/time].
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
          * @brief Evaluate the Voce hardening rate and its derivative w.r.t. the
          * hardening state, for the scalar-hardness SNLS solve:
          *
          *   ḣ = h₀ · γ̇_eff · [(τ_s - h) / (τ_s - τ_si)]^m′
          *
          * (linear Voce when `nonlinear` is false, since m′ = `m_xmprime` == 1 in that
          * case). If the saturation stress τ_s (`sv_sat`) has fallen to or below τ_si
          * (`m_tausi`), the hardening rate is pinned to zero. See the file-level
          * documentation in ECMech_kinetics_VocePL.h for the full symbol-to-code
          * mapping.
          * @param[out] sdot Hardening rate ḣ.
          * @param[out] dsdot_ds Derivative of `sdot` w.r.t. `h`.
          * @param h Current CRSS, h.
          * @param[in] evolVals Values from getEvolVals: `evolVals[0]` = γ̇_eff
          * (effective shear rate), `evolVals[1]` = τ_s (saturation stress).
          * @param tkelv Temperature [Kelvin]; not currently used by this model.
          */
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
