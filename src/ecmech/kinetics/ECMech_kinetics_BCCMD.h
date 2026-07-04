/**
 * @file ECMech_kinetics_BCCMD.h
 * @brief BCC mobile-dislocation-density (MD) pencil-glide kinetics model.
 *
 * Couples a Taylor-hardening-type slip law (CRSS set by the total dislocation content)
 * with a Kocks-Mecking-style dislocation multiplication/recovery hardening law, tracking
 * one mobile dislocation density per slip system. Meant to be paired with the
 * stress-dependent pencil-glide BCC slip geometry (SlipGeomBCCPencil,
 * slipgeom/ECMech_slipgeom_bcc.h), which packs an extra per-slip-system angle
 * ("chi", the MRSSP orientation within the pencil-glide zone) alongside the resolved
 * shear stress; this model uses that angle to capture the twinning/anti-twinning
 * asymmetry characteristic of BCC pencil glide (both in the slip-rate law's reference
 * rate and in the hardening law's multiplication rate).
 *
 * **CRSS (Taylor hardening)**, shared by all slip systems (see KineticsBCCMD::getVals):
 *
 *   g = α · μ · b · √(Σᵢ ρᵢ)
 *
 * where g = `crss` (the shared CRSS), α = `m_alpha` (Taylor hardening coefficient), μ =
 * `m_mu` (shear modulus), b = `m_bmag` (Burgers vector magnitude), and ρᵢ =
 * `h_state[i]` (dislocation density on slip system i).
 *
 * **Slip-rate law** (see KineticsBCCMD::evalGdot): the Peierls stress and reference slip
 * rate both depend on the pencil-glide angle χ,
 *
 *   τ_p_eff = τ_p / cos(χ - α_p)
 *   τ_eff = max(|τ| - τ_p_eff, 0)
 *
 *   γ̇_w0(χ) = γ̇_th + (3/π)(χ + π/6)(γ̇_at - γ̇_th),   γ̇_at = 0.1 · γ̇_th
 *
 * (linearly interpolating between a "thermal" reference rate γ̇_th and an
 * "anti-thermal" reference rate γ̇_at over χ ∈ [-π/6, π/6], capturing the
 * twinning/anti-twinning asymmetry), then the Orowan relation and a phonon-drag
 * velocity limit give
 *
 *   γ̇_w = ρ · b · γ̇_w0(χ)
 *   γ̇_max = ρ · b · v_max · (1 - e^(-τ_eff/τ_drag))
 *
 * and finally a power law smoothly capped at γ̇_max:
 *
 *   γ̇_pl = γ̇_w · (τ_eff/g)^(1/m) · sign(τ)
 *   γ̇ = γ̇_pl / (1 + (|γ̇_pl|/γ̇_max)^a_c)^(1/a_c)
 *
 * where τ = `tau` (resolved shear stress), χ = `chi` (pencil-glide MRSSP angle), τ_p =
 * `m_tau_p` (Peierls stress), α_p = `m_alpha_p` (Peierls angular offset), γ̇_th =
 * `m_gam_w0` (thermal reference slip rate), ρ = `rho` (mobile dislocation density for
 * this slip system), b = `m_bmag` (Burgers vector magnitude), v_max = `m_vmax`
 * (maximum phonon-drag-limited dislocation velocity), τ_drag = `m_tau_drag` (drag
 * stress scale), g = `crss` (current CRSS), m = `m_xm` (rate-sensitivity exponent), a_c
 * = 10 (hard-coded smoothing sharpness), and γ̇ = `gdot` (slip rate).
 *
 * **Hardening law** (see KineticsBCCMD::getSdotN): per slip system i, a
 * Kocks-Mecking-style multiplication/recovery/relaxation ODE, solved implicitly in log
 * space (d(ln ρᵢ)/dt = (1/ρᵢ)·dρᵢ/dt) so the always-positive dislocation densities
 * cannot be driven negative by the nonlinear solve (see #updateHN):
 *
 *   dρᵢ/dt = k1(χᵢ)·|γ̇ᵢ|·√ρᵢ   -   k2·|γ̇ᵢ|·ρᵢ   -   f(|γ̇ᵢ|)·k_relax·kr(ρᵢ)·ρᵢ
 *            \_____________/       \___________/     \_______________________/
 *             multiplication          recovery                relaxation
 *
 *   k1(χ) = k1_0 · (1 + a_k / cos(χ - α_p))
 *   k2 = k2_0 · ln(γ̇_tot/γ̇_0) · ln(T/T_0)
 *   f(x) = 1 - 1 / (1 + exp(-A·(x/γ̇_tot - t)))
 *   kr(ρ) = 1 - exp(-(ρ - ρ_min)/ρ_min)
 *
 * where f is a smooth logistic gate (steep around x/γ̇_tot = t, with A=100, t=0.01)
 * that suppresses relaxation for slip systems contributing negligibly to the total
 * effective shear rate, and kr is a smooth ramp that suppresses relaxation as the
 * density approaches its floor ρ_min.
 *
 * Symbol-to-code mapping: ρᵢ = `h[i]` (mobile dislocation density on slip system i, in
 * linear space), γ̇ᵢ = `evolVals[i]` (= `|gdot_i|`, absolute slip rate on slip system i),
 * γ̇_tot = `gamma` (`evolVals[nEvolVals-1]`, total effective shear rate across all slip
 * systems), k1_0 = `m_k1` (multiplication-rate coefficient), a_k = `m_ak` (angular
 * multiplication-rate enhancement factor), k2_0 = `m_k2` (recovery-rate coefficient),
 * γ̇_0 = `m_gdot_0` (reference shear rate for the recovery term), T = `tkelv`
 * (temperature), T_0 = `m_tkelv0` (reference temperature for the recovery term),
 * k_relax = `m_krelax` (relaxation-rate coefficient), and ρ_min = `m_hdn_min`
 * (dislocation density floor).
 *
 * (The multiplication term also involves a forest-interaction matrix, currently
 * hard-coded to the identity here -- see KineticsBCCMD::getSdotN.)
 *
 * @see cases/ECMech_cases_bcc_defs.h for `Kin_BCC_MD`, the concrete alias pairing this
 * model with SlipGeomBCCPencil
 * @see ECMech_kinetics.h for the kinetics model interface contract this class
 * implements
 */

// -*-c++-*-

#ifndef ECMECH_KINETICS_BCCMD_H
#define ECMECH_KINETICS_BCCMD_H

#include <cassert>
#include <cmath>

/**
 * @def ECMECH_NN_INDX(p, q, nDim)
 * @brief Row-major flattening of a 2D `(p, q)` index into an `nDim x nDim` matrix.
 * Already defined identically by ECMech_util.h (included transitively via
 * ECMech_kinetics.h before this file); redefined here as well so this header remains
 * self-contained if ever included on its own.
 */
#define ECMECH_NN_INDX(p, q, nDim) (p) * (nDim) + (q)

namespace ecmech {
   /**
    * @brief BCC mobile-dislocation-density (MD) pencil-glide slip and hardening
    * kinetics.
    *
    * See the file-level documentation in ECMech_kinetics_BCCMD.h for the governing
    * equations. Expects to be paired with a slip geometry that supplies a
    * stress-dependent MRSSP angle ("chi") alongside the resolved shear stress in its
    * `tau` array -- in practice, SlipGeomBCCPencil.
    *
    * @tparam SlipGeom Slip geometry class; determines the number of slip systems
    * (`SlipGeom::nslip`) and is expected to provide the pencil-glide angle convention
    * this model relies on.
    *
    * @ingroup ECMech_kinetics
    */
   template<class SlipGeom>
   class KineticsBCCMD
   {
      public:
         /**
          * @brief Number of hardening state variables: one mobile dislocation density
          * per slip system (`SlipGeom::nslip`).
          */
         static constexpr int nH = SlipGeom::nslip;
         /** @brief Number of slip systems, mirrored from `SlipGeom::nslip` for readability elsewhere in this class. */
         static constexpr int m_num_slip = SlipGeom::nslip;
         /** @brief Number of parameters needed to construct/configure the model: 8 power-law/pencil-glide + 4 + 1 + 3 hardening-related (7 hardening parameters plus the initial dislocation density) = 16. */
         static constexpr int nParams = 8+4+1+3;
         /** @brief Number of kinetic values precomputed by getVals and reused by evalGdots: per-slip-system CRSS and dislocation density (2 * nslip), plus temperature (1). */
         static constexpr int nVals = 2 * SlipGeom::nslip + 1;
         /** @brief Number of intermediate values precomputed by getEvolVals and reused by getSdotN: per-slip-system absolute slip rate (nH), plus the total effective shear rate (1). */
         static constexpr int nEvolVals = nH + 1;

         /**
          * @brief Construct with a given number of slip systems (unused directly; the
          * slip system count is fixed by `SlipGeom::nslip`). Parameters must be set
          * separately via setParams.
          */
         __ecmech_hdev__
         KineticsBCCMD(int) {}
         /** @brief Destructor (default; no owned resources). */
         ~KineticsBCCMD() = default;

         /**
          * @brief Construct and immediately set parameters.
          * @param params Parameter array; see setParams for the expected order.
          */
         __ecmech_hdev__
         KineticsBCCMD(const double* const params, int)
         {
            setParams(params);
         }

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
          * Expected order: `mu` (shear modulus), `bmag` (Burgers vector magnitude),
          * `xm` (power-law exponent), `gam_w0` (reference slip rate), `tau_p` (Peierls
          * stress), `alpha_p` (Peierls angular offset), `vmax` (drag-limited shear
          * velocity), `tau_drag` (drag stress) -- then `alpha` (Taylor hardening
          * coefficient), `k1` (multiplication-rate coefficient), `k2` (recovery-rate
          * coefficient), `krelax` (relaxation-rate coefficient), `gdot_0` (reference
          * shear rate for the recovery term), `tkelv0` (reference temperature for the
          * recovery term), `ak` (angular multiplication-rate enhancement factor) -- then
          * `hdn_init` (initial dislocation density for every slip system). Also derives
          * the power-law exponent helpers (`m_xnn`, `m_xn`), the overflow/underflow
          * stress-ratio thresholds (`m_t_min`, `m_t_max`), and the dislocation-density
          * floor `m_hdn_min = 1e-4 * hdn_init`.
          * @param params Flat parameter array of length #nParams.
          */
         __ecmech_hdev__
         inline
         void setParams(const double* const params) {
            const double* parsIt = params;

            //////////////////////////////
            // power-law stuff
            // shear modulus in case the model uses it
            m_mu = *parsIt; ++parsIt;
            // Burgers vector magnitude
            m_bmag = *parsIt; ++parsIt;
            // This would be the power law exponent term
            m_xm = *parsIt; ++parsIt;
            // This would be the references slip rate term
            m_gam_w0 = *parsIt; ++parsIt;
            // Peierls stress
            m_tau_p = *parsIt; ++parsIt;
            // alpha Peierls
            m_alpha_p = *parsIt; ++parsIt;
            // Shear velocity
            m_vmax = *parsIt; ++parsIt;
            // Drag stress
            m_tau_drag = *parsIt; ++parsIt;

            // These are terms that are constant during the simulation and we don't
            // really need to calculate them every time we call slip kinetics portion
            // of the class
            m_xnn = one / m_xm;
            m_xn = m_xnn - one;
            //
            // CALL set_t_min_max(pl)
            // For numerics, we define a minimum and maximum (rss / crss) value
            // that translates to either a slip rate that is essentially zero
            // or slip rate that is going off to infinity but we really want to
            // cap it to some large number
            m_t_min = pow(ecmech::gam_ratio_min, m_xm);
            m_t_max = pow(ecmech::gam_ratio_ovf, m_xm);

            //////////////////////////////
            // Hardening parameters
            m_alpha = *parsIt; ++parsIt;
            m_k1 = *parsIt; ++parsIt;
            m_k2 = *parsIt; ++parsIt;
            m_krelax = *parsIt; ++parsIt;
            m_gdot_0 = *parsIt; ++parsIt;
            m_tkelv0 = *parsIt; ++parsIt;
            m_ak = *parsIt; ++parsIt;
            //////////////////////////////
            // nH
            // All the terms related to our hardening state
            // You'll often see the initial state provided to us saved off as well.
            // This is just so when you call getParam the initial state can be provided back
            // Also, if you'd like to based on the initial state you could provide a lower
            // bound for example related to the dislocation content that you don't want the
            // model to go under when you start try to update your hardening state.
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
            params.push_back(m_bmag);
            params.push_back(m_xm);
            params.push_back(m_gam_w0);
            params.push_back(m_tau_p);
            params.push_back(m_alpha_p);
            params.push_back(m_vmax);
            params.push_back(m_tau_drag);

            //////////////////////////////
            // hardening stuff
            params.push_back(m_alpha);
            params.push_back(m_k1);
            params.push_back(m_k2);
            params.push_back(m_krelax);
            params.push_back(m_gdot_0);
            params.push_back(m_tkelv0);
            params.push_back(m_ak);
            //////////////////////////////
            // nH

            params.push_back(m_hdn_init);

            //////////////////////////////
#ifdef ECMECH_DEBUG
            assert((params.size() - paramsStart) == nParams);
#endif
         }

         /**
          * @brief Describe the per-slip-system mobile dislocation density history
          * variables ("rho_0", "rho_1", ...) for history-array bookkeeping.
          * @param[out] names Appended with `"rho_" + iSlip` for each slip system.
          * @param[out] init Appended with #m_hdn_init for each slip system.
          * @param[out] plot Appended with `true` for each slip system.
          * @param[out] state Appended with `true` for each slip system.
          */
         __ecmech_host__
         void getHistInfo(std::vector<std::string> & names,
                          std::vector<double>       & init,
                          std::vector<bool>        & plot,
                          std::vector<bool>        & state) const {
            
            for (int iSlip = 0; iSlip < SlipGeom::nslip; iSlip++) {
               names.push_back("rho_" + std::to_string(iSlip));
               init.push_back(m_hdn_init);
               plot.push_back(true);
               state.push_back(true);
            }
         }

      private:

         //////////////////////////////
         // Power-law stuff

         /** @brief Shear modulus `mu` [stress units] and Burgers vector magnitude `bmag`. */
         double m_mu, m_bmag;
         /** @brief Power-law rate-sensitivity exponent `m`. */
         double m_xm;
         /** @brief Reference slip rate `gam_w0` [1/time]; also returned by getFixedRefRate. */
         double m_gam_w0;
         /** @brief Peierls stress `tau_p` [stress units] and its angular offset `alpha_p` used in the pencil-glide angular dependence `tau_p / cos(chi - alpha_p)`. */
         double m_tau_p, m_alpha_p;
         /** @brief Maximum (phonon-drag-limited) dislocation velocity `vmax` and the drag stress scale `tau_drag` controlling how quickly that limit is approached. */
         double m_vmax, m_tau_drag;

         // derived from parameters
         /** @brief Derived overflow/underflow thresholds on the stress ratio (see #gam_ratio_min, #gam_ratio_ovf) and power-law exponent helpers `xnn = 1/xm`, `xn = xnn - 1`. */
         double m_t_max, m_t_min, m_xn, m_xnn;

         //////////////////////////////
         // Hardening
         /** @brief Taylor hardening coefficient: CRSS = `alpha * mu * bmag * sqrt(total dislocation density)`. */
         double m_alpha;
         /** @brief Dislocation multiplication-rate coefficient `k1`, recovery-rate coefficient `k2`, and relaxation-rate coefficient `krelax` used in getSdotN's hardening ODE. */
         double m_k1, m_k2, m_krelax;
         /** @brief Reference shear rate `gdot_0` and reference temperature `tkelv0` used to scale the recovery-rate coefficient `k2`, and the angular multiplication-rate enhancement factor `ak` (see `k1_func` in getSdotN). */
         double m_gdot_0, m_tkelv0, m_ak;
         //////////////////////////////
         // nH
         /** @brief Initial mobile dislocation density for every slip system. */
         double m_hdn_init;
         /** @brief Floor on the dislocation density (`1e-4 * m_hdn_init`), used to keep the log-space hardening update and the relaxation-rate ramp (`k_relax_func` in getSdotN) well-defined. */
         double m_hdn_min;

      public:

         /**
          * @brief Reference slip rate used for scaling elsewhere in the solve (e.g. by
          * the evptn elastic-strain/rotation solver). For models where the reference
          * rate is slip-system-dependent, a maximum over slip systems would be computed
          * instead (see KineticsOrowanD::getVals for that pattern).
          * @return #m_gam_w0.
          */
         __ecmech_hdev__
         inline double getFixedRefRate(const double* const // vals, not used
                                       ) const
         {
            return m_gam_w0;
         }

         /**
          * @brief Precompute the kinetic values used by evalGdots: a Taylor-hardening
          * CRSS shared by all slip systems, each slip system's own dislocation density,
          * and temperature -- see the "CRSS" equation in the file-level documentation
          * in ECMech_kinetics_BCCMD.h.
          *
          * The single isotropic CRSS is duplicated across all `vals[0..nslip-1]`
          * entries so evalGdots can index into it uniformly. `vals[nslip..2*nslip-1]`
          * holds each slip system's own density (passed through from `h_state`, used by
          * evalGdot for the mobile-dislocation-density slip-rate scaling), and
          * `vals[2*nslip]` holds `tkelv`.
          * @param[out] vals Kinetic values, length #nVals; see above for the layout.
          * @param p Pressure; not currently used by this model.
          * @param tkelv Temperature [Kelvin], copied to `vals[2*nslip]`.
          * @param[in] h_state Current hardening state: dislocation density per slip
          * system.
          * @return The average CRSS across all slip systems (here just the single
          * shared CRSS value, since it's isotropic).
          */
         __ecmech_hdev__
         inline
         double
         getVals(double* const vals,
                 double, // p, not currently used
                 double tkelv,
                 const double* const h_state
                 ) const
         {
            double crss = ecmech::zero;
            for (int iSlip = 0; iSlip < m_num_slip; ++iSlip) {
               crss += h_state[iSlip];
               assert(h_state[iSlip] > zero);
            }
            crss = m_alpha * m_mu * m_bmag * sqrt(crss);
            
            double mVals = ecmech::zero;
            for (int iSlip = 0; iSlip < m_num_slip; ++iSlip) {
               vals[iSlip] = crss;
               vals[m_num_slip + iSlip] = h_state[iSlip];
               mVals += vals[iSlip];
               assert(vals[iSlip] > zero);
            }
            constexpr double inv_nslip = 1.0 / m_num_slip;
            mVals *= inv_nslip;
            
            vals[2*m_num_slip] = tkelv;

            return mVals;
         }

         /**
          * @brief Evaluate the slip rate and its derivative w.r.t. resolved shear stress
          * on every slip system.
          *
          * `tau` is expected to carry two pieces of per-slip-system data back to back:
          * the resolved shear stress in `tau[0..nslip-1]` and the pencil-glide MRSSP
          * angle ("chi") in `tau[nslip..2*nslip-1]` (as packed by SlipGeomBCCPencil).
          * @param[out] gdot Slip rate on each slip system [1/time].
          * @param[out] dgdot_dtau Derivative of slip rate w.r.t. resolved shear stress on
          * each slip system.
          * @param[in] tau Resolved shear stress (`tau[0..nslip-1]`) followed by the
          * pencil-glide angle `chi` (`tau[nslip..2*nslip-1]`).
          * @param[in] vals Kinetic values from getVals.
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
            double tkelv = 0.0;

            for (int iSlip = 0; iSlip < m_num_slip; ++iSlip) {
               bool l_act;
               double taua = tau[iSlip];
               double chia = tau[SlipGeom::nslip + iSlip];
               
               double crss = vals[iSlip];
               double rhoa = vals[SlipGeom::nslip + iSlip];
               // traditionally we have a separate function that will calculate everything
               // for only one slip system
               this->evalGdot(gdot[iSlip], l_act, dgdot_dtau[iSlip],
                              crss, rhoa, taua, chia, tkelv);
            }
         }

         /**
          * @brief Evaluate the pencil-glide power-law slip rate, smoothly capped by a
          * phonon-drag velocity limit, for a single slip system -- see the "Slip-rate
          * law" equations and symbol table in the file-level documentation in
          * ECMech_kinetics_BCCMD.h.
          *
          * The stress ratio τ_eff/g is compared against thresholds derived from
          * #gam_ratio_min/#gam_ratio_ovf (via `m_t_min`/`m_t_max`): below `m_t_min` the
          * system is inactive; above `m_t_max` the power-law rate is capped at
          * `gam_ratio_ovffx * gam_w` rather than evaluated directly (with derivatives
          * left at zero). If the angularly-interpolated γ̇_w0(χ) comes out negative,
          * γ̇_w = |γ̇_w0(χ)| is used directly, bypassing the density/Burgers-vector
          * scaling.
          * @param[out] gdot Slip rate γ̇ [1/time].
          * @param[out] l_act `true` if the slip system is active (stress ratio above the
          * underflow threshold).
          * @param[out] dgdot_dtau Derivative of `gdot` w.r.t. resolved shear stress.
          * @param crss Current CRSS g (shared across all slip systems; from getVals).
          * @param rho Current mobile dislocation density ρ for this slip system.
          * @param tau Resolved shear stress τ [stress units].
          * @param chi Pencil-glide MRSSP angle χ for this slip system.
          * @param tkelv Temperature [Kelvin]; not currently used by this model (always
          * called with `0.0` from evalGdots).
          */
         __ecmech_hdev__
         inline
         void
         evalGdot(
            double & gdot,
            bool  & l_act,
            double & dgdot_dtau, // wrt resolved shear stress
            double   crss,
            double   rho,
            double   tau,
            double   chi,
            double   /*tkelv*/
            ) const
         {
            // zero things so that can more easily just return in inactive
            //// gdot_w = zero; gdot_r = zero; ! not used by l_linear or l_pl
            gdot = zero;
            //
            dgdot_dtau = zero;            
            l_act = false;

            double tau_p = m_tau_p / cos(chi-m_alpha_p);

            double t_eff = fmax(fabs(tau) - tau_p, 0.0);

            double xnn = m_xnn;
            double xn = m_xn;
            double gam_w0 = m_gam_w0;

            double v0_T = m_gam_w0;
            double v0_AT = v0_T * 0.1;
            gam_w0 = v0_T + (M_PI/6.0 + chi) * 3.0/M_PI * (v0_AT - v0_T);

            double g_i = one / crss; // assume have checked gIn>0 elsewhere
            double t_frac = t_eff * g_i;
            t_frac = copysign(t_frac, tau); // has sign of tau
            double at = fabs(t_frac);
            
            double gam_w;
            if (gam_w0 < 0.0) gam_w = fabs(gam_w0);
            else gam_w = rho * m_bmag * gam_w0;
            
            double gmax = rho * m_bmag * m_vmax * (1.0-exp(-t_eff/m_tau_drag));

            if (at > m_t_min) {
               l_act = true;
               if (at > m_t_max) {
                  // ierr = IERR_OVF_p
                  // set gdot big, evpp may need this for recovery
                  gdot = ecmech::gam_ratio_ovffx * gam_w;
                  gdot = copysign(gdot, tau);
                  // do not set any of deriviatives (they are, in truth, zero)
               }
               else {
                  double abslog = log(at);
                  double blog = xn * abslog;
                  double temp = gam_w * exp(blog);

                  gdot = temp * t_frac;

                  dgdot_dtau = temp * xnn * g_i; // note: always positive, = xnn * gdot/t
                  // Smooth capping to gmax with Lorentz-like factor
                   double ac = 10.0;
                   double gfrac_a = pow(fabs(gdot)/gmax, ac);
                   double dfact = pow(1.0 + gfrac_a, -(ac+1.0)/ac);
                   
                   // Correct derivative accounting for gmax(tau)
                   double dgmax = rho * m_bmag * m_vmax / m_tau_drag * exp(-t_eff/m_tau_drag);
                   double A = fabs(gdot) * dgmax * gfrac_a;
                   dgdot_dtau = dfact * (A + gmax * dgdot_dtau) / gmax;
                   gdot = gdot / pow(1.0 + gfrac_a, 1.0/ac);
               }
            }
         } // evalGdot

         /**
          * @brief Advance the per-slip-system mobile dislocation densities one time
          * step, solving implicitly in log space so the densities cannot go negative.
          *
          * Called externally by the elastic-strain/lattice-rotation update at the
          * beginning of the time step (not iteratively), so all inputs are
          * beginning-of-step values. Precomputes the evolution inputs via getEvolVals
          * and passes them through to the shared vector-hardness solve (#updateHN, with
          * `relaxed_solver = true` since this log-space, angularly-coupled ODE can be
          * numerically stiff).
          * @param[out] hs_u Updated (end-of-step) dislocation densities [SlipGeom::nslip].
          * @param[in] hs_o Beginning-of-step dislocation densities [SlipGeom::nslip];
          * floored at #m_hdn_min before taking the log.
          * @param dt Time step size.
          * @param[in] gdot Beginning-of-step slip rates on all slip systems.
          * @param[in] hvals Per-slip-system pencil-glide angles ("chi"), forwarded to
          * getSdotN.
          * @param tkelv Temperature [Kelvin].
          * @param outputLevel Optional; passed through to the SNLS solver for logging.
          * @return Function-evaluation count from updateHN, or a negative value if the
          * solve failed to converge (even with the relaxed retry).
          * @see updateHN in ECMech_kinetics.h
          */
         __ecmech_hdev__
         inline
         int
         updateH(double* const hs_u,
                 const double* const hs_o,
                 double dt,
                 const double* const gdot,
                 const double* const hvals,
                 double tkelv,
                 int outputLevel = 0) const
         {
            double log_hs_u[SlipGeom::nslip];
            double log_hs_o[SlipGeom::nslip];
            double evolVals[nEvolVals] = {};
            
            for(int islip = 0; islip < SlipGeom::nslip; islip++) {
               log_hs_o[islip] = log(fmax(hs_o[islip], m_hdn_min));
            }
            getEvolVals(evolVals, gdot);
            // If the equation is incredibly  stiff it's possible this won't solve
            int nFEvals = updateHN<KineticsBCCMD, true>(this,
                                                  log_hs_u, log_hs_o, dt, evolVals, hvals, tkelv,
                                                  outputLevel);

            for(int islip = 0; islip < SlipGeom::nslip; islip++) {
               hs_u[islip] = exp(log_hs_u[islip]);
            }
            return nFEvals;
         }

         /**
          * @brief Precompute the per-slip-system absolute slip rates and the total
          * effective shear rate used by getSdotN.
          * @param[out] evolVals `evolVals[0..nslip-1]` = `|gdot[i]|` for each slip
          * system; `evolVals[nEvolVals-1]` **accumulates** (`+=`, not `=`) the sum of all
          * `|gdot[i]|` (the total effective shear rate, "gamma"), so callers must
          * zero-initialize `evolVals` before calling this.
          * @param[in] gdot Slip rates on all slip systems [1/time].
          * @note updateH (below) zero-initializes its `evolVals` before calling this
          * directly, then passes that already-summed array through to #updateHN's
          * `gdot` parameter -- which calls this function a second time, into its own
          * *non*-zero-initialized local `evolVals` (see `ECMech_kinetics.h`). That
          * second call's `evolVals[nEvolVals-1]` therefore accumulates onto
          * uninitialized memory rather than starting from zero, which looks like an
          * unintentional bug in the interaction between this model's accumulation
          * pattern and updateHN's generic buffer.
          */
         __ecmech_hdev__
         inline
         void
         getEvolVals(double* const evolVals,
                     const double* const gdot
                     ) const
         {
            for (int i = 0; i < m_num_slip; i++) {
                evolVals[i] = abs(gdot[i]);
                evolVals[nEvolVals - 1] += abs(gdot[i]);
            }
         }

         /**
          * @brief Evaluate the Kocks-Mecking-style dislocation multiplication/recovery
          * hardening rate and its Jacobian, in log space, for the vector-hardness SNLS
          * solve (#updateHN) -- see the "Hardening law" equations and symbol table in
          * the file-level documentation in ECMech_kinetics_BCCMD.h.
          *
          * Code-to-symbol correspondence for the local helper lambdas: `k1_func(xi)` is
          * k1(χ), `k2_func()` is k2, `f_func(abs_gamma_dot)` is f(x), and
          * `k_relax_func(h)` is kr(ρ). The multiplication term uses a forest-interaction
          * matrix `m_a_mat`, but it is currently hard-coded to the identity matrix here
          * (a simplification noted in the surrounding code as a candidate for a future,
          * fully anisotropic interaction matrix, as used in KineticsOrowanD).
          *
          * @param[out] sdot Log-space hardening rate d(ln ρᵢ)/dt for each slip system.
          * @param[out] dsdot_ds Log-space Jacobian of `sdot` w.r.t. `h`
          * [nDimSys x nDimSys], row-major via #ECMECH_NN_INDX; left untouched if
          * `nullptr`.
          * @param[in] h_i Current hardening state, in log space (`h = exp(h_i)` = ρ,
          * floored at #idp_eps).
          * @param[in] evolVals Values from getEvolVals.
          * @param[in] hvals Per-slip-system pencil-glide angles ("chi", χᵢ).
          * @param tkelv Temperature [Kelvin], T.
          */
         __ecmech_hdev__
         inline
         void
         getSdotN(double *sdot,
                  double *dsdot_ds,
                  const double* const h_i,
                  const double* const evolVals,
                  const double* const hvals,
                  double tkelv
                ) const
         {
            constexpr bool LOGFORM = true;
            constexpr size_t nslip = SlipGeom::nslip;
            constexpr size_t JDIM = 2;
            constexpr size_t nDimSys = SlipGeom::nslip;
            constexpr size_t h_content = (LOGFORM) ? nslip : 1;

            double hexp[h_content];
            if (LOGFORM) {
               for (size_t iDD = 0; iDD < nslip; iDD++) {
                  // Prevent divide by 0 errors...
                  hexp[iDD] = fmax(exp(h_i[iDD]), ecmech::idp_eps);
               }
            }
            const double* const h = (LOGFORM) ? &hexp[0] : h_i;

            const double gamma = evolVals[nEvolVals - 1];
            auto k1_func = [=](const double xi) -> double {
               /*
               Current version does things a bit different from the paper.
               double chia = hvals[islip];
               // If we are in the AT zone, then we need to increase k1
               // to account for the fact that dislocations do take
               // a longer path and thus are likely to multiply more
               double amin = 0.95;
         
               //double a = fmin(1.0 + (amin - 1.0) * xi * 6.0 / M_PI, 1.0);
               //a = 1.0 / a;
         
               double a = 1.0/(1.0/cos(M_PI/6.0 - m_alpha_p)-1.0) * (1.0/amin - 1.0);
               a = 1.0 + a * (1.0 / cos(xi - m_alpha_p) - 1.0);
         
               return m_k1 * a;
               */
               // Pure paper implementation of things
               return m_k1 * (1.0 +(m_ak / (cos(xi - m_alpha_p))));
            };

            auto k2_func = [=] () -> double {
               const double gamma_ratio = (gamma > ecmech::gam_ratio_min) ? (gamma / m_gdot_0) : ecmech::gam_ratio_min;
               return m_k2 * log(gamma_ratio) * log(tkelv / m_tkelv0);
            };

            auto f_func = [=] (const double abs_gamma_dot) -> double {
               constexpr double A = 100.0;
               constexpr double t = 0.01;
               const double gamma_ratio = (gamma > ecmech::gam_ratio_min) ? (abs_gamma_dot / gamma) : gam_ratio_max;
               const double exp_inner = -A * (gamma_ratio - t);
               return 1.0 - (1.0 / ( 1.0 + exp(exp_inner)));
            };

            auto k_relax_func = [=] (const double h) -> double {
               const double relax_term = 1.0 - exp(-(h - m_hdn_min) / m_hdn_min);
               return relax_term;
            };

            // From the paper if sqrt(a_ij * rho_j) does not correspond to
            // sqrt(I_ij * rho_j) where I is the identity matrix
            // then we'd something like the below
            double m_a_mat[nslip * nslip] = {};
            {
               RAJA::View<double, RAJA::Layout<JDIM> > amat(&m_a_mat[0], nslip, nslip);
               for (size_t islip = 0; islip < nslip; islip++) {
                  amat(islip, islip) = 1.0;
               }
            }
            double amat_rho[SlipGeom::nslip] = {};
            vecsVMa<SlipGeom::nslip>(&amat_rho[0], &m_a_mat[0], &h[0]);

            for (int islip = 0; islip < SlipGeom::nslip; islip++) {
               const double k1   = k1_func(hvals[islip]) * evolVals[islip];
               const double k2   = k2_func() * evolVals[islip];
               const double fval = f_func(evolVals[islip]) * m_krelax * k_relax_func(h[islip]);
               amat_rho[islip] = sqrt(amat_rho[islip]);
               sdot[islip] = (k1 * amat_rho[islip] - k2 * h[islip]) - fval * h[islip];
            }

            if (LOGFORM) {
               for (int iDD = 0; iDD < SlipGeom::nslip; iDD++) {
                  sdot[iDD] *= (1.0 / h[iDD]);
               }
            }
            if (dsdot_ds)
            {
               // zero out dsdot_ds matrix
               for (size_t i = 0; i < nDimSys * nDimSys; i++) {
                  dsdot_ds[i] = ecmech::zero;
               }
               RAJA::View<const double, RAJA::Layout<JDIM> > amat(&m_a_mat[0], nslip, nslip);
               RAJA::View<double, RAJA::Layout<JDIM> > dsdot_ds_view(dsdot_ds, nDimSys, nDimSys);
               for (size_t islip = 0; islip < nslip; islip++) {
                  const double k1   = k1_func(hvals[islip]) * evolVals[islip];
                  const double k2   = k2_func() * evolVals[islip];
                  const double fval = f_func(evolVals[islip]) * m_krelax;
                  const double ratio = h[islip] / m_hdn_min;
                  // Don't want this getting too big or else we get an infinity error later on...
                  const double ratio_max = fmin(ratio, 80.0);
                  // Copied this from Wolfram alpha for the f_func() * k_relax * k_relax_func() * h
                  const double fval_der = fval - (fval * exp(1 - ratio_max) * (m_hdn_min - h[islip]))/m_hdn_min;
                  dsdot_ds_view(islip, islip) += -(k2 + fval_der);
                  for (size_t jslip = 0; jslip < nslip; jslip++) {
                     dsdot_ds_view(islip, jslip) += k1 * amat(islip, jslip) * 0.5 / amat_rho[jslip];
                  }
               }
               // Generic solution to transform over into log space
               if (LOGFORM) {
                  for (size_t iDD = 0; iDD < nslip; iDD++) {
                     dsdot_ds_view(iDD, iDD) -= sdot[iDD];
                  }
                  for (size_t iDD = 0; iDD < nslip; iDD++) {
                     for (size_t jDD = 0; jDD < nslip; jDD++) {
                        dsdot_ds_view(iDD, jDD) *= h[jDD] / h[iDD];
                     }
                  }
               }
            }
         }
   }; // class KineticsBCCMD
} // namespace ecmech

#endif // ECMECH_KINETICS_BCCMD_H
