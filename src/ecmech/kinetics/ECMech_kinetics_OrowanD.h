/**
 * @file ECMech_kinetics_OrowanD.h
 * @brief Mobile/total dislocation-density (Orowan) hardening law, paired with balanced
 * thermally-activated MTS-like slip kinetics and phonon drag.
 *
 * This model tracks two dislocation-density state variables per slip system: a
 * **mobile** density (dislocations currently able to glide) and a **total** density
 * (mobile + immobile/forest). Both are **absolute** densities with physical units
 * (contrast with KineticsKMBalD's single, dimensionless relative density). The
 * slip-rate law shares the same balanced (bidirectional), thermally-activated-plus-drag
 * form as KineticsKMBalD, but is evaluated per slip system from that slip system's own
 * mobile density via the Orowan relation γ̇ = ρ_m·b·v; the hardening law is its own,
 * extending beyond a single Kocks-Mecking density with explicit dislocation
 * multiplication, trapping, and annihilation terms (see "Hardening law" below).
 *
 * **Template parameters**:
 * - `withGAthermal`, `pOne`, `qOne`: select the athermal-floor/MTS-normalizing-stress
 *   role assignment (between the per-slip-system CRSS ĝ and the reference stress τ_a =
 *   `m_tau_a`) and whether the MTS exponents p / q are fixed at 1 (see "Slip-rate law"
 *   below).
 * - `isotropic`: when `true`, a single scalar (`m_inter_mat[0]`) is used as the
 *   forest-interaction strength between every pair of slip systems, rather than a full
 *   per-pair interaction matrix (see "Kinetic values" below).
 * - `perSS`: when `true`, the per-group parameters (`m_c_1`, `m_c_2`, `m_berg_mag`) are
 *   given one-per-slip-system (`nVPer` must equal the slip system count) rather than
 *   shared across all slip systems (`nVPer` must be 1).
 * - `nVPer`: number of parameter groups; 1 if `!perSS`, else the slip system count.
 * - `SlipGeom`: the slip geometry class; also used once, at `setParams` time, to
 *   compute the forest-interaction matrix used by the hardening law (see "Hardening
 *   law" below).
 * - `LOGFORM`: when `true`, `updateH` solves for both densities in log space (so they
 *   cannot be driven negative by the nonlinear solve); when `false`, negative results
 *   are instead caught and retried by substepping (see `updateH`).
 *
 * **Kinetic values** (see KineticsOrowanD::getVals): per slip system i, a
 * Taylor-hardening-type CRSS driven by the forest (total) dislocation density, plus a
 * combined thermal-activation/phonon-drag reference rate:
 *
 *   forestᵢ = Σⱼ A_interᵢⱼ · qTⱼ      (isotropic: forestᵢ = a_inter · Σⱼ qTⱼ)
 *   ĝᵢ = c_2ᵢ · √forestᵢ
 *   rateᵢ = 1 / (1/γ̇_wᵢ + 1/γ̇_rᵢ),      γ̇_wᵢ = (L̄/b)·f_D / √qMᵢ,      γ̇_rᵢ = γ̇_r0 · qMᵢ
 *   c_tᵢ = C1ᵢ / T
 *
 * where:
 * - forestᵢ = internal only (not stored in `vals`) -- forest dislocation density seen
 *   by slip system i
 * - A_interᵢⱼ / a_inter = `m_inter_mat` -- user-supplied forest-interaction strength
 *   between slip systems i and j (a single scalar `a_inter` if `isotropic`, else one
 *   entry per pair); this is a *different* matrix from the one the hardening law uses
 *   (`m_a_mat`, see "Hardening law" below) -- the two coincide only if `ORO_USE_INTERMAT`
 *   is defined
 * - qTⱼ = `h_state[nslip+j]` -- total dislocation density on slip system j
 * - ĝᵢ = `vals[1+i]` -- per-slip-system CRSS; plays the role of either the athermal
 *   floor or the MTS-normalizing stress in the slip-rate law below, depending on
 *   `withGAthermal`
 * - c_2ᵢ = `m_c_2[i]` -- Taylor hardening coefficient
 * - γ̇_wᵢ -- thermally-activated reference rate for slip system i (reused in the
 *   slip-rate law below)
 * - L̄/b = `m_lbar_b` -- mean forest spacing over Burgers vector
 * - f_D = `m_fD` -- thermal attempt-frequency factor
 * - qMᵢ = `vals[1+nslip+i]` (= `h_state[i]`) -- mobile dislocation density on slip
 *   system i
 * - γ̇_rᵢ -- phonon-drag reference rate for slip system i (reused in the slip-rate law
 *   below)
 * - γ̇_r0 = `m_gam_ro` -- reference-rate coefficient for the drag-limited branch
 * - rateᵢ -- per-slip-system representative reference rate; its maximum over all slip
 *   systems becomes `vals[0]`, returned by getFixedRefRate as the model's overall
 *   reference rate
 * - c_tᵢ = `vals[1+2*nslip+i]` -- thermal energy scale
 * - C1ᵢ = `m_c_1[i]` -- thermal-activation energy scale numerator
 * - T = `tkelv` -- temperature
 *
 * The function returns the average of ĝᵢ across all slip systems.
 *
 * **MTS thermal-activation energy function** (see KineticsOrowanD::get_mts_dG), a
 * Kocks-Argon-Ashby-style activation-energy profile (identical in form to
 * KineticsKMBalD::get_mts_dG):
 *
 *   E(t) = -c_e · [1 - sign(t)·|t|^p]^q
 *
 * so that a thermally-activated rate is `(reference rate) · exp(E(t))` -- `E` pegs to 0
 * once `t ≥ 1` (barrier fully overcome by stress) and grows more negative (stronger
 * suppression) as `t` decreases. `t` itself is where the resolved shear stress and CRSS
 * actually enter the thermal-activation rate:
 *
 *   t = (σ - g_ath) / g_MTS
 *
 * for a signed driving-stress term σ (the forward evaluation uses σ = |τ|; the
 * reverse/balancing evaluation uses σ = -|τ|; see "Slip-rate law" below) and the
 * athermal-floor/MTS-normalizing-stress pair g_ath/g_MTS (assigned from ĝ and τ_a
 * depending on `withGAthermal`, also in "Slip-rate law" below). Here:
 * - E = `exp_arg` -- the log of the thermal-activation rate factor
 * - c_e = `c_e` -- thermal energy prefactor (= c_t·μ)
 * - p = `m_p`, q = `m_q` -- MTS activation-energy exponents
 * - t = `t_frac` -- the dimensionless MTS argument passed in
 *
 * **Slip-rate law** (see KineticsOrowanD::evalGdot): depending on `withGAthermal`, the
 * athermal stress floor g_ath and MTS-normalizing stress g_MTS are assigned from the
 * per-slip-system CRSS ĝ and the reference stress τ_a in one of two ways:
 *
 *   withGAthermal:   g_ath = ĝ,     g_MTS = τ_a
 *   !withGAthermal:  g_ath = τ_a,   g_MTS = ĝ
 *
 * The thermally-activated (balanced/bidirectional) rate and the drag-limited rate are
 *
 *   γ̇_th = γ̇_w · [exp(E(t_fwd)) - exp(E(t_rev))]
 *   t_fwd = (|τ| - g_ath) / g_MTS,   t_rev = -(|τ| + g_ath) / g_MTS
 *
 *   γ̇_drag = γ̇_r · (1 - exp(-(|τ| - g_ath)/w_rD))
 *
 * (the reverse-hop term is dropped once negligible), plus a high-stress power-law tail
 * once the clamped forward argument `max(0, t_fwd)` exceeds an underflow threshold
 * t_min (mirroring #gam_ratio_min/#gam_ratio_ovf, via an effective rate-sensitivity
 * exponent derived below):
 *
 *   γ̇_pl = 10 · γ̇_w · max(0, t_fwd)^xnn
 *
 * and finally the thermal and drag branches are combined by harmonic mean (as
 * resistances combine in series), since either mechanism being much slower than the
 * other dominates the net rate:
 *
 *   γ̇ = 1 / (1/(γ̇_th + γ̇_pl) + 1/γ̇_drag) · sign(τ)
 *
 * where:
 * - g_ath, g_MTS -- athermal stress floor and MTS-normalizing stress for this slip
 *   system (assigned above)
 * - τ_a = `m_tau_a` -- reference stress
 * - γ̇_w, γ̇_r -- the per-slip-system reference rates already computed by getVals (see
 *   "Kinetic values" above)
 * - τ = `tau` -- resolved shear stress
 * - w_rD = `m_wrD` -- drag stress scale
 * - xnn -- power-law exponent helper, `1/xm` (see below)
 * - γ̇ = `gdot` -- slip rate
 *
 * Once the clamped forward argument exceeds an overflow threshold t_max, the
 * thermally-activated part is treated as having overflowed and the rate is purely
 * drag-limited (γ̇ = γ̇_drag).
 *
 * An effective power-law rate-sensitivity exponent is derived once, at `setParams`
 * time, from the reference MTS parameters, so the power-law tail is asymptotically
 * consistent with the MTS thermal part at high stress ratios (`xnn = 1/xm`, `xn = xnn -
 * 1`, matching the `t_min`/`t_max` convention used by KineticsKMBalD):
 *
 *   xmᵢ = 1 / (2 · (C1ᵢ/T_ref) · μ_ref · p · q)
 *
 * where T_ref = `m_tkelv_ref` and μ_ref = `m_mu_ref` (reference temperature and shear
 * modulus).
 *
 * **Hardening law** (see KineticsOrowanD::getEvolVals / KineticsOrowanD::getSdotN): per
 * slip system i, the mobile and total dislocation densities evolve by multiplication,
 * trapping (mobile becoming immobile/forest), and annihilation, solved implicitly
 * (optionally in log space if `LOGFORM`, so densities cannot go negative by
 * construction; see `updateH` for the non-log-space fallback otherwise):
 *
 *   dqMᵢ/dt = c_mult·√fᵢ·qMᵢ·vᵢ   -   c_trap·√fᵢ·qMᵢ·vᵢ   -   c_ann·d_ann·qMᵢ²·vᵢ
 *   dqTᵢ/dt = c_mult·√fᵢ·qMᵢ·vᵢ                            -   c_ann·d_ann·qMᵢ²·vᵢ
 *   fᵢ = Σⱼ A_matᵢⱼ · qTⱼ
 *   vᵢ = |γ̇ᵢ| / (qMᵢ·bᵢ)
 *
 * where multiplication and trapping share the same rate expression (they only differ
 * by coefficient, since both represent existing mobile dislocations sweeping through
 * the forest at velocity v), and total density is unaffected by trapping (which only
 * moves density from the mobile to the immobile/forest bucket, not out of the total
 * count). Here:
 * - qMᵢ, qTᵢ -- mobile and total dislocation density on slip system i (the two halves
 *   of the hardening state vector `h`)
 * - c_mult = `m_c_mult` -- multiplication-rate coefficient
 * - c_trap = `m_c_trap` -- mobile-to-forest trapping-rate coefficient
 * - c_ann = `m_c_ann`, d_ann = `m_d_ann` -- annihilation-rate coefficient and capture
 *   distance
 * - fᵢ -- forest density seen by slip system i
 * - A_matᵢⱼ = `m_a_mat` -- forest-interaction matrix used specifically by the
 *   hardening law. Unless `ORO_USE_INTERMAT` is defined (in which case it's simply
 *   copied from `m_inter_mat`, see "Kinetic values" above), it is instead computed
 *   geometrically at `setParams` time from each pair of slip systems' normals (m) and
 *   directions (s):
 *
 *     A_matᵢⱼ = (1/2)·(|mᵢ·sⱼ| + |mᵢ·(mⱼ×sⱼ)|)
 *
 *   combining the Schmid-type overlap between systems with an out-of-plane/cross term.
 *   Building this requires the slip system m/s vectors, which are otherwise not owned
 *   by this class -- so a temporary `SlipGeom` instance is constructed in `setParams`
 *   from a second copy of the slip-geometry parameters appended to the end of this
 *   model's own parameter array (hence the `+ SlipGeom::nParams` term in #nParams).
 * - vᵢ = `evolVals[i]` (= `nu[i]` in `updateH`) -- dislocation glide velocity on slip
 *   system i (Orowan relation), precomputed by `updateH` from the beginning-of-step
 *   slip rate and passed through unchanged by getEvolVals
 * - γ̇ᵢ -- slip rate on slip system i (the `gdot` passed in to `updateH`)
 * - bᵢ = `m_berg_mag[i]` -- Burgers vector magnitude
 *
 * @see ECMech_kinetics.h for the kinetics model interface contract this class
 * implements, and KineticsVocePL/KineticsKMBalD/KineticsBCCMD for the other available
 * kinetics models
 */

// -*-c++-*-

#ifndef ECMECH_KINETICS_OROWAND_H
#define ECMECH_KINETICS_OROWAND_H

#include <cassert>
#include <cmath>

#include <string>
#include <vector>

#include "RAJA/RAJA.hpp"

namespace ecmech {
   /**
    * @brief Mobile/total dislocation-density (Orowan) hardening law with balanced,
    * thermally-activated MTS-like slip kinetics and phonon drag.
    *
    * See the file-level documentation in ECMech_kinetics_OrowanD.h for the governing
    * equations and template-parameter meanings.
    *
    * @tparam withGAthermal Selects which of the per-slip-system CRSS or the reference
    * stress `m_tau_a` is the athermal floor vs. the MTS-normalizing stress.
    * @tparam pOne When `true`, the MTS exponent p is fixed at 1.
    * @tparam qOne When `true`, the MTS exponent q is fixed at 1.
    * @tparam isotropic When `true`, the forest-interaction matrix collapses to a single
    * scalar shared by every slip-system pair.
    * @tparam perSS When `true`, the per-group parameters vary per slip system.
    * @tparam nVPer Number of parameter groups (1 if `!perSS`, else the slip system count).
    * @tparam SlipGeom Slip geometry class.
    * @tparam LOGFORM When `true`, solve the hardening update in log space.
    *
    * @ingroup ECMech_kinetics
    */
   template<bool withGAthermal,
            bool pOne, // l_p_1
            bool qOne, // l_q_1
            bool isotropic, // H^{\alpha\beta} = 1 so isotropic interaction matrix
            bool perSS, // If varying params per SS usually used for non-cubic materials
            int nVPer, // If perSS then nVPer should equal nslip
            class SlipGeom,
            bool LOGFORM = false> // LOGFORM dictates whether or not we use a logrithmic form for our hardness update
   class KineticsOrowanD
   {
      public:
         /** @brief Number of hardening state variables: one mobile plus one total dislocation density per slip system. */
         static constexpr int nH = 2 * SlipGeom::nslip; // Number of mobile and total dislocation density
         /** @brief Number of independent forest-interaction-matrix entries: 1 if `isotropic`, else one per slip-system pair. */
         static constexpr int nIH = isotropic ? 1 : (SlipGeom::nslip * SlipGeom::nslip); // Number of params in interaction matrix
         /** @brief Number of parameters: 12 MTS/power-law + 4 per group + nH initial densities + nIH interaction-matrix entries + a second copy of the slip geometry's own parameters (see "Hardening law" in the file-level documentation). */
         static constexpr int nParams = 12 + 4 * nVPer + nH + nIH + SlipGeom::nParams;
         /** @brief Number of kinetic values precomputed by getVals and reused by evalGdots: reference slip rate (1), per-slip-system CRSS and mobile density (2 * nslip), plus C1/T per group (nVPer). */
         static constexpr int nVals = 1 + nVPer + 2 * SlipGeom::nslip; //Our ref_slip_rate, CRSS, C1/T, and b*q_m params
         /** @brief Number of intermediate values precomputed by getEvolVals and reused by getSdotN: one dislocation velocity per slip system. */
         static constexpr int nEvolVals = SlipGeom::nslip; // We really don't need to evolve anything here
         /**
          * @brief Construct with a given number of slip systems; parameters must be set
          * separately via setParams.
          * @param _nslip Number of slip systems; must equal `SlipGeom::nslip`, and must
          * equal `nVPer` if `perSS` (otherwise `nVPer` must be 1).
          */
         __ecmech_hdev__
         KineticsOrowanD(int _nslip) : nslip(_nslip) {
            assert(nslip == SlipGeom::nslip);
            if (perSS) {
               assert(nslip == nVPer);
            }
            else {
               assert(nVPer == 1);
            }
         }
         /** @brief Destructor (default; no owned resources). */
         ~KineticsOrowanD() = default;

         /**
          * @brief Construct with a given number of slip systems and immediately set
          * parameters.
          * @param params Parameter array; see setParams for the expected order.
          * @param _nslip Number of slip systems; must equal `SlipGeom::nslip`, and must
          * equal `nVPer` if `perSS` (otherwise `nVPer` must be 1).
          */
         __ecmech_hdev__
         KineticsOrowanD(const double* const params, int _nslip) :
         nslip(_nslip)
         {
            assert(nslip == SlipGeom::nslip);
            if (perSS) {
               assert(nslip == nVPer);
            }
            else {
               assert(nVPer == 1);
            }
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
          * Expected order: `mu_ref` (μ_ref), `tkelv_ref` (T_ref), `berg_mag[nVPer]`
          * (Burgers vector magnitude bᵢ), `lbar_b` (mean forest spacing over Burgers
          * vector, L̄/b), `gam_ro` (γ̇_r0), `wrD` (w_rD), `fD` (thermal attempt-frequency
          * factor), `c_1[nVPer]` (C1ᵢ), `tau_a` (τ_a), `p`, `q`, `c_2[nVPer]` (Taylor
          * hardening coefficient per group), `inter_mat[nIH]` (forest-interaction values
          * used directly by getVals; also used in place of the geometrically-computed
          * matrix by getSdotN if `ORO_USE_INTERMAT` is defined -- see the "Kinetic
          * values" and "Hardening law" sections in the file-level documentation) -- then
          * `c_ann`, `d_ann`, `c_trap`,
          * `c_mult` (dislocation evolution coefficients) -- then `qM[nslip]` (initial
          * mobile densities), `qT[nslip]` (initial total densities) -- then a second
          * copy of the slip geometry's own parameters, used only to reconstruct the m/s
          * vectors needed for the forest-interaction matrix (see the file-level
          * documentation in ECMech_kinetics_OrowanD.h). Also derives, per group, the
          * effective power-law rate-sensitivity exponent and the overflow/underflow
          * stress-ratio thresholds `m_xnn`/`m_xn`/`m_t_min`/`m_t_max` (same construction
          * as KineticsKMBalD), and the dislocation-density floor `m_hdn_min` (`1e-4`
          * times the smallest initial mobile density).
          * @param params Flat parameter array of length #nParams.
          */
         __ecmech_hdev__
         inline
         void setParams(const double* const params) {
            const double* parsIt = params;

            //////////////////////////////
            // power-law stuff

            m_mu_ref = *parsIt; ++parsIt;
            m_tkelv_ref = *parsIt; ++parsIt;
            for (int iVal = 0; iVal < nVPer; ++iVal) {
               m_berg_mag[iVal] = *parsIt; ++parsIt;
            }

            m_lbar_b = *parsIt; ++parsIt;

            m_gam_ro = *parsIt; ++parsIt;
            m_wrD = *parsIt; ++parsIt;

            // thermal activation params
            m_fD = *parsIt; ++parsIt;
            for (int iVal = 0; iVal < nVPer; ++iVal) {
               m_c_1[iVal] = *parsIt; ++parsIt;
            }

            m_tau_a = *parsIt; ++parsIt;
            m_p = *parsIt; ++parsIt;
            m_q = *parsIt; ++parsIt;
            for (int iVal = 0; iVal < nVPer; ++iVal) {
               m_c_2[iVal] = *parsIt; ++parsIt;
            }

            for (int iVal = 0; iVal < nIH; ++iVal) {
               m_inter_mat[iVal] = *parsIt; ++parsIt;
            }

            if (withGAthermal) {
               assert(m_tau_a > zero);
            }
            if (pOne) {
               assert(m_p == one);
            }
            if (qOne) {
               assert(m_q == one);
            }

            // Figure out what the equivalent from the KMBalD is for the down below
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
               // These factors are the same from the balanced-MTS kinetic mobility law
               // and if the ratio tau/crss is below m_t_min we aren't moving at all
               // and if the ratio is above m_t_max we're strictly in phonon drag mobility
               m_t_min[iVal] = pow(ecmech::gam_ratio_min, xm);
               m_t_max[iVal] = pow(ecmech::gam_ratio_ovf, xm);

            }

            //////////////////////////////
            // Dislocation evolution stuff

            m_c_ann = *parsIt; ++parsIt;
            m_d_ann = *parsIt; ++parsIt;
            m_c_trap = *parsIt; ++parsIt;
            m_c_mult = *parsIt; ++parsIt;

            //////////////////////////////
            // Dislocation Densities
            for (int iVal = 0; iVal < SlipGeom::nslip; iVal++) {
               m_qM[iVal] = *parsIt; ++parsIt;
            }

            for (int iVal = 0; iVal < SlipGeom::nslip; iVal++) {
               m_qT[iVal] = *parsIt; ++parsIt;
            }

            m_hdn_min = m_qM[0];
            // The mobile dd should be the smallest so find the smallest
            // here and base our m_hdn_min on that.
            for (int iVal = 0; iVal < SlipGeom::nslip; iVal++) {
               if (m_hdn_min > m_qM[iVal]) {
                  m_hdn_min = m_qM[iVal];
               }
            }
            // Might want to make this smaller if provided large initial DD value?
            m_hdn_min *= 1.0e-4;

            //////////////////////////////
            // Initialize slip system matrix
            {
               // Unfortunately, it looks like the simplest way to have this work for various
               // slip systems is by passing in the SlipGeom's params in twice...
               SlipGeom slipgeom(parsIt);
               parsIt += SlipGeom::nParams;
               const double* mref = slipgeom.getM();
               const double* sref = slipgeom.getS();
               // our forest interaction matrix has the following calculation:
               // A^{\alpha\beta} = 1/2 * (|m^alpha \cdot s^alpha| + |m^alpha \cdot (m^beta \cross s^beta)|)
               RAJA::View<const double, RAJA::Layout<2> > mView(mref, SlipGeom::nslip, ecmech::ndim);
               RAJA::View<const double, RAJA::Layout<2> > sView(sref, SlipGeom::nslip, ecmech::ndim);
               RAJA::View<double, RAJA::Layout<2> > aView(&m_a_mat[0], SlipGeom::nslip, SlipGeom::nslip);

               for (int alpha = 0; alpha < SlipGeom::nslip; alpha++) {
                  for (int beta = 0; beta < SlipGeom::nslip; beta++) {
#ifndef ORO_USE_INTERMAT
                     const double mds = mView(alpha, 0) * sView(beta, 0) +
                                        mView(alpha, 1) * sView(beta, 1) +
                                        mView(alpha, 2) * sView(beta, 2);
                     const double mdmxs = mView(alpha, 0) * (mView(beta, 1) * sView(beta, 2) - mView(beta, 2) * sView(beta, 1)) +
                                          mView(alpha, 1) * (mView(beta, 2) * sView(beta, 0) - mView(beta, 0) * sView(beta, 2)) +
                                          mView(alpha, 2) * (mView(beta, 0) * sView(beta, 1) - mView(beta, 1) * sView(beta, 0));
                     aView(alpha, beta) = 1.0 / 2.0 * (std::abs(mds) + std::abs(mdmxs));
#else
                     // use interaction matrix
                     aView(alpha, beta) = m_inter_mat[alpha * SlipGeom::nslip + beta];
#endif
                  }
               }
            }

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
          *
          * @note Unlike setParams, this does **not** append a trailing copy of the slip
          * geometry's own parameters (setParams consumes `SlipGeom::nParams` extra
          * values at the end to reconstruct a temporary `SlipGeom` for the
          * forest-interaction matrix; see the file-level documentation). So the number
          * of values pushed here is `#nParams - SlipGeom::nParams` (for the `nVPer ==
          * 1` case that all current `cases/` configurations use), which only equals
          * `#nParams` when `SlipGeom::nParams == 0`. For a slip geometry with its own
          * parameters -- e.g. `SlipGeomBCCNonSchmid` (`nParams == 3`), used by
          * `Kin_OroD_Aniso_BCC_NS`/`"evptn_BCC_E"` -- this mismatch means the
          * `ECMECH_DEBUG` assertion below will actually fail, and a round trip through
          * `setParams(getParams())` will not reproduce the original parameters. Not
          * fixed here since this is a documentation-only pass; flagged for follow-up.
          * @param[in,out] params Parameter vector to append to (not cleared first).
          */
         __ecmech_host__
         void getParams(std::vector<double> & params
                        ) const {
#ifdef ECMECH_DEBUG
            // do not clear params in case adding to an existing set
            int paramsStart = params.size();
#endif
            params.push_back(m_mu_ref);
            params.push_back(m_tkelv_ref);
            for (int iVal = 0; iVal < nVPer; ++iVal) {
               params.push_back(m_berg_mag[iVal]);
            }

            params.push_back(m_lbar_b);
            // phonon drag params
            params.push_back(m_gam_ro);
            params.push_back(m_wrD);

            // thermal activation params
            params.push_back(m_fD);
            for (int iVal = 0; iVal < nVPer; ++iVal) {
               params.push_back(m_c_1[iVal]);
            }

            params.push_back(m_tau_a);
            params.push_back(m_p);
            params.push_back(m_q);
            for (int iVal = 0; iVal < nVPer; ++iVal) {
               params.push_back(m_c_2[iVal]);
            }

            for (int iVal = 0; iVal < nIH; ++iVal) {
               params.push_back(m_inter_mat[iVal]);
            }

            //////////////////////////////
            // Dislocation evolution stuff

            params.push_back(m_c_ann);
            params.push_back(m_d_ann);
            params.push_back(m_c_trap);
            params.push_back(m_c_mult);

            //////////////////////////////
            // Dislocation Densities
            for (int iVal = 0; iVal < SlipGeom::nslip; iVal++) {
               params.push_back(m_qM[iVal]);
            }

            for (int iVal = 0; iVal < SlipGeom::nslip; iVal++) {
               params.push_back(m_qT[iVal]);
            }

            //////////////////////////////
#ifdef ECMECH_DEBUG
            assert((params.size() - paramsStart) == nParams);
#endif
         }

         /**
          * @brief Describe the per-slip-system mobile and total dislocation density
          * history variables ("rho_dd_mobile_N", "rho_dd_total_N") for history-array
          * bookkeeping.
          * @param[out] names Appended with `"rho_dd_mobile_" + iSlip` for each slip
          * system, then `"rho_dd_total_" + iSlip` for each slip system.
          * @param[out] init Appended with the corresponding initial density (#m_qM /
          * #m_qT).
          * @param[out] plot Appended with `true` for each entry.
          * @param[out] state Appended with `true` for each entry.
          */
         __ecmech_host__
         void getHistInfo(std::vector<std::string> & names,
                          std::vector<double>       & init,
                          std::vector<bool>        & plot,
                          std::vector<bool>        & state) const {
            for (int iSlip = 0; iSlip < SlipGeom::nslip; iSlip++) {
               names.push_back("rho_dd_mobile_" + std::to_string(iSlip));
               init.push_back(m_qM[iSlip]);
               plot.push_back(true);
               state.push_back(true);
            }

            for (int iSlip = 0; iSlip < SlipGeom::nslip; iSlip++) {
               names.push_back("rho_dd_total_" + std::to_string(iSlip));
               init.push_back(m_qT[iSlip]);
               plot.push_back(true);
               state.push_back(true);
            }
         }

      private:

         /** @brief Number of slip systems handled by this instance. */
         const int nslip; // could template on this if there were call to do so

         //////////////////////////////
         // MTS-like stuff

         // parameters
         /** @brief Mean forest spacing over Burgers vector, L̄/b, used in the thermal-activation reference rate γ̇_w (see evalGdot). Shared across all slip systems (not yet made per-group). */
         double m_lbar_b; // We might need to make this per SS as well
         /** @brief Shear modulus at reference conditions, μ_ref [stress units]. */
         double m_mu_ref;
         /** @brief Reference temperature, T_ref [Kelvin]. */
         double m_tkelv_ref;
         /** @brief Thermal attempt-frequency factor f_D used in the thermal-activation reference rate γ̇_w. */
         double m_fD;
         /** @brief Burgers vector magnitude per group, bᵢ. */
         double m_berg_mag[nVPer];
         /** @brief Per-group thermal-activation energy scale numerator C1ᵢ. */
         double m_c_1[nVPer];
         /** @brief Reference stress τ_a; plays the role of either the athermal stress floor or the MTS-normalizing stress in the slip-rate law, depending on `withGAthermal` (same role-swap as KineticsKMBalD). */
         double m_tau_a;
         /** @brief Per-group Taylor hardening coefficient c_2ᵢ, scaling the forest-density term into a CRSS contribution (see getVals). */
         double m_c_2[nVPer];
         /** @brief MTS activation-energy exponent p; only used if `pOne` is false (otherwise p is taken to be exactly 1). */
         double m_p; // only used if pOne is false
         /** @brief MTS activation-energy exponent q; only used if `qOne` is false (otherwise q is taken to be exactly 1). */
         double m_q; // only used if qOne is false
         /** @brief User-supplied forest-interaction values, used directly (unconditionally, regardless of `ORO_USE_INTERMAT`) by getVals to compute the CRSS forest-hardening contribution; a single scalar if `isotropic`, else one entry per slip-system pair. Distinct from #m_a_mat, the (usually geometrically-computed) matrix that getSdotN uses for the dislocation-density evolution ODE -- see the file-level documentation for how the two relate. */
         double m_inter_mat[nIH]; // symmetric matrix

         /** @brief Reference-rate coefficient γ̇_r0 for the drag-limited branch. */
         double m_gam_ro;
         /** @brief Drag stress scale w_rD. */
         double m_wrD;

         // derived from parameters
         /** @brief Per-group overflow/underflow stress-ratio thresholds (mirroring #gam_ratio_min, #gam_ratio_ovf) and power-law exponent helpers `xnn = 1/xm`, `xn = xnn - 1`, where the effective rate-sensitivity exponent xm is derived from the reference MTS parameters (same construction as KineticsKMBalD). */
         double m_t_max[nVPer], m_t_min[nVPer], m_xn[nVPer], m_xnn[nVPer];

         //////////////////////////////
         // Dislocation evolution stuff

         /** @brief Dislocation annihilation-rate coefficient c_ann. */
         double m_c_ann;
         /** @brief Annihilation capture distance d_ann. */
         double m_d_ann;
         /** @brief Mobile-to-forest trapping-rate coefficient c_trap. */
         double m_c_trap;
         /** @brief Dislocation multiplication-rate coefficient c_mult. */
         double m_c_mult;
         // stored c-style
         /** @brief Forest-interaction matrix A^αβ used by getSdotN for the dislocation-density evolution ODE. Computed geometrically at setParams time from the slip system m/s vectors (see the file-level documentation) unless `ORO_USE_INTERMAT` is defined, in which case it's copied from #m_inter_mat instead. Distinct from #m_inter_mat, which getVals always reads directly for the CRSS calculation regardless of this matrix. */
         double m_a_mat[SlipGeom::nslip * SlipGeom::nslip]; // Forest interaction matrix

         //////////////////////////////
         // Initial dislocation densities
         // so _hdn_init in other models

         /** @brief Initial mobile dislocation density per slip system. */
         double m_qM[SlipGeom::nslip];
         /** @brief Initial total dislocation density per slip system. */
         double m_qT[SlipGeom::nslip];
         /** @brief Floor on both dislocation densities (`1e-4` times the smallest initial mobile density). */
         double m_hdn_min;

      public:

         /**
          * @brief Reference slip rate used for scaling elsewhere in the solve (e.g. by
          * the evptn elastic-strain/rotation solver): the maximum, over all slip
          * systems, of the per-slip-system thermal-activation-plus-phonon reference
          * rate computed by getVals.
          * @param vals Kinetic values from getVals; `vals[0]` is the precomputed maximum
          * rate.
          * @return `vals[0]`.
          */
         __ecmech_hdev__
         inline
         double
         getFixedRefRate(const double* const vals) const
         {
            return vals[0];
         }

         /**
          * @brief Precompute the kinetic values used by evalGdots: the per-slip-system
          * CRSS from forest hardening, the mobile density, the thermal energy scale,
          * and a representative reference rate -- see the "Kinetic values" equations in
          * the file-level documentation in ECMech_kinetics_OrowanD.h.
          *
          * Could eventually bring in additional pressure and temperature dependence
          * through the dependence of `m_mu_ref` on such conditions.
          * @param[out] vals Kinetic values, length #nVals: `vals[0]` = max reference
          * rate, `vals[1+i]` = ĝᵢ (CRSS), `vals[1+nslip+i]` = qMᵢ (mobile density),
          * `vals[1+2*nslip+i]` = c_tᵢ.
          * @param p Pressure; not currently used by this model.
          * @param tkelv Temperature [Kelvin], T.
          * @param[in] h_state Current hardening state: `h_state[0..nslip-1]` mobile
          * densities qM, `h_state[nslip..2*nslip-1]` total densities qT.
          * @return The average CRSS ĝ across all slip systems.
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
            double const nVPerInv = 1.0 / nslip;

            double maxRefRate = 0.0;
            double hdnScale = 0.;
            for (int iVal = 0; iVal < nslip; ++iVal) {
               const double int_q = isotropic ? sqrt(m_inter_mat[0] * vecsssumabs<SlipGeom::nslip>(&h_state[nslip])) :
                                                sqrt(vecsyadotb<SlipGeom::nslip>(&m_inter_mat[iVal * nslip], &h_state[nslip]));
               const double hdnI = perSS ? (m_c_2[iVal] * int_q) : (m_c_2[0] * int_q);
               hdnScale += hdnI;
               vals[1 + iVal] = hdnI;
               vals[1 + nslip + iVal] = h_state[iVal];
               // Thermal activation + phonon ref slip rate = (1/(f_D * \bar{L}/b * sqrt(qM_0)/sqrt(qM)) + 1/(gammadot_r0 * qM))^-1
               const double isqrth = 1.0 / sqrt(vals[1 + nslip + iVal]);
               const double rate = 1.0 / ((1.0 / (m_lbar_b * m_fD * isqrth)) + (1.0 / (m_gam_ro * vals[1 + nslip + iVal])));
               if (rate > maxRefRate) {
                  maxRefRate = rate;
               }
            }

            // average flow strength across all slip systems
            hdnScale = hdnScale * nVPerInv;
            vals[0] = maxRefRate;

            for (int iVal = 0; iVal < nVPer; ++iVal) {
               vals[1 + 2 * nslip + iVal] = m_c_1[iVal] / tkelv; // _c_t
            }

            return hdnScale;
         }

         /**
          * @brief Evaluate the slip rate and its derivative w.r.t. resolved shear stress
          * on every slip system.
          * @param[out] gdot Slip rate on each slip system γ̇ [1/time].
          * @param[out] dgdot_dtau Derivative of slip rate w.r.t. resolved shear stress on
          * each slip system.
          * @param[in] tau Resolved shear stress on each slip system τ [stress units].
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
          * @brief Evaluate the MTS activation-energy function E(t) = -c_e·[1 -
          * sign(t)·|t|^p]^q and its derivative factor, for a single (signed)
          * dimensionless MTS argument -- see the "MTS thermal-activation energy
          * function" section in the file-level documentation in
          * ECMech_kinetics_OrowanD.h.
          *
          * Handles three regimes: `t` near zero (linearizes to avoid a 0/0 in the
          * derivative when `pOne` is false), `q_arg = 1 - p_func` at or below zero
          * (barrier fully overcome -- "pegged" to `E = 0`), and the general case.
          * @param[out] exp_arg E(t), the log of the thermal-activation rate factor.
          * @param[out] mts_dfac Derivative helper factor; combined with the caller's own
          * chain-rule terms to get `d(exp(E))/dτ`.
          * @param c_e Thermal energy prefactor c_e (= c_t·μ).
          * @param denom_i 1/g_MTS, the reciprocal of whichever stress is currently
          * playing the MTS-normalizing role.
          * @param t_frac The dimensionless MTS argument t.
          */
         __ecmech_hdev__
         inline
         void
         get_mts_dG(double &exp_arg,
                    double &mts_dfac,
                    const double c_e, const double denom_i,
                    const double t_frac) const {
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
          * @brief Evaluate the balanced thermally-activated + power-law-tail slip rate,
          * combined by harmonic mean with the drag-limited rate, for a single slip
          * system -- see the "Slip-rate law" equations and symbol table in the
          * file-level documentation in ECMech_kinetics_OrowanD.h (`gam_w`/`gam_r` here
          * are γ̇_w/γ̇_r, computed by getVals from this slip system's own mobile
          * density).
          * @param[out] gdot Slip rate γ̇ [1/time].
          * @param[out] l_act `true` if the slip system is active.
          * @param[out] dgdot_dtau Derivative of `gdot` w.r.t. resolved shear stress.
          * @param[in] vals Kinetic values from getVals.
          * @param iSlip Slip system index (selects which parameter group to use when
          * `perSS`).
          * @param tau Resolved shear stress τ [stress units].
          * @param mu Shear modulus, used to form the thermal energy prefactor c_e.
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

            const double gIn = vals[1 + iSlip];
            const double qm = vals[1 + nslip + iSlip];
            const double xn = perSS ? m_xn[iSlip] : m_xn[0];
            const double xnn = perSS ? m_xnn[iSlip] : m_xnn[0];
            const double t_min = perSS ? m_t_min[iSlip] : m_t_min[0];
            const double t_max = perSS ? m_t_max[iSlip] : m_t_max[0];
            const double c_t = perSS ? vals[1 + 2 * nslip + iSlip] : vals[1 + 2 * nslip];
            const double gam_w = m_lbar_b * m_fD / sqrt(qm);
            const double gam_r = m_gam_ro * qm;

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
            double gdot_r, dgdot_r_dtau;
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
               dgdot_r_dtau = gam_r * temp / m_wrD;
            }
            //
            if (at_0 > t_max) {
               // have overflow of thermally activated kinetics, purely drag limited

               gdot = gdot_r;

               dgdot_dtau = dgdot_r_dtau;
               gdot = copysign(gdot, tau);

               l_act = true;
               return;
            }

            double gdot_w, dgdot_w_dtau;
            double dgdot_w_dg; // only used if !withGAthermal
            //
            // calculate thermally activated kinetics
            {
               double c_e = c_t * mu;
               //
               const double t_frac = (fabs(tau) - gAth) * g_i;
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
               dgdot_w_dtau = mts_dfac * gdot_w;
               if (!withGAthermal) {
                  dgdot_w_dg = dgdot_w_dtau * t_frac;
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
                  dgdot_w_dtau = dgdot_w_dtau - contrib; // sign used to be the other way, but suspect that was a bug
                  if (!withGAthermal) {
                     dgdot_w_dg = dgdot_w_dg - contrib * t_frac_m;
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
               dgdot_w_dtau = dgdot_w_dtau + contrib;
               if (!withGAthermal) {
                  dgdot_w_dg = dgdot_w_dg + contrib * at_0;
               }
            }

            l_act = true;
            //
            {
               // gdot = gdot_r;
               gdot = one / (one / gdot_w + one / gdot_r);
               const double gdrdiv2 = one / (gdot_r * gdot_r);
               const double gdwdiv2 = one / (gdot_w * gdot_w);

               // dgdot_dtau = dgdot_r_dtau;
               dgdot_dtau = (gdot * gdot) * (dgdot_w_dtau * gdwdiv2 + dgdot_r_dtau * gdrdiv2);
               //
            }

            gdot = copysign(gdot, tau);
         } // evalGdot

         /**
          * @brief Advance the per-slip-system mobile and total dislocation densities
          * one time step.
          *
          * First computes each slip system's dislocation glide velocity via the Orowan
          * relation, `nu[i] = |gdot[i]| / (max(hs_o[i], m_hdn_min) · bᵢ)`, then drives
          * the shared vector-hardness SNLS solve (#updateHN, with `relaxed_solver =
          * true`).
          *
          * If `LOGFORM`, the solve is done in log space (both densities floored at
          * #m_hdn_min before taking the log), guaranteeing positivity by construction,
          * and the result is exponentiated back. Otherwise, the raw (non-log) solve has
          * no such guarantee, so the result is checked for negative densities; if any
          * are found, the step is retried by substepping: `dt` is split into 10 equal
          * substeps, and each substep re-derives `nu` from the *original* `gdot` (the
          * slip rate is assumed constant across the whole time step, per the comment
          * below) but the *current* (most recently updated) density estimate, then
          * re-solves for that substep. If densities are still negative after all 10
          * substeps, this fails hard via `ECMECH_FAIL` (this substep fallback is
          * acknowledged in the code as "pretty ad-hoc" but reported to work well enough
          * in practice for simple test cases).
          * @param[out] hs_u Updated (end-of-step) densities [2*SlipGeom::nslip]: mobile
          * then total.
          * @param[in] hs_o Beginning-of-step densities [2*SlipGeom::nslip]; floored at
          * #m_hdn_min.
          * @param dt Time step size.
          * @param[in] gdot Beginning-of-step slip rates on all slip systems.
          * @param[in] hvals Additional per-slip-system values, forwarded to getSdotN
          * (currently unused there).
          * @param tkelv Temperature [Kelvin].
          * @param outputLevel Verbosity passed through to the SNLS solver.
          * @return Function-evaluation count (summed across any substeps), or a
          * negative value if the (non-substepped) solve failed to converge.
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

            double ihs_o[SlipGeom::nslip * 2];
            double nu[SlipGeom::nslip];
            for (int i = 0; i < nslip * 2; i++) {
               if (i < nslip) {
                  const double div = perSS ? fmax(hs_o[i], m_hdn_min) * m_berg_mag[i] :
                                     fmax(hs_o[i], m_hdn_min) * m_berg_mag[0];
                  nu[i] = abs(gdot[i]) / (div);
               }
               ihs_o[i] = fmax(hs_o[i], m_hdn_min);
               if (LOGFORM) {
                  ihs_o[i] = log(ihs_o[i]);
               }
            }

            int nFEvals = updateHN<KineticsOrowanD>(this,
                                                   &hs_u[0], &ihs_o[0], dt, nu, hvals, tkelv,
                                                   outputLevel);
            if (LOGFORM) {
               for (int i = 0; i < nslip * 2; i++) {
                  hs_u[i] = exp(hs_u[i]);
               }
            }
            else
            {
               // We need to check that none of our solutions became negative
               // If we did obtain something negative then we should abort
               // It means our time step was too large for this step.
               // If this is not desirable / possible then we should probably
               // do a terrible hack and cut the dt by some factor resolve things by
               // assuming a constant slip rate during the time step, and then
               // evolve the dd content. We would get a solution, but it wouldn't necessarily
               // be correct.
               bool flag = false;
               for (int i = 0; i < 2 * nslip; i++) {
                  if(hs_u[i] < zero) {
                     flag = true;
                     break;
                  }
               }
               if (flag)
               {
                  ECMECH_WARN(__func__, "Solver returned negative dislocation values trying again by substepping through the solution");
                  // This is pretty ad-hoc but it seems to work fairly well for a number of simple test cases.
                  // It's definitely not the best way to probably do things though...
                  const double dtnew = dt / 10.0;
                  double hs_temp[2 * SlipGeom::nslip];

                  for (int iSlip = 0; iSlip < 2 * SlipGeom::nslip; iSlip++) {
                     hs_u[iSlip] = fmax(hs_o[iSlip], m_hdn_min);
                  }

                  for (int i = 0; i < 10; i++)
                  {
                     for (int iSlip = 0; iSlip < 2 * SlipGeom::nslip; iSlip++) {
                        hs_temp[iSlip] = fmax(hs_u[iSlip], m_hdn_min);
                        if (iSlip < nslip)
                        {
                           const double div = perSS ? fmax(hs_temp[iSlip], m_hdn_min) * m_berg_mag[iSlip] :
                           fmax(hs_temp[iSlip], m_hdn_min) * m_berg_mag[0];
                           nu[iSlip] = abs(gdot[iSlip]) / (div);
                        }
                     }
                     nFEvals += updateHN<KineticsOrowanD>(this,
                                                         &hs_u[0], hs_temp, dtnew, nu, hvals, tkelv,
                                                         outputLevel);
                     flag = false;
                     for (int iSlip = 0; iSlip < 2 * nslip; iSlip++) {
                        if(hs_u[iSlip] < zero) {
                           flag = true;
                           break;
                        }
                     }
                  }

                  if (flag)
                  {
                     for (int iSlip = 0; iSlip < 2 * SlipGeom::nslip; iSlip++) {
                        printf("dd[%d]: %lf ", iSlip, hs_u[iSlip]);
                     }
                     printf("\n");
                     ECMECH_FAIL(__func__, "Solver returned negative dislocation values!");
                  }
               }
            }

            return nFEvals;
         }

         /**
          * @brief Pass the per-slip-system dislocation velocities through unchanged, as
          * the "evolution values" getSdotN needs.
          *
          * Unlike KineticsBCCMD's `getEvolVals` (which accumulates into its last
          * element and therefore requires the caller to zero-initialize `evolVals`
          * first), this is a plain copy with no accumulation, so it's safe regardless of
          * how the caller's buffer was initialized.
          * @param[out] evolVals Copy of `nu`, length #nEvolVals.
          * @param[in] nu Per-slip-system dislocation velocities, as computed by
          * `updateH` via the Orowan relation.
          */
         __ecmech_hdev__
         inline
         void
         getEvolVals(double* const evolVals,
                     const double* const nu
                     ) const
         {
            // We're not really evolving anything here at this point in time
            // so we can just return...
            // The Jacobian doesn't take in gdots / nu so we're left with this.
            for (int i = 0; i < SlipGeom::nslip; i++) {
               evolVals[i] = nu[i];
            }
         }

         /**
          * @brief Evaluate the mobile/total dislocation-density evolution rates
          * (multiplication, trapping, annihilation) and their Jacobian, in log space if
          * `LOGFORM`, for the vector-hardness SNLS solve (#updateHN) -- see the
          * "Hardening law" equations and symbol table in the file-level documentation
          * in ECMech_kinetics_OrowanD.h.
          *
          * `h` holds the mobile densities in `h[0..nslip-1]` and total densities in
          * `h[nslip..2*nslip-1]` (exponentiated from `h_i` first if `LOGFORM`). The
          * per-slip-system forest density (`forest_dis`) is the forest-interaction
          * matrix `m_a_mat` applied to the total-density half of `h`.
          * @param[out] sdot Density evolution rate for each of the 2*nslip components
          * (mobile then total), in log space if `LOGFORM`.
          * @param[out] dsdot_ds Jacobian of `sdot` w.r.t. `h`
          * [nDimSys x nDimSys, nDimSys = 2*nslip], row-major via #ECMECH_NN_INDX; left
          * untouched if `nullptr`.
          * @param[in] h_i Current hardening state (mobile then total densities), in log
          * space if `LOGFORM`.
          * @param[in] evolVals Per-slip-system dislocation velocities from
          * getEvolVals.
          * @param hvals Unused by this model.
          * @param tkelv Temperature [Kelvin]; not currently used by this model.
          */
         __ecmech_hdev__
         inline
         void
         getSdotN( double* sdot,
                   double* dsdot_ds,
                   const double* const h_i,
                   const double* const evolVals,
                   const double* const /*hvals*/,
                   double /*tkelv*/                   ) const
         {
            // Hopefully, the compiler is pretty smart here and is able to optimize these
            // loops as if we're using the templated values. Since, this is essentially
            // a constexpr for these variables.
            const int nslip = SlipGeom::nslip;
            const int JDIM = 2;
            const int nDimSys = 2 * SlipGeom::nslip;

            double forest_dis[nslip];
            constexpr int h_content = (LOGFORM) ? 2 * SlipGeom::nslip : 1;
            double hexp[h_content];
            if (LOGFORM) {
               for (int iDD = 0; iDD < 2 * nslip; iDD++) {
                  hexp[iDD] = exp(h_i[iDD]);
               }
            }
            const double* const h = (LOGFORM) ? &hexp[0] : h_i;
            vecsVMa<SlipGeom::nslip>(&forest_dis[0], &m_a_mat[0], &h[nslip]);

            for (int iM = 0; iM < nslip; iM++) {
               const double sqrt_fD = sqrt(forest_dis[iM]);
               const double q_dmult = m_c_mult * sqrt_fD * h[iM] * evolVals[iM];
               const double q_dtrap = m_c_trap * sqrt_fD * h[iM] * evolVals[iM];
               // This could become a very large number and could become problematic
               // later on. Do we want to cap it at some large value?
               // Although, it might be that this is only a problem if q and qM are defined
               // with units 1/m^2 rather than 1/mm^2 or 1/micron^2
               const double q_dann = m_c_ann * m_d_ann * h[iM] * h[iM] * evolVals[iM];
               // mobile dislocation density rate of change
               sdot[iM] = q_dmult - q_dtrap - q_dann;
               // total dislocation density rate of change
               sdot[iM + nslip] = q_dmult - q_dann;
            }

            if (LOGFORM) {
               for (int iDD = 0; iDD < 2 * nslip; iDD++) {
                  sdot[iDD] *= (1.0 / h[iDD]);
               }
            }
            // The dsdot_ds calculation for our nonlinear solve
            if (dsdot_ds) {
               // zero out dsdot_ds matrix
               for (int i = 0; i < nDimSys * nDimSys; i++) {
                  dsdot_ds[i] = ecmech::zero;
               }

               RAJA::View<double, RAJA::Layout<JDIM> > dsdot_ds_view(dsdot_ds, nDimSys, nDimSys);
               // dqM/dqM portion of dsdot_ds
               for (int iM = 0; iM < nslip; iM++) {
                  const double sqrt_fD = sqrt(forest_dis[iM]);
                  const double q_dmult_dtrap = (m_c_mult - m_c_trap) * sqrt_fD;
                  // Although, it might be that this is only a problem if q and qM are defined
                  // with units 1/m^2 rather than 1/mm^2 or 1/micron^2
                  const double q_dann = 2 * m_c_ann * m_d_ann * h[iM];
                  dsdot_ds_view(iM, iM) = evolVals[iM] * (q_dmult_dtrap - q_dann);
               }

               // dq/dqM portion of dsdot_ds
               for (int iT = 0; iT < nslip; iT++) {
                  const double sqrt_fD = sqrt(forest_dis[iT]);
                  const double q_dmult = m_c_mult * sqrt_fD;
                  // This could become a very large number and could become problematic
                  // later on. Do we want to cap it at some large value?
                  // Although, it might be that this is only a problem if q and qM are defined
                  // with units 1/m^2 rather than 1/mm^2 or 1/micron^2
                  const double q_dann = 2 * m_c_ann * m_d_ann * h[iT];
                  dsdot_ds_view(iT + nslip, iT) =  evolVals[iT] * (q_dmult - q_dann);
               }

               RAJA::View<const double, RAJA::Layout<JDIM> > amat(&m_a_mat[0], nslip, nslip);
               // dq/dq portion of dsdot_ds
               for (int iT = 0; iT < nslip; iT++) {
                  for (int jT = 0; jT < nslip; jT++) {
                     // First, terms found only on the diagonal of this submatrix
                     const double ifact = ecmech::onehalf / sqrt(forest_dis[iT]);
                     const double q_dmult = m_c_mult * amat(iT, jT) * ifact;

                     dsdot_ds_view(iT + nslip, jT + nslip) = h[iT] * evolVals[iT] * q_dmult;
                  }
               }

               // dqM/dq portion of dsdot_dt
               for (int iT = 0; iT < nslip; iT++) {
                  for (int jT = 0; jT < nslip; jT++) {
                     const double ifact = ecmech::onehalf / sqrt(forest_dis[iT]);
                     const double q_dmult_dtrap = (m_c_mult - m_c_trap) * amat(iT, jT) * ifact;

                     dsdot_ds_view(iT, jT + nslip) = h[iT] * evolVals[iT] * (q_dmult_dtrap);
                  }
               }

               if (LOGFORM) {
                  for (int iDD = 0; iDD < 2 * nslip; iDD++) {
                     dsdot_ds_view(iDD, iDD) -= sdot[iDD];
                  }
                  for (int iDD = 0; iDD < 2 * nslip; iDD++) {
                     for (int jDD = 0; jDD < 2 * nslip; jDD++) {
                        dsdot_ds_view(iDD, jDD) *= h[jDD] / h[iDD];
                     }
                  }
               }
            } // if dsdot_ds
         }
   }; // class KineticsOrowanD
} // namespace ecmech

#endif // ECMECH_KINETICS_OROWAND_H
