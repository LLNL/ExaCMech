/**
 * @file ECMech_kinetics.h
 * @brief Kinetics model interface contract, plus the shared nonlinear-solve
 * infrastructure used to update a model's hardening state variable(s) each time step.
 *
 * A "kinetics model" bundles two physically distinct things into one class:
 * 1. **Slip kinetics**: the slip-rate law, i.e. how fast each slip system shears (γ̇)
 *    given the resolved shear stress (τ) and the current flow strength/hardening state.
 * 2. **Hardness evolution kinetics**: the ODE governing how the hardening state
 *    variable(s) (CRSS, dislocation density, etc.) evolve with time and slip rate.
 *
 * These two could be separated to increase modularity, but keeping them grouped in a
 * single class makes it easier to keep their interactions both flexible and simple to
 * manage. The concrete model classes (`kinetics/ECMech_kinetics_VocePL.h`,
 * `_BCCMD.h`, `_KMBalD.h`, `_OrowanD.h`) are `#include`d at the bottom of this file, and
 * are used as template parameters wherever `evptn::matModel` (ECMech_evptnWrap.h) is
 * instantiated (see `cases/ECMech_cases_*_defs.h`).
 *
 * **Required traits** (every kinetics model must provide these as `static constexpr
 * int`s):
 * - `nH`: number of hardening state variables
 * - `nParams`: number of parameters needed to construct/configure the model
 * - `nVals`: number of "kinetic values" precomputed once per point by `getVals` and
 *   reused across all slip systems by `evalGdots`
 * - `nEvolVals`: number of intermediate values precomputed once per point by
 *   `getEvolVals` and reused by the hardening-state solve
 *
 * **Required member functions**:
 *
 * - `constructor(int nslip)`
 *   Construct with a given number of slip systems.
 *
 * - `void setParams(const std::vector<double> & params)`
 *   Configure from a flat parameter vector (order is model-specific; see each concrete
 *   model).
 *
 * - `void getParams(std::vector<double> & params) const`
 *   Append the model's current parameters to `params`, in the same order `setParams`
 *   expects them.
 *
 * - `void getHistInfo(names, init, plot, state) const`
 *   Describe each hardening history variable (name, initial value, whether it should be
 *   plotted, whether it's a state variable) for history-array bookkeeping.
 *
 * - `double getVals(double* const vals, double p, double tkelv, const double* const h_state) const`
 *   Precompute the `nVals` kinetic quantities (e.g. reference slip rate, CRSS) from
 *   pressure `p`, temperature `tkelv`, and hardening state `h_state`; returns a
 *   representative/average flow-strength-like value.
 *
 * - `void evalGdots(gdot, dgdot_dtau, tau, vals) const`
 *   Evaluate the slip rate and its derivative w.r.t. resolved shear stress on every slip
 *   system, given the resolved shear stresses `tau` and the `vals` from `getVals`.
 *
 * - `int updateH(hs_u, hs_o, dt, gdot, hvals, tkelv, outputLevel = 0) const`
 *   Advance the hardening state from `hs_o` (start of step) to `hs_u` (end of step) over
 *   time step `dt`, given the slip rates `gdot`; `hvals` holds any extra variables the
 *   internal solver needs. Returns the nonlinear solver's function-evaluation count, or
 *   a **negative value to signal a failed solve** (rather than throwing, since
 *   exceptions are unavailable on the GPU — see ECMech_port.h).
 *
 * - `void getEvolVals(evolVals, gdot) const`
 *   Precompute the `nEvolVals` values the hardening-rate evaluation needs from the
 *   current slip rates `gdot`.
 *
 * - Exactly one of:
 *
 *   - `void getSdot1(double &sdot, double &dsdot_ds, double h, const double* const evolVals, double tkelv) const`
 *     For models with a single hardening state variable (paired with #updateH1). Both
 *     `sdot` (hardening rate) and `dsdot_ds` (its derivative w.r.t. `h`) are taken by
 *     reference and are always computed; the caller decides whether to use `dsdot_ds`.
 *
 *   - `void getSdotN(double* sdot, double* dsdot_ds, const double* const h, const double* const evolVals, const double* const hvals, double tkelv) const`
 *     For models with multiple hardening state variables (paired with #updateHN). Here
 *     `dsdot_ds` (the Jacobian of the rates w.r.t. `h`) is a pointer that may be
 *     `nullptr` when the caller doesn't need it, so implementations must check for
 *     `nullptr` and skip the Jacobian calculation in that case.
 *
 * @see updateH1, updateHN for the shared trust-region-dogleg solves that a model's
 * `updateH` typically delegates to
 */

// -*-c++-*-

#ifndef ECMECH_kinetics_H
#define ECMECH_kinetics_H

#include "SNLS_TrDLDenseG.h"

#include "ECMech_core.h"
#include "ECMech_util.h"

namespace ecmech {
   /**
    * @brief SNLS problem formulation for updating a single scalar hardening state
    * variable implicitly over a time step.
    *
    * Solves `h_delta - sdot(h_o + h_delta) * dt = 0` for `h_delta` using SNLS's
    * trust-region dogleg solver (see #updateH1, which drives this), where `sdot` is
    * supplied by the kinetics model's `getSdot1`. The unknown is scaled by
    * `m_x_scale = max(h_o, 1)` so that the solver operates on an O(1) quantity
    * regardless of the hardening variable's absolute magnitude.
    *
    * @tparam Kinetics Concrete kinetics model type providing `getSdot1` and `nH == 1`.
    */
   template<class Kinetics>
   class Kinetics_H1Problem
   {
      public:
         /** @brief System dimension for the SNLS solve; always 1 for this scalar-hardness problem (mirrors `Kinetics::nH`). */
         static constexpr int nDimSys = Kinetics::nH;

         /**
          * @brief Construct the scalar hardness-update subproblem for one time step.
          * @param kinetics Kinetics model providing `getSdot1`.
          * @param h_o Hardening state value at the start of the step.
          * @param dt Time step size.
          * @param evolVals Precomputed values from the model's `getEvolVals`, forwarded
          * to `getSdot1` unchanged.
          * @param tkelv Temperature [Kelvin].
          */
         __ecmech_hdev__
         Kinetics_H1Problem(const Kinetics* const kinetics,
                            double h_o,
                            double dt,
                            const double* const evolVals,
                            double tkelv) :
            m_kinetics(kinetics), m_h_o(h_o), m_dt(dt), m_evolVals(evolVals), m_tkelv(tkelv)
         {
            m_x_scale = fmax(m_h_o, 1.0); // TO_DO -- generalize this to not max with 1
            // NOTE : see comment below about changing Jacobian calculation if m_res_scale != one / s_scale
            m_res_scale = one / m_x_scale;
         }

         /** @brief Destructor (default; no owned resources). */
         __ecmech_hdev__
         ~Kinetics_H1Problem() {}

         /**
          * @brief Recover the physical (unscaled) hardening state from the solver's
          * scaled solution vector.
          * @param x Solver's scaled solution vector, `x[0]`.
          * @return `h_o + x[0] * m_x_scale`, the end-of-step hardening state.
          */
         __ecmech_hdev__
         inline
         double getHn(const double* const x) const {
            return m_h_o + x[0] * m_x_scale;
         }

         /**
          * @brief SNLS residual/Jacobian evaluation for the implicit backward-Euler
          * update `h_delta = sdot(h_o + h_delta) * dt`.
          * @param[out] resid Scaled residual, `resid[0]`.
          * @param[out] Jacobian Scaled Jacobian, `Jacobian[0]`; left untouched if
          * `nullptr` (Jacobian not requested by the solver).
          * @param[in] x Solver's current scaled solution vector, `x[0]`.
          * @return Always `true` (this problem has no notion of an invalid evaluation).
          */
         __ecmech_hdev__
         inline
         bool computeRJ(double* const resid,
                        double* const Jacobian,
                        const double* const x) {
            bool doComputeJ = (Jacobian != nullptr);

            double h_delta = x[0] * m_x_scale;
            double h = m_h_o + h_delta;

            double sdot, dsdot_ds;
            m_kinetics->getSdot1(sdot, dsdot_ds, h, m_evolVals, m_tkelv);

            resid[0] = (h_delta - sdot * m_dt) * m_res_scale;

            if (doComputeJ) {
               // The below is based on the assumption that m_res_scale = 1/m_x_scale
               // if this were to change in the future version than this would need to become
               // Jacobian[0] = (one - dsdot_ds * m_dt) * m_res_scale * m_x_scale;
               Jacobian[0] = (one - dsdot_ds * m_dt);
            }

            return true;
         } // computeRJ

      private:
         /** @brief Kinetics model supplying the hardening-rate law via `getSdot1`. */
         const Kinetics* m_kinetics;
         /** @brief Hardening state at the start of the step (`m_h_o`) and the time step size (`m_dt`). */
         const double m_h_o, m_dt;
         /** @brief Precomputed values from the model's `getEvolVals`, passed through to `getSdot1`. */
         const double* const m_evolVals;
         /** @brief Temperature [Kelvin], passed through to `getSdot1`. */
         const double m_tkelv;
         /** @brief Scaling factors so the solver operates on an O(1) unknown/residual: `m_x_scale = max(h_o, 1)`, `m_res_scale = 1/m_x_scale`. */
         double m_x_scale, m_res_scale;
   }; // class Kinetics_H1Problem

   /**
    * @brief Advance a single scalar hardening state variable one time step, for kinetics
    * models with `nH == 1` that implement `getSdot1` (e.g. KineticsVocePL,
    * KineticsKMBalD).
    *
    * Builds a Kinetics_H1Problem and drives it to convergence with SNLS's trust-region
    * dogleg solver (`snls::SNLSTrDlDenseG`), starting from the scaled solution `x = 0`
    * (i.e. `h_delta = 0`).
    *
    * @tparam Kinetics Concrete kinetics model type; must implement `getSdot1`.
    * @tparam relaxed_solver If `true` and the first solve (tolerance `1e-10`, 100
    * iterations) fails to converge, retry once from scratch with a relaxed tolerance
    * (`1e-9`) before giving up. Used by callers willing to accept a slightly less
    * converged answer rather than fail the time step outright.
    * @param kinetics Kinetics model instance.
    * @param[out] hs_n End-of-step hardening state.
    * @param hs_o Start-of-step hardening state.
    * @param dt Time step size.
    * @param gdot Slip rates on all slip systems, used to compute the evolution inputs
    * via `kinetics->getEvolVals`.
    * @param tkelv Temperature [Kelvin].
    * @param outputLevel Verbosity passed through to the SNLS solver's diagnostic output.
    * @return The number of function evaluations used by the (possibly relaxed) solve, or
    * **-1 if the solve failed to converge** even after the relaxed retry (if enabled) --
    * signaled this way rather than by throwing, since exceptions are unavailable on the
    * GPU.
    */
   template<class Kinetics, bool relaxed_solver = false>
   __ecmech_hdev__
   inline
   int
   updateH1(const Kinetics* const kinetics,
            double &hs_n,
            double hs_o,
            double dt,
            const double* const gdot,
            double tkelv,
            int outputLevel = 0)
   {
      double evolVals[Kinetics::nEvolVals];
      kinetics->getEvolVals(evolVals, gdot);

      Kinetics_H1Problem<Kinetics> prob(kinetics, hs_o, dt, evolVals, tkelv);
      snls::SNLSTrDlDenseG<Kinetics_H1Problem<Kinetics> > solver(prob);

      snls::TrDeltaControl deltaControl;
      deltaControl._deltaInit = 1e0;
      {
         int maxIter = 100;
         double tolerance = 1e-10;
         solver.setupSolver(maxIter, tolerance, &deltaControl, outputLevel);
      }

      for (int iX = 0; iX < prob.nDimSys; ++iX) {
         solver._x[iX] = 0e0;
      }

      snls::SNLSStatus_t status = solver.solve( );

      int nFevals = solver.getNFEvals();
      if (status != snls::converged) {
         snls::SNLSStatus_t status2 = status;
         if constexpr(relaxed_solver) {
            {
               int maxIter = 100;
               double tolerance = 1e-9;
               solver.setupSolver(maxIter, tolerance, &deltaControl, outputLevel);
            }
            for (int iX = 0; iX < prob.nDimSys; ++iX) {
               solver._x[iX] = 0e0;
            }
            status2= solver.solve();
            nFevals = solver.getNFEvals();
         }
         if (status2 != snls::converged) {
            nFevals = -1;
         }
      }

      hs_n = prob.getHn(solver._x);

      return nFevals;
   } // updateH1

   /**
    * @brief SNLS problem formulation for updating a vector of hardening state variables
    * implicitly over a time step.
    *
    * Vector-valued analog of Kinetics_H1Problem: solves `h_delta - sdot(h_o + h_delta) *
    * dt = 0` for `h_delta` using SNLS's trust-region dogleg solver (see #updateHN, which
    * drives this), where `sdot` and its Jacobian are supplied by the kinetics model's
    * `getSdotN`. Each component of the unknown is independently scaled by
    * `m_x_scale[i] = max(h_o[i], 1)`.
    *
    * @tparam Kinetics Concrete kinetics model type providing `getSdotN`.
    */
   template<class Kinetics>
   class Kinetics_HNProblem
   {
      public:
         /** @brief System dimension for the SNLS solve; equal to `Kinetics::nH`, the number of hardening state variables. */
         static constexpr int nDimSys = Kinetics::nH;

         /**
          * @brief Construct the vector-valued hardness-update subproblem for one time
          * step.
          * @param kinetics Kinetics model providing `getSdotN`.
          * @param h_o Hardening state vector at the start of the step [nDimSys].
          * @param dt Time step size.
          * @param evolVals Precomputed values from the model's `getEvolVals`, forwarded
          * to `getSdotN` unchanged.
          * @param hvals Additional per-slip-system values the model's internal solve
          * needs, forwarded to `getSdotN` unchanged.
          * @param tkelv Temperature [Kelvin].
          */
         __ecmech_hdev__
         Kinetics_HNProblem(const Kinetics* const kinetics,
                            const double* const h_o,
                            double dt,
                            const double* const evolVals,
                            const double* const hvals,
                            double tkelv) :
            m_kinetics(kinetics), m_h_o(h_o), m_dt(dt), m_evolVals(evolVals), m_hvals(hvals), m_tkelv(tkelv)
         {
            for (int i = 0; i < nDimSys; i++) {
               m_x_scale[i] = fmax(m_h_o[i], 1.0); // TO_DO -- generalize this to not max with 1
               // NOTE : see comment below about changing Jacobian calculation if m_res_scale != one / s_scale
               m_res_scale[i] = one / m_x_scale[i];
            }
         }

         /** @brief Destructor (default; no owned resources). */
         __ecmech_hdev__
         ~Kinetics_HNProblem() {}

         /**
          * @brief Recover the physical (unscaled) hardening state vector from the
          * solver's scaled solution vector.
          * @param[out] h End-of-step hardening state vector [nDimSys].
          * @param[in] x Solver's scaled solution vector [nDimSys].
          */
         __ecmech_hdev__
         inline
         void getHn(double* h, const double* const x) const {
            for (int i = 0; i < nDimSys; i++) {
               h[i] = m_h_o[i] + x[i] * m_x_scale[i];
            }
         }

         /**
          * @brief SNLS residual/Jacobian evaluation for the implicit backward-Euler
          * update of the hardening state vector.
          * @param[out] resid Scaled residual vector [nDimSys].
          * @param[out] Jacobian Scaled Jacobian matrix [nDimSys x nDimSys], row-major via
          * #ECMECH_NN_INDX; left untouched if `nullptr` (Jacobian not requested by the
          * solver). Note this is filled in two stages: `getSdotN` first writes the raw
          * `d(sdot)/d(h)` Jacobian into this buffer, which is then rescaled in place and
          * has the scaled identity term added to account for the `x -> h` change of
          * variables.
          * @param[in] x Solver's current scaled solution vector [nDimSys].
          * @return Always `true` (this problem has no notion of an invalid evaluation).
          */
         __ecmech_hdev__
         inline
         bool computeRJ(double* const resid,
                        double* const Jacobian,
                        const double* const x) {
            double h[nDimSys];
            for (int i = 0; i < nDimSys; i++) {
               h[i] = m_h_o[i] + x[i] * m_x_scale[i];
            }

            double sdot[nDimSys];
            // The dsdot_ds portion of the Jacobian is set in here if it was provided
            m_kinetics->getSdotN(sdot, Jacobian, h, m_evolVals, m_hvals, m_tkelv);

            for (int i = 0; i < nDimSys; i++) {
               resid[i] = (x[i] * m_x_scale[i] - sdot[i] * m_dt) * m_res_scale[i];
            }

            if (Jacobian) {
               // Multiply dsdot_ds terms by the negative outer product of x_scale and res_scale and dt
               for (int i = 0; i < nDimSys; i++) {
                  for (int j = 0; j < nDimSys; j++) {
                     Jacobian[ECMECH_NN_INDX(i, j, nDimSys)] *= -m_x_scale[j] * m_res_scale[i] * m_dt;
                  }
               }

               // Now add in the identity term
               // The below is based on the assumption that m_res_scale = 1/m_x_scale
               // if this were to change in the future version than this would need to become
               // Jacobian[ECMECH_NN_INDX(i, i, nDimSys)] += ecmech::one * m_x_scale[i] * _r_scale[i]
               for (int i = 0; i < nDimSys; i++) {
                  Jacobian[ECMECH_NN_INDX(i, i, nDimSys)] += ecmech::one;
               }
            } // if Jacobian

            return true;
         } // computeRJ

      private:
         /** @brief Kinetics model supplying the hardening-rate law via `getSdotN`. */
         const Kinetics* m_kinetics;
         /** @brief Hardening state vector at the start of the step [nDimSys]. */
         const double* const m_h_o;
         /** @brief Time step size. */
         const double m_dt;
         /** @brief Precomputed values from the model's `getEvolVals`, passed through to `getSdotN`. */
         const double* const m_evolVals;
         /** @brief Additional per-slip-system values the model's internal solve needs, passed through to `getSdotN`. */
         const double* const m_hvals;
         /** @brief Temperature [Kelvin], passed through to `getSdotN`. */
         const double m_tkelv;
         /** @brief Per-component scaling factors so the solver operates on O(1) unknowns/residuals: `m_x_scale[i] = max(h_o[i], 1)`, `m_res_scale[i] = 1/m_x_scale[i]`. */
         double m_x_scale[nDimSys], m_res_scale[nDimSys];
   }; // class Kinetics_HNProblem

   /**
    * @brief Advance a vector of hardening state variables one time step, for kinetics
    * models with `nH > 1` that implement `getSdotN` (e.g. KineticsOrowanD,
    * KineticsBCCMD).
    *
    * Vector-valued analog of #updateH1: builds a Kinetics_HNProblem and drives it to
    * convergence with SNLS's trust-region dogleg solver (`snls::SNLSTrDlDenseG`),
    * starting from the scaled solution `x = 0` (i.e. `h_delta = 0`).
    *
    * @tparam Kinetics Concrete kinetics model type; must implement `getSdotN`.
    * @tparam relaxed_solver If `true` and the first solve (tolerance `1e-10`, 100
    * iterations) fails to converge, retry once from scratch with a relaxed tolerance
    * (`1e-9`) before giving up.
    * @param kinetics Kinetics model instance.
    * @param[out] hs_n End-of-step hardening state vector [Kinetics::nH].
    * @param hs_o Start-of-step hardening state vector [Kinetics::nH].
    * @param dt Time step size.
    * @param gdot Slip rates on all slip systems, used to compute the evolution inputs
    * via `kinetics->getEvolVals`.
    * @param hvals Additional per-slip-system values the model's internal solve needs
    * (forwarded to `getSdotN`).
    * @param tkelv Temperature [Kelvin].
    * @param outputLevel Verbosity passed through to the SNLS solver's diagnostic output.
    * @return The number of function evaluations used by the (possibly relaxed) solve, or
    * **-1 if the solve failed to converge** even after the relaxed retry (if enabled) --
    * signaled this way rather than by throwing, since exceptions are unavailable on the
    * GPU.
    */
   template<class Kinetics, bool relaxed_solver = false>
   __ecmech_hdev__
   inline
   int
   updateHN(const Kinetics* const kinetics,
            double* hs_n,
            const double* const hs_o,
            double dt,
            const double* const gdot,
            const double* const hvals,
            double tkelv,
            int outputLevel = 0)
   {
      double evolVals[Kinetics::nEvolVals];
      kinetics->getEvolVals(evolVals, gdot);

      Kinetics_HNProblem<Kinetics> prob(kinetics, hs_o, dt, evolVals, hvals, tkelv);
      snls::SNLSTrDlDenseG<Kinetics_HNProblem<Kinetics> > solver(prob);

      snls::TrDeltaControl deltaControl;
      deltaControl._deltaInit = 1e0;
      {
         int maxIter = 100;
         double tolerance = 1e-10;
         solver.setupSolver(maxIter, tolerance, &deltaControl, outputLevel);
      }

      for (int iX = 0; iX < prob.nDimSys; ++iX) {
         solver._x[iX] = 0e0;
      }

      snls::SNLSStatus_t status = solver.solve( );

      int nFevals = solver.getNFEvals();
      if (status != snls::converged) {
         snls::SNLSStatus_t status2 = status;
         if constexpr(relaxed_solver) {
            {
               int maxIter = 100;
               double tolerance = 1e-9;
               solver.setupSolver(maxIter, tolerance, &deltaControl, outputLevel);
            }
            for (int iX = 0; iX < prob.nDimSys; ++iX) {
               solver._x[iX] = 0e0;
            }
            status2= solver.solve();
            nFevals = solver.getNFEvals();
         }
         if (status2 != snls::converged) {
            nFevals = -1;
         }
      }

      prob.getHn(hs_n, solver._x);

      return nFevals;
   } // updateHN
} // namespace ecmech

// Concrete kinetics model implementations, satisfying the interface contract documented
// at the top of this file.
#include "ECMech_kinetics_KMBalD.h"
#include "ECMech_kinetics_VocePL.h"
#include "ECMech_kinetics_OrowanD.h"
#include "ECMech_kinetics_BCCMD.h"

#endif // ECMECH_KINETICS_H
