/**
 * @file ECMech_evptn.h
 * @brief Coupled elasto-viscoplastic problem formulations for crystal plasticity.
 * 
 * This file provides the top-level problem classes that combine:
 * - Elastic lattice strain evolution (EvptnLatticeStrainProblem)
 * - Lattice rotation evolution (EvptnLatticeRotationProblem)
 * - Slip system kinetics and hardening
 * 
 * into unified nonlinear systems suitable for implicit SNLS solvers.
 * 
 * **Problem formulations**:
 * 
 * 1. **EvptnUpdstProblem**: Coupled elastic strain + rotation update
 *    - Primary formulation for most crystal plasticity simulations
 *    - Solves for both elastic strain and lattice orientation
 *    - System size: ntvec + nwvec = 8 unknowns (5 strain + 3 rotation)
 *    - Includes hardening update as sub-iteration
 * 
 * 2. **EvptnNRUpdstProblem**: Elastic strain only (fixed rotation)
 *    - Simplified formulation for cases with negligible lattice rotation
 *    - Solves only for elastic strain (orientation extrapolated)
 *    - System size: ntvec = 5 unknowns
 *    - Faster convergence but less accurate for large rotations
 *    - Enabled via ECMECH_EXTRA_SOLVERS compile flag
 * 
 * 3. **RotUpdProblem**: Rotation only (fixed elastic strain)
 *    - Specialized formulation for pure rotation updates
 *    - Useful in operator-split schemes
 *    - System size: nwvec = 3 unknowns
 *    - Enabled via ECMECH_EXTRA_SOLVERS compile flag
 * 
 * **Key features**:
 * - **SNLS interface**: Provides computeRJ() for residual and Jacobian evaluation
 * - **Automatic scaling**: Normalizes unknowns for better conditioning
 * - **Material tangent**: Optional stiffness matrix computation for FEM codes
 * - **GPU compatibility**: All methods marked __ecmech_hdev__ for device execution
 * - **Template flexibility**: Works with any SlipGeom, Kinetics, ThermoElastN combination
 * 
 * **Integration workflow**:
 * ```cpp
 * // 1. Create problem instance
 * EvptnUpdstProblem prob(slipGeom, kinetics, elastN, prob_state);
 * 
 * // 2. Create SNLS solver
 * snls::SNLSTrDlDenseG solver(prob);
 * 
 * // 3. Solve nonlinear system
 * solver.solve();
 * 
 * // 4. Extract solution
 * prob.stateFromX(elast_new, quat_new, solver._x);
 * ```
 * 
 * **Numerical approach**:
 * - Implicit backward Euler time integration
 * - Newton-Raphson solution via SNLS trust region solver
 * - Analytical Jacobian for quadratic convergence
 * - Adaptive scaling based on loading rate magnitude
 * 
 * @see ECMech_evptnSngl.h for wrapper functions that handle full time integration
 * @see ECMech_base_classes.h for subproblem formulations
 * @see ECMech_base_fcns.h for kinematic utility functions
 * @see SNLS_TrDLDenseG.h for nonlinear solver interface
 */
// -*-c++-*-

#ifndef ECMECH_EVPTN_H
#define ECMECH_EVPTN_H

#include <cassert>

#include "ECMech_core.h"
#include "ECMech_util.h"
#include "evptn/ECMech_base_classes.h"
#include "evptn/ECMech_base_fcns.h"

namespace ecmech {
   namespace evptn {

      /**
       * @brief Coupled elastic strain and lattice rotation update problem.
       * 
       * This class formulates the primary nonlinear system for crystal plasticity
       * time integration, solving simultaneously for:
       * - Elastic strain increment in crystal frame (5 DOFs)
       * - Lattice rotation increment (3 DOFs)
       * - Total system: 8 unknowns
       * 
       * **Governing equations**:
       * ```
       * R_ε: Δε_e = Δt * (D_sample - D_plastic(ε_e, ω))
       * R_ω: Δω = Δt * (W_sample - W_plastic(ε_e, ω) - W_elastic(ε_e, ω))
       * ```
       * 
       * **Coupling**:
       * - Elastic strain ε_e → stress τ → slip rates γ̇ → plastic rates D_p, W_p
       * - Lattice rotation ω → crystal frame → transformed kinematics → slip activity
       * - Strong two-way coupling requires simultaneous solution
       * 
       * **Scaling strategy**:
       * - Strain unknowns scaled by e_scale ≈ 5e-4 (typical elastic strain magnitude)
       * - Rotation unknowns scaled by r_scale ≈ 0.01 (typical rotation increment)
       * - Residuals scaled by inverse reference rates for O(1) conditioning
       * - Adaptive scaling based on ‖D_sample‖ to handle variable loading rates
       * 
       * **Template parameters**:
       * @tparam SlipGeom Slip geometry class (e.g., SlipGeomFCC)
       * @tparam Kinetics Kinetics class (e.g., KineticsKMBalD)
       * @tparam ThermoElastN Thermoelastic model (e.g., ThermoElastNCubic)
       * @tparam ProblemState State container (typically ProblemState<...>)
       * 
       * **Interface requirements** (SNLS compatibility):
       * - `static constexpr int nDimSys`: System dimension (8)
       * - `bool computeRJ(resid, Jacobian, x)`: Evaluate residual and Jacobian
       * - `void stateFromX(...)`: Extract physical state from solution vector
       * 
       * **Material tangent option**:
       * - If m_mtan_sI != nullptr, accumulate stiffness contributions
       * - Used by FEM codes for consistent tangent operator
       * - Computed via finite difference of Jacobian (implicit function theorem)
       * 
       * @see EvptnUpdstProblem::computeRJ() for the core evaluation method
       * @see getResponseSngl() for high-level usage
       */
      template<class SlipGeom, class Kinetics, class ThermoElastN, class ProblemState>
      class EvptnUpdstProblem
      {
         public:
           /**
            * @brief Total system dimension: elastic strain + rotation.
            * Value: ntvec + nwvec = 5 + 3 = 8 unknowns
            */
            static constexpr int nDimSys = ecmech::ntvec + ecmech::nwvec;

           /**
            * @brief Constructor: Initialize coupled update problem.
            * 
            * Sets up the nonlinear system for one time step by:
            * 1. Storing references to constitutive models
            * 2. Initializing elastic strain and rotation subproblems
            * 3. Extracting kinetic values at current state
            * 4. Computing adaptive scaling factors
            * 5. Storing prescribed kinematics
            * 
            * **Scaling computation**:
            * ```
            * reference_rate = max(‖D_sample‖, fixed_ref_rate)
            * epsdot_scale_inv = 1 / reference_rate
            * rotincr_scale_inv = (1/Δt) * epsdot_scale_inv
            * ```
            * This makes residuals O(1) for better Newton convergence.
            * 
            * @param[in] slipGeom Slip geometry defining slip systems and Schmid tensors
            * @param[in] kinetics Kinetics model for rate-dependent strength
            * @param[in] thermoElastN Thermoelastic constitutive law
            * @param[in] prob_state Reference to material point state container
            *                       Contains beginning-of-step values and storage for updates
            * 
            * **State dependencies**:
            * - prob_state.dt: Time step size
            * - prob_state.rel_vol_new: Volume for finite strain corrections
            * - prob_state.energy_new: Internal energy for thermoelasticity
            * - prob_state.pressure_EOS: Pressure from equation of state
            * - prob_state.tkelv: Temperature for kinetics
            * - prob_state.elast_d5_n: Beginning-of-step elastic strain
            * - prob_state.quat_n: Beginning-of-step orientation
            * - prob_state.def_rate_d5_sample: Prescribed deformation rate
            * - prob_state.spin_vec_sample: Prescribed spin
            * - prob_state.h_state_u: Hardening state (updated in constructor)
            * 
            * **Kinetic value extraction**:
            * - m_kinetics.getVals() populates m_kin_vals array
            * - Returns hardening scale factor (m_hdn_scale)
            * - Kinetic values include slip resistances, mobilities, etc.
            * 
            * @note Constructor computes derived quantities but does NOT solve the system
            * @note Must be followed by solver.solve() call
            */
            __ecmech_hdev__
            EvptnUpdstProblem(const SlipGeom& slipGeom,
                              const Kinetics& kinetics,
                              const ThermoElastN& thermoElastN,
                              ProblemState& prob_state
                              )
               : m_slipGeom(slipGeom),
               m_kinetics(kinetics),
               m_lattice_strain_prob(thermoElastN, prob_state.dt, prob_state.rel_vol_new, prob_state.energy_new, prob_state.pressure_EOS, prob_state.tkelv, prob_state.elast_d5_n),
               m_lattice_rot_prob(prob_state.dt, prob_state.quat_n),
               m_def_rate_d5_sample(prob_state.def_rate_d5_sample), // vel_grad_sm%d_vecds
               m_spin_vec_sample(prob_state.spin_vec_sample), // vel_grad_sm%w_veccp
               m_mtan_sI(nullptr)
            {
               m_hdn_scale = m_kinetics.getVals(m_kin_vals, prob_state.pressure_EOS, prob_state.tkelv, prob_state.h_state_u);

               double adots_ref = m_kinetics.getFixedRefRate(m_kin_vals);
               double eff = vecNorm<ntvec>(m_def_rate_d5_sample); // do not worry about factor of sqrt(twothird)
               if (eff < epsdot_scl_nzeff * adots_ref) {
                  m_epsdot_scale_inv = one / adots_ref;
               }
               else {
                  m_epsdot_scale_inv = fmin(one / eff, 1e6 * m_lattice_strain_prob.m_dt);
               }
               //
               m_rotincr_scale_inv = m_lattice_strain_prob.m_inv_dt * m_epsdot_scale_inv;
            }

            /**
             * @brief Destructor (default).
             */
            __ecmech_hdev__
            ~EvptnUpdstProblem() {}

           /**
            * @brief Enable material tangent stiffness computation.
            * 
            * Sets pointer to array where tangent stiffness should be accumulated.
            * When non-null, computeRJ() will compute and store ∂σ/∂ε contributions.
            * 
            * @param[in,out] mtan_sI Pointer to 6×6 tangent stiffness matrix (Voigt form)
            *                        Modified during computeRJ() calls
            * 
            * @see computeTangentStiffness() in ECMech_evptnSngl.h for post-processing
            */
            __ecmech_hdev__
            inline
            void provideMTan(double* mtan_sI) { m_mtan_sI = mtan_sI; }

           /**
            * @brief Disable material tangent stiffness computation.
            * 
            * Resets tangent pointer to nullptr, disabling stiffness accumulation.
            * Should be called after tangent computation is complete.
            */
            __ecmech_hdev__
            inline
            void clearMTan( ) { m_mtan_sI = nullptr; }

           /**
            * @brief Get inverse time step (1/Δt).
            * 
            * Used for converting strain increments to rates in post-processing.
            * 
            * @return Inverse time step [1/time]
            */
            __ecmech_hdev__
            inline
            double getDtRi() const { return m_lattice_strain_prob.m_inv_dt; }

           /**
            * @brief Get hardening scale factor.
            * 
            * Returns the current flow strength scale from kinetics.
            * Used in flow stress calculations.
            * 
            * @return Hardening scale factor [stress units]
            */
            __ecmech_hdev__
            inline
            double getHdnScale() const { return m_hdn_scale; }

           /**
            * @brief Reconstruct physical state from solution vector.
            * 
            * Decodes the scaled solution vector x into:
            * - End-of-step elastic strain (deviatoric 5-vector)
            * - End-of-step orientation quaternion (unit quaternion)
            * 
            * **Inverse operations**:
            * ```
            * elast_d5 = elast_d5_n + e_scale * x[0:4]
            * quat = exp(r_scale * x[5:7] / 2) ∘ quat_n
            * ```
            * 
            * @param[out] elast_dev_press_vec Elastic strain ε_e [ntvec = 5 components]
            * @param[out] quat Orientation quaternion Q [qdim = 4 components]
            * @param[in] x Solution vector [nDimSys = 8 components]
            * 
            * **Safety notes**:
            * - Not safe if output arrays alias input or internal storage
            * - Quaternion normalized during extraction for numerical stability
            * 
            * @see EvptnLatticeStrainProblem::stateFromX() for elastic strain extraction
            * @see EvptnLatticeRotationProblem::stateFromX() for quaternion extraction
            */
            __ecmech_hdev__
            inline
            void stateFromX(double* const elast_dev_press_vec,
                            double* const quat,
                            const double* const x) {
               m_lattice_strain_prob.stateFromX(elast_dev_press_vec,  &(x[m_i_sub_e]));
               m_lattice_rot_prob.stateFromX(quat,  &(x[m_i_sub_r]));
            }

           /**
            * @brief Convert elastic strain to Cauchy stress.
            * 
            * Convenience method for stress evaluation without needing to access
            * lattice strain subproblem directly.
            * 
            * @param[out] cauchy_xtal Cauchy stress in crystal frame [nsvec = 6 components]
            * @param[in] elast_d5_f Elastic strain [ntvec = 5 components]
            */
            __ecmech_hdev__
            inline
            void elastNEtoC(double* const cauchy_xtal, // nsvec
                            const double* const elast_d5_f // ntvec
                            ) const {
               m_lattice_strain_prob.elast_strain_to_cauchy_stress(cauchy_xtal, elast_d5_f);
            }

           /**
            * @brief Evaluate residual vector and Jacobian matrix.
            * 
            * This is the core method called by SNLS solver at each iteration.
            * Computes the nonlinear residual R(x) and optionally its Jacobian J(x).
            * 
            * **Residual structure**:
            * ```
            * R[0:4]  : Elastic strain residual R_ε
            * R[5:7]  : Lattice rotation residual R_ω
            * ```
            * 
            * **Jacobian structure** (8×8 matrix):
            * ```
            * J = | ∂R_ε/∂ε_e   ∂R_ε/∂ω |  [5×5]  [5×3]
            *     | ∂R_ω/∂ε_e   ∂R_ω/∂ω |  [3×5]  [3×3]
            * ```
            * 
            * **Computation sequence**:
            * 1. Decode solution vector x → elastic strain, rotation
            * 2. Transform kinematics to crystal frame
            * 3. Compute stress from elastic strain
            * 4. Resolve stress onto slip systems
            * 5. Evaluate slip rates and derivatives
            * 6. Compute plastic deformation and spin
            * 7. Evaluate elastic strain residual
            * 8. Evaluate rotation residual (with higher-order corrections)
            * 9. If Jacobian requested:
            *    a. Compute slip rate derivatives w.r.t. strain
            *    b. Assemble elastic-elastic block
            *    c. Assemble elastic-rotation coupling
            *    d. Assemble rotation-elastic coupling
            *    e. Assemble rotation-rotation block
            *    f. Apply all scaling transformations
            * 
            * **Scaling in Jacobian**:
            * Each block scaled by appropriate combination of:
            * - e_scale (strain unknown scaling)
            * - r_scale (rotation unknown scaling)
            * - epsdot_scale_inv (residual conditioning)
            * - rotincr_scale_inv (rotation residual conditioning)
            * 
            * **Material tangent**:
            * If m_mtan_sI != nullptr, accumulates contributions needed for
            * ∂σ/∂ε tangent stiffness (used by FEM codes).
            * 
            * @param[out] resid Residual vector R(x) [nDimSys = 8 components]
            * @param[out] Jacobian Jacobian matrix J(x) [nDimSys×nDimSys = 8×8]
            *                      Set to nullptr to skip Jacobian computation
            * @param[in] x Solution vector (scaled unknowns) [nDimSys = 8 components]
            * 
            * @return true if evaluation successful, false if computation failed
            * 
            * **Return value**:
            * - true: Residual/Jacobian computed successfully
            * - false: Numerical issue (overflow, invalid state, etc.)
            * - Currently always returns true; future versions may add failure detection
            * 
            * @note Jacobian zeroed at start if requested (no need to pre-initialize)
            * @note Residual always evaluated; Jacobian only if pointer non-null
            * @note This method is called many times per time step (Newton iterations)
            * 
            * @see SNLS solver documentation for interface requirements
            */
            __ecmech_hdev__
            inline
            bool computeRJ(double* const resid,
                           double* const Jacobian,
                           const double* const x) {
               bool doComputeJ = (Jacobian != nullptr);

               if (doComputeJ) {
                  // zero the Jacobian so that do not need to worry about zero
                  // entries in the midst of other things later
                  //
                  for (size_t ijJ = 0; ijJ<m_nXnDim; ++ijJ) {
                     Jacobian[ijJ] = 0.0;
                  }
               }
               //
               for (int iR = 0; iR<nDimSys; ++iR) {
                  resid[iR] = 0.0;
               }

               //////////////////////////////
               // PULL VALUES out of x, with scalings
               //
               double elast_dt_d5[ecmech::ntvec];
               vecsVxa<ntvec>(elast_dt_d5, ecmech::e_scale, &(x[m_i_sub_e]) ); // elast_dt_d5 is now the delta, _not_ yet elast_dt_d5
               // elast_d5_f is end-of-step
               double elast_d5_f[ntvec];
               vecsVapb<ntvec>(elast_d5_f, elast_dt_d5, m_lattice_strain_prob.m_elast_d5_n);
               vecsVsa<ntvec>(elast_dt_d5, m_lattice_strain_prob.m_inv_dt); // _now_ elast_dt_d5 has elast_dt_d5
               //
               double xi_f[nwvec];
               vecsVxa<nwvec>(xi_f, ecmech::r_scale, &(x[m_i_sub_r]) );

               double A_quat[ecmech::qdim];
               emap_to_quat(A_quat, xi_f);
               //
               double xtal_ori_quat[ecmech::qdim];
               get_c_quat(xtal_ori_quat, A_quat, m_lattice_rot_prob.m_xtal_ori_quat_n);
               //
               double xtal_rmat[ecmech::ndim * ecmech::ndim];
               quat_to_tensor(xtal_rmat, xtal_ori_quat);
               //
               double rmat_5x5_sample2xtal[ecmech::ntvec * ecmech::ntvec];
               get_rot_mat_vecd(rmat_5x5_sample2xtal, xtal_rmat);

               double def_rate_d5_xtal[ecmech::ntvec];
               double spin_vec_xtal[ecmech::nwvec]; // assumes nwvec = ndim

               get_xtal_frame_vel_grad_terms(def_rate_d5_xtal, spin_vec_xtal, m_def_rate_d5_sample,  m_spin_vec_sample, xtal_rmat, rmat_5x5_sample2xtal);

               //////////////////////////////
               // CALCULATIONS

               double kirchoff[ecmech::nsvec];
               m_lattice_strain_prob.elast_strain_to_kirchoff_stress(kirchoff, elast_d5_f);

               double dgdot_dtau[SlipGeom::nslip] = { 0.0 }; // crys%tmp2_slp
               double plastic_def_rate_d5[ecmech::ntvec] = { 0.0 };
               double plastic_spin_vec[ecmech::nwvec] = { 0.0 }; // \pcDhat

               get_slip_rate_terms(dgdot_dtau, plastic_def_rate_d5, plastic_spin_vec, kirchoff, m_kin_vals, m_slipGeom, m_kinetics);

               // Higher-order terms related to the elasticity stuff that's used in the residuals and jacobian calculation
               double A_e_M35[ecmech::nwvec * ecmech::ntvec];
               double ee_wvec[ecmech::nwvec];
               double ee_fac;
               elasticity_higher_order_terms(A_e_M35, ee_wvec, ee_fac, m_lattice_strain_prob.m_inv_a_vol, elast_d5_f, elast_dt_d5);

               // Residual Calculations
               m_lattice_strain_prob.get_elast_strain_residual(resid, m_epsdot_scale_inv, elast_dt_d5, plastic_def_rate_d5, def_rate_d5_xtal);
               m_lattice_rot_prob.get_omega_residual(resid, m_rotincr_scale_inv, ee_fac, xi_f, spin_vec_xtal, plastic_spin_vec, ee_wvec);

               //////////////////////////////////////////////////////////////////////
               // JACOBIAN, fixed hardness and temperature
               //
               if (doComputeJ) {
                  // use RAJA::View machinery to simplify indexing for blocks in the Jacobian matrix ;
                  // can always swap this out later if it ends up being too heavyweight ;
                  // RAJA defaults to "row-major" -- final dimension indexing the fastest
                  //
                  // preliminaries
                  //
                  double dpl_deps_symm[ ecmech::ntvec * ecmech::ntvec ] = { 0.0 };
                  double dpl_deps_skew[ ecmech::nwvec * ecmech::ntvec ] = { 0.0 };
                  get_slip_rate_deriv_terms(dpl_deps_symm, dpl_deps_skew, dgdot_dtau, m_lattice_strain_prob.m_inv_a_vol, m_slipGeom, m_lattice_strain_prob.m_thermo_elast_n);

                  // derivatives with respect to lattice orientation changes
                  double dxtal_ori_quat_dxi_T[ ecmech::nwvec * ecmech::qdim ];
                  double dDsm_dxi[ ecmech::ntvec * ecmech::nwvec ];
                  double dWsm_dxi[ ecmech::nwvec * ecmech::nwvec ];
                  eval_d_dxi_impl_quat(dxtal_ori_quat_dxi_T, dDsm_dxi, dWsm_dxi,
                                       m_def_rate_d5_sample, m_spin_vec_sample,
                                       xi_f, 
                                       m_lattice_rot_prob.m_xtal_ori_quat_n,
                                       xtal_rmat, xtal_ori_quat);

                  // d(B_S)/d(elast_d5_f)
                  //
                  m_lattice_strain_prob.template get_deriv_elast_strain_wrt_elast_strain<nDimSys>(Jacobian, dpl_deps_symm);
                  // d(B_S)/d(xi_f)
                  //
                  // jacob_er = -dDsm_dxi(:,:)
                  m_lattice_rot_prob.template get_deriv_elast_strain_wrt_omega<nDimSys>(Jacobian, dDsm_dxi);
                  // d(B_xi)/d(elast_dev_press_vecs_f)
                  //
                  m_lattice_strain_prob.template get_deriv_omega_wrt_elast_strain<nDimSys, m_i_sub_r>(Jacobian, elast_dt_d5, ee_fac, dpl_deps_skew, A_e_M35);
                  // d(B_xi)/d(xi_f)
                  //
                  m_lattice_rot_prob.template get_deriv_omega_wrt_omega<nDimSys>(Jacobian, dWsm_dxi);

                  if (m_mtan_sI) {
                     double cauchy_stress_lattice[ ecmech::nsvec ];
                     m_lattice_strain_prob.m_thermo_elast_n.getCauchy(cauchy_stress_lattice, kirchoff, m_lattice_strain_prob.m_inv_det_v_e);
                     get_material_tangent_stiffness<ThermoElastN, nDimSys, m_i_sub_r>
                     (m_mtan_sI, Jacobian,
                     dxtal_ori_quat_dxi_T, rmat_5x5_sample2xtal,
                     xtal_ori_quat, xtal_rmat,
                     cauchy_stress_lattice,
                     m_lattice_strain_prob.m_inv_det_v_e,
                     m_lattice_strain_prob.m_inv_a_vol,
                     m_lattice_strain_prob.m_thermo_elast_n);
                  }

                  // SCALING
                  {
                     double scaleFactorJ;
                     for (size_t iJ = 0; iJ<m_i_sub_r; ++iJ) {
                        // Jacobian(i_sub_e:i_sup_e,i_sub_e:i_sup_e) = jacob_ee * epsdot_scale_inv  * e_scale ! resid, x
                        scaleFactorJ = m_epsdot_scale_inv * ecmech::e_scale;
                        for (size_t jJ = 0; jJ<m_i_sub_r; ++jJ) { // <=_i_sup_e
                           int ijJ = ECMECH_NN_INDX(iJ, jJ, nDimSys);
                           Jacobian[ ijJ ] *= scaleFactorJ;
                        }

                        // Jacobian(i_sub_e:i_sup_e,i_sub_r:i_sup_r) = jacob_er * epsdot_scale_inv  * r_scale
                        scaleFactorJ = m_epsdot_scale_inv * ecmech::r_scale;
                        for (int jJ = m_i_sub_r; jJ<nDimSys; ++jJ) { // <_i_sup_r
                           int ijJ = ECMECH_NN_INDX(iJ, jJ, nDimSys);
                           Jacobian[ ijJ ] *= scaleFactorJ;
                        }
                     }

                     for (int iJ = m_i_sub_r; iJ<nDimSys; ++iJ) {
                        // Jacobian(i_sub_r:i_sup_r,i_sub_e:i_sup_e) = jacob_re * rotincr_scale_inv * e_scale
                        scaleFactorJ = m_rotincr_scale_inv * ecmech::e_scale;
                        for (size_t jJ = 0; jJ<m_i_sub_r; ++jJ) { // <=_i_sup_e
                           int ijJ = ECMECH_NN_INDX(iJ, jJ, nDimSys);
                           Jacobian[ ijJ ] *= scaleFactorJ;
                        }

                        // Jacobian(i_sub_r:i_sup_r,i_sub_r:i_sup_r) = jacob_rr * rotincr_scale_inv * r_scale
                        scaleFactorJ = m_rotincr_scale_inv * ecmech::r_scale;
                        for (size_t jJ = m_i_sub_r; jJ<nDimSys; ++jJ) { // <_i_sup_r
                           int ijJ = ECMECH_NN_INDX(iJ, jJ, nDimSys);
                           Jacobian[ ijJ ] *= scaleFactorJ;
                        }
                     }
                  } // SCALING
               }
               return true;
            } // computeRJ

           /**
            * @brief Calculate slip-related quantities from current elastic strain.
            * 
            * Post-processing method that computes:
            * - Plastic dissipation rate
            * - Effective shear rate
            * - Individual slip rates
            * 
            * Used for output and diagnostics after convergence.
            * 
            * @param[out] pl_disipation_rate Plastic power [stress/time]
            * @param[out] effective_shear_rate Effective γ̇ [1/time]
            * @param[out] gdot Slip rates [1/time, nslip components]
            * @param[in] elast_strain Elastic strain [ntvec = 5 components]
            * 
            * @see get_slip_contributions() for implementation details
            */
            __ecmech_hdev__
            inline
            void get_slip_contribution(double& pl_disipation_rate,
                                       double& effective_shear_rate,
                                       double* const gdot,
                                       const double* const elast_strain
                                      )
            {
               get_slip_contributions(pl_disipation_rate, effective_shear_rate, gdot,
                                      m_lattice_strain_prob.m_inv_det_v_e, elast_strain, m_kin_vals,
                                      m_slipGeom, m_kinetics, m_lattice_strain_prob);
            }
                              

         private:
            /** @brief Reference to slip geometry object */
            const SlipGeom &m_slipGeom;
            /** @brief Reference to kinetics object */
            const Kinetics &m_kinetics;
            /** @brief Elastic strain subproblem instance */
            const EvptnLatticeStrainProblem<ThermoElastN> m_lattice_strain_prob;
            /** @brief Lattice rotation subproblem instance */
            const EvptnLatticeRotationProblem<ecmech::ntvec> m_lattice_rot_prob;
            /** @brief Hardening scale factor [stress units] */
            double m_hdn_scale;
            /** @brief Inverse strain rate scaling [time] */
            double m_epsdot_scale_inv;
            /** @brief Inverse rotation increment scaling [time] */
            double m_rotincr_scale_inv;
            /** @brief Kinetic state values array [Kinetics::nVals components] */
            double m_kin_vals[Kinetics::nVals];
            /** @brief Prescribed deformation rate in sample frame [ntvec] */
            const double* const m_def_rate_d5_sample;
            /** @brief Prescribed spin in sample frame [nwvec] */
            const double* const m_spin_vec_sample;
            /** @brief Jacobian size (nDimSys squared) */
            static constexpr size_t m_nXnDim = nDimSys * nDimSys;
            /** @brief Starting index for elastic strain unknowns */
            static constexpr size_t m_i_sub_e = 0;
            /** @brief Starting index for rotation unknowns */
            static constexpr size_t m_i_sub_r = ecmech::ntvec;
            /** @brief Pointer to material tangent stiffness (nullptr if not computing) */
            double* m_mtan_sI;
      }; // class EvptnUpdstProblem

#if defined(ECMECH_EXTRA_SOLVERS)

      /**
       * @brief Lattice rotation-only update problem for crystal plasticity.
       * 
       * This class formulates a simplified integration problem where only the lattice
       * rotation (crystal orientation) is updated, while elastic strain is held fixed.
       * This decoupled approach can converge when the fully coupled EvptnUpdstProblem
       * fails, trading some accuracy for robustness.
       * 
       * **Physical interpretation**:
       * Updates crystal orientation Q due to:
       * - Prescribed spin in sample frame (rotated to crystal frame)
       * - Plastic spin from slip (computed from fixed elastic strain)
       * - Elastic spin corrections (geometric terms)
       * 
       * **When this is used**:
       * - Backup solver in getResponseRetry() after primary solver failure
       * - Optional R* (rotation rate) solve in preprocessing
       * - Debugging/testing decoupled rotation integration
       * 
       * **Simplifications vs. EvptnUpdstProblem**:
       * - Updates: Lattice rotation (3 DOFs)
       * - Fixed: Elastic strain (uses beginning-of-step values)
       * - Fixed: Hardening state (uses beginning-of-step values)
       * - Result: Simpler 3×3 system instead of 8×8 coupled system
       * 
       * **System dimensions**:
       * - nDimSys = nwvec = 3 (rotation DOFs only)
       * - Unknowns: ω = [ω₁₂, ω₂₃, ω₁₃] (incremental rotation vector)
       * 
       * **Residual equation**:
       * ```
       * R_ω = ω - Δt * (W_lat - W_p + ee_fac * W_e)
       * ```
       * where all spins are evaluated with fixed elastic strain.
       * 
       * **Template parameters**:
       * @tparam SlipGeom Slip geometry class (e.g., SlipGeomFCC)
       * @tparam ThermoElastN Thermoelastic model (e.g., ThermoElastNCubic)
       * @tparam ProblemState State container (ProblemState<...>)
       * 
       * **Compile-time flags**:
       * - Requires: ECMECH_EXTRA_SOLVERS defined
       * - Used with: ECMECH_DEBUG for additional checks
       * 
       * @see EvptnUpdstProblem for fully coupled integration
       * @see getResponseRetry() for usage in convergence recovery
       * @see preprocess() for optional R* solve usage
       */
      template<class SlipGeom, class ThermoElastN, class ProblemState>
      class RotUpdProblem
      {
         public:
        /**
         * @brief System dimension (number of DOFs).
         * 
         * Value: nwvec = 3 (rotation vector components)
         */
         static constexpr int nDimSys = ecmech::nwvec;


        /**
         * @brief Constructor initializes rotation-only update problem.
         * 
         * Sets up the rotation integration problem using fixed elastic strain
         * and slip system state from the beginning of the time step.
         * 
         * **Initialization sequence**:
         * 1. Store slip geometry and elastic state
         * 2. Create lattice strain and rotation subproblems
         * 3. Compute adaptive scaling factors for conditioning
         * 4. Evaluate plastic deformation and spin from beginning-of-step slip rates
         * 
         * **Scaling strategy**:
         * The solver uses adaptive scaling to improve conditioning:
         * ```
         * adots_ref = ‖gdot‖  (reference slip rate magnitude)
         * eff = ‖D_sample‖     (deformation rate magnitude)
         * 
         * if eff < epsdot_scl_nzeff * adots_ref:
         *     epsdot_scale_inv = 1 / adots_ref
         * else:
         *     epsdot_scale_inv = min(1 / eff, 1e6 * dt)
         * 
         * rotincr_scale_inv = (1/dt) * epsdot_scale_inv
         * ```
         * 
         * This ensures residuals and Jacobian entries are O(1) for better Newton convergence.
         * 
         * **Plastic quantities from fixed state**:
         * Using beginning-of-step slip rates gdot (from prob_state):
         * ```
         * D_p = Σ_α γ̇^α P^α  (plastic deformation rate)
         * W_p = Σ_α γ̇^α Q^α  (plastic spin)
         * ```
         * where P^α and Q^α are slip system Schmid tensors.
         * 
         * @param[in] slipGeom Slip system geometry (Schmid tensors, etc.)
         *                     Provides P (symmetric) and Q (skew) for each slip system
         * @param[in] thermoElastN Thermoelastic model for stress evaluation
         *                          Needed to construct lattice strain subproblem
         * @param[in] prob_state Problem state container with:
         *                       - dt: Time step size
         *                       - gdot: Beginning-of-step slip rates
         *                       - elast_d5_n: Beginning-of-step elastic strain
         *                       - quat_n: Beginning-of-step orientation
         *                       - def_rate_d5_sample: Prescribed deformation rate
         *                       - spin_vec_sample: Prescribed spin
         *                       - rel_vol_new, energy_new, pressure_EOS, tkelv: Thermodynamic state
         * 
         * **Member initialization**:
         * - m_lattice_strain_prob: Elastic state (fixed, for stress evaluation)
         * - m_lattice_rot_prob: Rotation evolution (active unknowns)
         * - m_plastic_def_rate_d5: D_p from beginning-of-step slip
         * - m_plastic_spin_vec: W_p from beginning-of-step slip
         * - m_epsdot_scale_inv, m_rotincr_scale_inv: Adaptive scaling factors
         * 
         * @note Elastic strain and hardening remain at beginning-of-step values
         * @note Plastic quantities computed once during construction (not updated)
         * @note Scaling factors ensure O(1) residuals for solver
         * 
         * @see EvptnLatticeStrainProblem for elastic strain handling
         * @see EvptnLatticeRotationProblem for rotation evolution
         */
         __ecmech_hdev__
         RotUpdProblem(const SlipGeom& slipGeom,
                        const ThermoElastN& thermoElastN,
                        ProblemState& prob_state
                        ) :
            m_slipGeom(slipGeom),
            m_lattice_strain_prob(thermoElastN, prob_state.dt, prob_state.rel_vol_new, prob_state.energy_new, prob_state.pressure_EOS, prob_state.tkelv, prob_state.elast_d5_n),
            m_lattice_rot_prob(prob_state.dt, prob_state.quat_n),
            m_def_rate_d5_sample(prob_state.def_rate_d5_sample), // vel_grad_sm%d_vecds
            m_spin_vec_sample(prob_state.spin_vec_sample) // vel_grad_sm%w_veccp
         {

            double adots_ref = vecNorm<SlipGeom::nslip>(prob_state.gdot);

            double eff = vecNorm<ecmech::ntvec>(m_def_rate_d5_sample); // do not worry about factor of sqrt(twothird)
            if (eff < epsdot_scl_nzeff * adots_ref) {
                  m_epsdot_scale_inv = one / adots_ref;
            }
            else {
                  m_epsdot_scale_inv = fmin(one / eff, 1e6 * m_lattice_strain_prob.m_dt);
            }
            //
            m_rotincr_scale_inv = m_lattice_strain_prob.m_inv_dt * m_epsdot_scale_inv;

            vecsVMa<ntvec, SlipGeom::nslip>(m_plastic_def_rate_d5, slipGeom.getP(), prob_state.gdot);
            vecsVMa<nwvec, SlipGeom::nslip>(m_plastic_spin_vec, slipGeom.getQ(), prob_state.gdot);

         }

        /**
         * @brief Destructor (default).
         */
         __ecmech_hdev__
         ~RotUpdProblem() {}

        /**
         * @brief Extract converged orientation quaternion from solution vector.
         * 
         * Converts the scaled solution vector (rotation increment in solver space)
         * back to physical quaternion representation.
         * 
         * **Conversion sequence**:
         * 1. Unscale: ω = r_scale * x
         * 2. Convert to quaternion increment: A = exp(ω/2)
         * 3. Compose with beginning-of-step: Q^{n+1} = A ∘ Q^n
         * 
         * **Mathematical details**:
         * - x: Scaled rotation vector (solver unknowns)
         * - ω: Physical rotation vector [radians]
         * - A: Quaternion representing incremental rotation
         * - Q^{n+1}: End-of-step orientation
         * 
         * @param[out] quat End-of-step orientation quaternion [qdim = 4]
         *                  Represents sample→crystal frame rotation
         *                  Normalized during conversion
         * @param[in] x Scaled solution vector [nDimSys = 3]
         *              Rotation increment in solver space
         * 
         * @note quat must be pre-allocated (4 doubles)
         * @note Output quaternion is automatically normalized
         * @note Safe to use same memory for quat as prob_state.quat_u
         * 
         * @see EvptnLatticeRotationProblem::stateFromX for implementation
         * @see emap_to_quat() for exponential map conversion
         * @see get_c_quat() for quaternion composition
         */
         __ecmech_hdev__
         inline
         void stateFromX(double* const quat,
                        const double* const x) {
            m_lattice_rot_prob.stateFromX(quat, x);
         }


        /**
         * @brief Compute residual and Jacobian for SNLS trust region solver.
         * 
         * Evaluates the rotation evolution equation and its derivatives with
         * respect to rotation unknowns. This is the core computation for the
         * implicit integration of crystal orientation.
         * 
         * **Residual formulation**:
         * ```
         * R_ω = (ω - Δt * (W_lat - W_p + ee_fac * W_e)) * rotincr_scale_inv
         * ```
         * At convergence: R_ω = 0, implying ω = Δt * (W_lat - W_p + ee_fac * W_e)
         * 
         * **Computation sequence**:
         * 
         * 1. **Extract rotation state** from solution vector x:
         *    - Unscale: ω = r_scale * x
         *    - Convert to quaternion: Q^{n+1} = exp(ω) ∘ Q^n
         *    - Generate rotation matrices: R (3×3), R₅ (5×5)
         * 
         * 2. **Transform velocity gradient** to crystal frame:
         *    - D_xtal = R₅^T · D_sample · R₅
         *    - W_lat = R^T · W_sample · R
         * 
         * 3. **Compute elastic strain rate** (from fixed elastic strain):
         *    - ε̇_e = (1/a_vol) * (D_xtal - D_p)
         *    - D_p computed in constructor from beginning-of-step slip rates
         * 
         * 4. **Higher-order elastic terms**:
         *    - A_e, ee_spin, ee_fac = elasticity_higher_order_terms(...)
         *    - Geometric corrections for finite strain
         * 
         * 5. **Evaluate rotation residual**:
         *    - R_ω = lattice_rot_prob.get_omega_residual(...)
         *    - Uses W_lat, W_p, ee_spin, ee_fac
         * 
         * 6. **Jacobian** (if requested):
         *    - Compute ∂W_lat/∂ω via eval_d_dxi_impl_quat
         *    - J = ∂R_ω/∂ω = I - Δt * ∂W_lat/∂ω
         *    - Apply scaling: rotincr_scale_inv * r_scale
         * 
         * **Jacobian structure** (3×3):
         * ```
         * [∂R₁/∂ω₁  ∂R₁/∂ω₂  ∂R₁/∂ω₃]
         * [∂R₂/∂ω₁  ∂R₂/∂ω₂  ∂R₂/∂ω₃]
         * [∂R₃/∂ω₁  ∂R₃/∂ω₂  ∂R₃/∂ω₃]
         * ```
         * 
         * **Scaling for conditioning**:
         * - Residual scaled by rotincr_scale_inv ≈ 1/‖W_sample‖
         * - Solution scaled by r_scale ≈ 0.01
         * - Jacobian scaled by product: rotincr_scale_inv * r_scale
         * - Result: O(1) entries for Newton solver
         * 
         * @param[out] resid Residual vector [nDimSys = 3]
         *                   Rotation evolution equation
         *                   Set to zero, then populated
         * @param[out] Jacobian Jacobian matrix [nDimSys × nDimSys = 9]
         *                      Row-major storage: J[i*nDimSys + j] = ∂R_i/∂x_j
         *                      Set to zero, then populated
         *                      If nullptr, Jacobian computation skipped
         * @param[in] x Solution vector [nDimSys = 3]
         *              Scaled rotation increment
         *              Solver iterates on this to find ω satisfying R(ω) = 0
         * 
         * @return true Always returns true (rotation update typically converges)
         * 
         * **Performance notes**:
         * - Jacobian computation is optional (nullptr check)
         * - Most expensive: eval_d_dxi_impl_quat (rotation derivatives)
         * - Vectorized operations where possible
         * - No dynamic allocation
         * 
         * @note resid and Jacobian zeroed at start of function
         * @note Jacobian only computed if Jacobian != nullptr
         * @note All operations use fixed elastic strain and hardening state
         * @note Thread-safe (no shared state modified)
         * 
         * @see EvptnLatticeRotationProblem::get_omega_residual for R_ω computation
         * @see eval_d_dxi_impl_quat() for rotation derivative computation
         * @see elasticity_higher_order_terms() for elastic spin corrections
         */
         __ecmech_hdev__
         bool computeRJ(double* const resid,
                        double* const Jacobian,
                        const double* const x) {
            bool doComputeJ = (Jacobian != nullptr);

            if (doComputeJ) {
                  // zero the Jacobian so that do not need to worry about zero
                  // entries in the midst of other things later
                  //
                  for (int ijJ = 0; ijJ< m_nXnDim; ++ijJ) {
                     Jacobian[ijJ] = 0.0;
                  }
            }
            //
            for (int iR = 0; iR<nDimSys; ++iR) {
                  resid[iR] = 0.0;
            }

            double xi_f[nwvec];
            vecsVxa<nwvec>(xi_f, ecmech::r_scale, x);

            double A_quat[ecmech::qdim];
            emap_to_quat(A_quat, xi_f);
            //
            double xtal_ori_quat[ecmech::qdim];
            get_c_quat(xtal_ori_quat, A_quat, m_lattice_rot_prob.m_xtal_ori_quat_n);
            //
            double xtal_rmat[ecmech::ndim * ecmech::ndim];
            quat_to_tensor(xtal_rmat, xtal_ori_quat);
            //
            double rmat_5x5_sample2xtal[ecmech::ntvec * ecmech::ntvec];
            get_rot_mat_vecd(rmat_5x5_sample2xtal, xtal_rmat);

            double def_rate_d5_xtal[ecmech::ntvec];
            double spin_vec_xtal[ecmech::nwvec]; // assumes nwvec = ndim
            get_xtal_frame_vel_grad_terms(def_rate_d5_xtal, spin_vec_xtal, m_def_rate_d5_sample,  m_spin_vec_sample, xtal_rmat, rmat_5x5_sample2xtal);

            double elast_dt_d5[ecmech::ntvec];
            // Calculate what this elast_dt_d5 term should be given the current
            // state information.
            for (int i = 0; i < ecmech::ntvec; i++)
            {
                  elast_dt_d5[i] = m_lattice_strain_prob.m_inv_a_vol * (def_rate_d5_xtal[i] - m_plastic_def_rate_d5[i]);
            }

            // Higher-order terms related to the elasticity stuff that's used in the residuals and jacobian calculation
            double A_e_M35[ecmech::nwvec * ecmech::ntvec];
            double ee_wvec[ecmech::nwvec];
            double ee_fac;
            elasticity_higher_order_terms(A_e_M35, ee_wvec, ee_fac, m_lattice_strain_prob.m_inv_a_vol, m_lattice_strain_prob.m_elast_d5_n, elast_dt_d5);

            // Residual Calculations
            m_lattice_rot_prob.get_omega_residual(resid, m_rotincr_scale_inv, ee_fac, xi_f, spin_vec_xtal, m_plastic_spin_vec, ee_wvec);

            //////////////////////////////////////////////////////////////////////
            // JACOBIAN, fixed hardness and temperature
            //
            if (doComputeJ) {

                  //
                  //
                  // derivatives with respect to lattice orientation changes
                  double dxtal_ori_quat_dxi_T[ ecmech::nwvec * ecmech::qdim ];
                  double dDsm_dxi[ ecmech::ntvec * ecmech::nwvec ];
                  double dWsm_dxi[ ecmech::nwvec * ecmech::nwvec ];
                  eval_d_dxi_impl_quat(dxtal_ori_quat_dxi_T, dDsm_dxi, dWsm_dxi,
                                    m_def_rate_d5_sample, m_spin_vec_sample,
                                    xi_f,
                                    m_lattice_rot_prob.m_xtal_ori_quat_n,
                                    xtal_rmat, xtal_ori_quat);

                  // d(B_xi)/d(xi_f)
                  //
                  m_lattice_rot_prob.template get_deriv_omega_wrt_omega<nDimSys>(Jacobian, dWsm_dxi);

                  const double scaleFactorJ = m_rotincr_scale_inv * ecmech::r_scale;
                  for (int iJ = 0; iJ<nDimSys; ++iJ) {
                     // Jacobian(i_sub_r:i_sup_r,i_sub_r:i_sup_r) = jacob_rr * rotincr_scale_inv * r_scale
                     for (int jJ = 0; jJ<nDimSys; ++jJ) { // <_i_sup_r
                        int ijJ = ECMECH_NN_INDX(iJ, jJ, nDimSys);
                        Jacobian[ ijJ ] *= scaleFactorJ;
                     }
                  }
            }
            return true;
         }

         private:

         /**
          * @brief Slip geometry reference (Schmid tensors, slip systems).
          */
         const SlipGeom &m_slipGeom;
         /**
          * @brief Lattice strain subproblem (fixed elastic state).
          * 
          * Provides stress evaluation and elastic properties but does not
          * update elastic strain during iteration.
          */
         const EvptnLatticeStrainProblem<ThermoElastN> m_lattice_strain_prob;
         /**
          * @brief Lattice rotation subproblem (active unknowns).
          * 
          * Template parameter ind_sub_r = 0 since rotation is the only unknown.
          */
         const EvptnLatticeRotationProblem<0> m_lattice_rot_prob;
         /**
          * @brief Inverse elastic strain rate scaling factor [time].
          * 
          * Conditioning factor: epsdot_scale_inv ≈ 1/‖D_sample‖
          */
         double m_epsdot_scale_inv;
         /**
          * @brief Inverse rotation increment scaling factor [time].
          * 
          * Conditioning factor: rotincr_scale_inv = (1/dt) * epsdot_scale_inv
          */
         double m_rotincr_scale_inv;
         /**
          * @brief Plastic deformation rate [1/time, ntvec = 5].
          * 
          * Computed once in constructor from beginning-of-step slip rates:
          * D_p = Σ_α γ̇^α P^α
          */
         double m_plastic_def_rate_d5[ecmech::ntvec];
         /**
          * @brief Plastic spin vector [rad/time, nwvec = 3].
          * 
          * Computed once in constructor from beginning-of-step slip rates:
          * W_p = Σ_α γ̇^α Q^α
          */
         double m_plastic_spin_vec[ecmech::nwvec];
         /**
          * @brief Prescribed deformation rate in sample frame [1/time, ntvec = 5].
          */
         const double* const m_def_rate_d5_sample;
         /**
          * @brief Prescribed spin in sample frame [rad/time, nwvec = 3].
          */
         const double* const m_spin_vec_sample;
         /**
          * @brief Jacobian array size (nDimSys × nDimSys = 9).
          */
         static constexpr int m_nXnDim = nDimSys * nDimSys;
      };

      /**
       * @brief Elastic strain-only update problem with fixed lattice rotation.
       * 
       * This class formulates a simplified integration problem where only the elastic
       * deviatoric strain is updated, while lattice rotation is held fixed. This
       * decoupled approach serves as a backup solver when the fully coupled
       * EvptnUpdstProblem fails to converge.
       * 
       * **Physical interpretation**:
       * Updates elastic strain ε_e due to:
       * - Prescribed deformation rate (rotated to crystal frame with fixed orientation)
       * - Plastic deformation from slip (stress-dependent, computed from current ε_e)
       * - Balance equation: ε̇_e = D_xtal - D_p
       * 
       * **When this is used**:
       * - Backup solver in getResponseRetry() after primary solver failure
       * - Alternative solver in getResponseNRSngl()
       * - Debugging/testing decoupled elastic integration
       * 
       * **Simplifications vs. EvptnUpdstProblem**:
       * - Updates: Elastic deviatoric strain (5 DOFs)
       * - Fixed: Lattice rotation (uses end-of-step values from prob_state.quat_u)
       * - Fixed: Hardening state (uses end-of-step values from prob_state.h_state_u)
       * - Result: Simpler 5×5 system instead of 8×8 coupled system
       * 
       * **System dimensions**:
       * - nDimSys = ntvec = 5 (elastic strain DOFs only)
       * - Unknowns: ε_e = deviatoric elastic strain in crystal frame
       * 
       * **Residual equation**:
       * ```
       * R_ε = ε̇_e - (1/a_vol) * (D_xtal - D_p)
       * ```
       * where D_xtal = R^T · D_sample · R uses fixed orientation R.
       * 
       * **Trade-offs**:
       * - More robust: simpler system often converges when coupled system fails
       * - Less accurate: neglects elastic-rotation coupling during step
       * - Still captures: stress-plasticity-hardening coupling
       * 
       * **Template parameters**:
       * @tparam SlipGeom Slip geometry class (e.g., SlipGeomFCC)
       * @tparam Kinetics Hardening kinetics class (e.g., KineticsKMBalD)
       * @tparam ThermoElastN Thermoelastic model (e.g., ThermoElastNCubic)
       * @tparam ProblemState State container (ProblemState<...>)
       * 
       * **Compile-time flags**:
       * - Requires: ECMECH_EXTRA_SOLVERS defined
       * - Optional: ECMECH_DEBUG for additional diagnostics
       * 
       * **Naming note**:
       * "NR" likely refers to "Newton-Raphson" for the nonlinear solve strategy,
       * distinguishing this from rotation-only (RotUpdProblem) or coupled updates.
       * 
       * @see EvptnUpdstProblem for fully coupled integration
       * @see RotUpdProblem for rotation-only updates
       * @see getResponseRetry() for usage in convergence recovery
       * @see getResponseNRSngl() for direct invocation
       */
      template<class SlipGeom, class Kinetics, class ThermoElastN, class ProblemState>
      class EvptnNRUpdstProblem
      {
         public:
           /**
            * @brief System dimension (number of DOFs).
            * 
            * Value: ntvec = 5 (deviatoric elastic strain components)
            */
            static constexpr int nDimSys = ecmech::ntvec;

           /**
            * @brief Constructor initializes elastic strain-only update problem.
            * 
            * Sets up the elastic strain integration problem using fixed orientation
            * from the end of the time step. Hardening state is updated implicitly
            * during iteration via the kinetics model.
            * 
            * **Initialization sequence**:
            * 1. Store references to constitutive models
            * 2. Create lattice strain and rotation subproblems
            * 3. Evaluate kinetic values at end-of-step hardening state
            * 4. Compute adaptive scaling factors for conditioning
            * 5. Store fixed orientation from prob_state.quat_u
            * 
            * **Key difference from EvptnUpdstProblem**:
            * Uses prob_state.quat_u (end-of-step, user-provided) instead of
            * solving for rotation. This assumes orientation is already known or
            * has been separately updated.
            * 
            * **Scaling strategy**:
            * Same adaptive approach as EvptnUpdstProblem:
            * ```
            * adots_ref = kinetics reference rate
            * eff = ‖D_sample‖
            * 
            * if eff < epsdot_scl_nzeff * adots_ref:
            *     epsdot_scale_inv = 1 / adots_ref
            * else:
            *     epsdot_scale_inv = min(1 / eff, 1e6 * dt)
            * 
            * rotincr_scale_inv = (1/dt) * epsdot_scale_inv
            * ```
            * Note: rotincr_scale_inv computed but not used (no rotation updates).
            * 
            * **Hardening state handling**:
            * Unlike RotUpdProblem, this problem can update hardening state during
            * iteration because hardening depends on slip rates, which depend on
            * stress, which depends on elastic strain (the active unknowns).
            * 
            * @param[in] slipGeom Slip system geometry (Schmid tensors)
            * @param[in] kinetics Hardening kinetics model
            *                     Provides getVals() for kinetic parameters
            * @param[in] thermoElastN Thermoelastic model
            * @param[in] prob_state Problem state container with:
            *                       - dt: Time step size
            *                       - h_state_u: End-of-step hardening state
            *                       - quat_u: End-of-step orientation (FIXED)
            *                       - elast_d5_n: Beginning-of-step elastic strain
            *                       - def_rate_d5_sample: Prescribed deformation rate
            *                       - spin_vec_sample: Prescribed spin (unused)
            *                       - Thermodynamic state: rel_vol_new, energy_new, pressure_EOS, tkelv
            * 
            * **Member initialization**:
            * - m_kin_vals: Kinetic parameters from end-of-step hardening
            * - m_hdn_scale: Hardening scaling factor
            * - m_lattice_strain_prob: Elastic state (active unknowns)
            * - m_lattice_rot_prob: Rotation (fixed, for reference)
            * - m_xtal_ori_quat: Fixed orientation from prob_state.quat_u
            * - m_epsdot_scale_inv, m_rotincr_scale_inv: Scaling factors
            * 
            * @note Orientation fixed at prob_state.quat_u throughout iteration
            * @note Hardening state can be updated via kinetics during solve
            * @note No plastic quantities pre-computed (evaluated during iteration)
            * 
            * @see EvptnLatticeStrainProblem for elastic strain handling
            * @see Kinetics::getVals() for kinetic parameter evaluation
            */
            __ecmech_hdev__
            EvptnNRUpdstProblem(const SlipGeom& slipGeom,
                                const Kinetics& kinetics,
                                const ThermoElastN& thermoElastN,
                                ProblemState& prob_state
                               )
                        :
                        m_slipGeom(slipGeom),
                        m_kinetics(kinetics),
                        m_thermoElastN(thermoElastN),
                        m_lattice_strain_prob(thermoElastN, prob_state.dt, prob_state.rel_vol_new, prob_state.energy_new, prob_state.pressure_EOS, prob_state.tkelv, prob_state.elast_d5_n),
                        m_lattice_rot_prob(prob_state.dt, prob_state.quat_n),
                        m_elast_d5_n(prob_state.elast_d5_n),
                        m_xtal_ori_quat(prob_state.quat_u),
                        m_def_rate_d5_sample(prob_state.def_rate_d5_sample), // vel_grad_sm%d_vecds
                        m_spin_vec_sample(prob_state.spin_vec_sample), // vel_grad_sm%w_veccp
                        m_mtan_sI(nullptr)
            {
               m_hdn_scale = m_kinetics.getVals(m_kin_vals, prob_state.pressure_EOS, prob_state.tkelv, prob_state.h_state_u);

               double adots_ref = m_kinetics.getFixedRefRate(m_kin_vals);
               double eff = vecNorm<ntvec>(m_def_rate_d5_sample); // do not worry about factor of sqrt(twothird)
               if (eff < (epsdot_scl_nzeff * adots_ref)) {
                  m_epsdot_scale_inv = one / adots_ref;
               }
               else {
                  m_epsdot_scale_inv = fmin(one / eff, 1e6 * m_lattice_strain_prob.m_dt);
               }
               //
               m_rotincr_scale_inv = m_lattice_strain_prob.m_inv_dt * m_epsdot_scale_inv;
            }

           /**
            * @brief Destructor (default).
            */
            __ecmech_hdev__
            ~EvptnNRUpdstProblem() {}

           /**
            * @brief Provide storage for material tangent stiffness computation.
            * 
            * Sets internal pointer to enable tangent stiffness calculation during
            * Jacobian evaluation. The tangent relates stress increments to strain
            * rate increments for finite element consistent linearization.
            * 
            * @param[in] mtan_sI Pointer to tangent stiffness array [nsvec² = 36]
            *                    Format: deviatoric-pressure (vecds) notation
            *                    Will be populated during computeRJ() if non-null
            * 
            * @see clearMTan() to reset pointer after tangent computation
            * @see computeTangentStiffness() for typical usage context
            * @see get_material_tangent_stiffness() for tangent evaluation
            */
            __ecmech_hdev__
            inline
            void provideMTan(double* mtan_sI) { m_mtan_sI = mtan_sI; }

           /**
            * @brief Clear material tangent pointer (disable tangent computation).
            * 
            * Resets internal pointer to nullptr, disabling tangent stiffness
            * computation in subsequent computeRJ() calls.
            * 
            * **Usage**: Call after tangent has been computed to avoid unnecessary
            * computation in subsequent residual-only evaluations.
            * 
            * @see provideMTan() to enable tangent computation
            */
            __ecmech_hdev__
            inline
            void clearMTan( ) { m_mtan_sI = nullptr; }

           /**
            * @brief Get inverse time step (1/Δt) for tangent scaling.
            */
            __ecmech_hdev__
            inline
            double getDtRi() const { return m_lattice_strain_prob.m_inv_dt; }

           /**
            * @brief Get hardening scaling factor.
            * 
            * Returns the scaling factor applied to hardening evolution, typically
            * related to temperature and pressure effects on kinetics.
            * 
            * @return Hardening scale factor [dimensionless]
            * 
            * @see Kinetics::getVals() for how this is computed
            */
            __ecmech_hdev__
            inline
            double getHdnScale() const { return m_hdn_scale; }

           /**
            * @brief Extract converged elastic strain from solution vector.
            * 
            * Converts the scaled solution vector (elastic strain increment in
            * solver space) back to physical elastic strain in crystal frame.
            * 
            * **Conversion sequence**:
            * 1. Unscale: Δε_e = e_scale * x
            * 2. Add to beginning-of-step: ε_e = Δε_e + ε_e^n
            * 
            * @param[out] elast_dev_press_vec End-of-step elastic strain [ntvec = 5]
            *                                  Deviatoric strain in crystal frame
            * @param[in] x Scaled solution vector [nDimSys = 5]
            *              Elastic strain increment in solver space
            * 
            * @note Safe to use same memory as prob_state.elast_d5_u
            * @note Only extracts elastic strain (rotation remains fixed)
            * 
            * @see EvptnLatticeStrainProblem::stateFromX for implementation
            */
            __ecmech_hdev__
            inline
            void stateFromX(double* const elast_dev_press_vec,
                            const double* const x) {
               m_lattice_strain_prob.stateFromX(elast_dev_press_vec,  &(x[m_i_sub_e]));
            }

           /**
            * @brief Compute Cauchy stress from elastic strain (convenience method).
            * 
            * Evaluates Cauchy stress in crystal frame given elastic strain,
            * accounting for finite deformation and thermodynamic state.
            * 
            * @param[out] cauchy_xtal Cauchy stress in crystal frame [nsvec = 6]
            *                          Voigt notation with pressure
            * @param[in] elast_d5_f Elastic deviatoric strain [ntvec = 5]
            *                        In crystal frame
            * 
            * **Stress evaluation**:
            * 1. Elastic law: σ_K = C : ε_e (Kirchhoff stress)
            * 2. Push-forward: σ = (1/J_e) * σ_K (Cauchy stress)
            * 3. Add EOS pressure: σ_vol = pressure_EOS
            * 
            * @see EvptnLatticeStrainProblem::elast_strain_to_cauchy_stress
            * @see ThermoElastN::eval for elastic law
            */
            __ecmech_hdev__
            inline
            void elastNEtoC(double* const cauchy_xtal, // nsvec
                            const double* const elast_d5_f // ntvec
                            ) const {
               m_lattice_strain_prob.elast_strain_to_cauchy_stress(cauchy_xtal, elast_d5_f);
            }

           /**
            * @brief Compute residual and Jacobian for SNLS trust region solver.
            * 
            * Evaluates the elastic strain evolution equation and its derivatives
            * with respect to elastic strain unknowns, using fixed orientation.
            * This is the core computation for the implicit integration of elastic strain.
            * 
            * **Residual formulation**:
            * ```
            * R_ε = (ε̇_e - (1/a_vol) * (D_xtal - D_p)) * epsdot_scale_inv
            * ```
            * where:
            * - ε̇_e = (ε_e - ε_e^n) / Δt (elastic strain rate)
            * - D_xtal = R^T · D_sample · R (deformation rate in crystal frame, R fixed)
            * - D_p = Σ_α γ̇^α P^α (plastic deformation rate from slip)
            * 
            * At convergence: R_ε = 0, implying ε̇_e = (1/a_vol) * (D_xtal - D_p)
            * 
            * **Computation sequence**:
            * 
            * 1. **Initialize arrays**:
            *    - Zero Jacobian (if requested)
            *    - Zero residual
            * 
            * 2. **Extract elastic strain state** from solution vector x:
            *    - Unscale and add to beginning-of-step
            *    - ε_e^{n+1} = e_scale * x + ε_e^n
            *    - Compute ε̇_e = (ε_e^{n+1} - ε_e^n) / Δt
            * 
            * 3. **Transform deformation rate** to crystal frame (fixed rotation):
            *    - Use m_xtal_ori_quat to generate rotation matrices
            *    - D_xtal = R₅^T · D_sample · R₅
            *    - W_xtal = R^T · W_sample · R (spin, unused for elastic-only)
            * 
            * 4. **Evaluate stress** from elastic strain:
            *    - σ_K = C : ε_e (Kirchhoff)
            *    - Via lattice_strain_prob.elast_strain_to_kirchoff_stress()
            * 
            * 5. **Compute slip rates** from stress and hardening:
            *    - Resolve stress onto slip systems: τ^α = P^α : σ_K
            *    - Kinetics: γ̇^α = f(τ^α, h, T)
            *    - Also compute ∂γ̇^α/∂τ^α for Jacobian
            * 
            * 6. **Plastic deformation rate**:
            *    - D_p = Σ_α γ̇^α P^α
            *    - W_p = Σ_α γ̇^α Q^α (plastic spin, unused)
            * 
            * 7. **Evaluate elastic strain residual**:
            *    - R_ε = lattice_strain_prob.get_elast_strain_residual(...)
            *    - R_ε = ε̇_e - (1/a_vol) * (D_xtal - D_p)
            * 
            * 8. **Jacobian** (if requested):
            *    - Compute chain rule derivatives:
            *      * ∂D_p/∂ε_e via ∂γ̇/∂τ and ∂τ/∂ε_e (elastic stiffness)
            *    - J = ∂R_ε/∂ε_e = (1/Δt) * I + (1/a_vol) * ∂D_p/∂ε_e
            *    - Apply scaling: epsdot_scale_inv * e_scale
            * 
            * 9. **Material tangent** (if m_mtan_sI != nullptr):
            *    - Compute dσ/dD via implicit function theorem
            *    - Uses converged Jacobian and rotation derivatives
            *    - See get_material_tangent_stiffness() for details
            * 
            * **Jacobian structure** (5×5):
            * ```
            * [∂R₁/∂ε₁  ∂R₁/∂ε₂  ∂R₁/∂ε₃  ∂R₁/∂ε₄  ∂R₁/∂ε₅]
            * [∂R₂/∂ε₁  ∂R₂/∂ε₂  ∂R₂/∂ε₃  ∂R₂/∂ε₄  ∂R₂/∂ε₅]
            * [∂R₃/∂ε₁  ∂R₃/∂ε₂  ∂R₃/∂ε₃  ∂R₃/∂ε₄  ∂R₃/∂ε₅]
            * [∂R₄/∂ε₁  ∂R₄/∂ε₂  ∂R₄/∂ε₃  ∂R₄/∂ε₄  ∂R₄/∂ε₅]
            * [∂R₅/∂ε₁  ∂R₅/∂ε₂  ∂R₅/∂ε₃  ∂R₅/∂ε₄  ∂R₅/∂ε₅]
            * ```
            * 
            * **Material tangent computation** (optional):
            * If m_mtan_sI != nullptr, additionally computes material tangent stiffness
            * relating stress to deformation rate. This requires:
            * - Converged state (small residual)
            * - Jacobian evaluation
            * - Higher-order elastic terms
            * - Rotation derivatives
            * 
            * The tangent is computed via get_material_tangent_stiffness() and stored
            * in the provided array for subsequent finite element assembly.
            * 
            * @param[out] resid Residual vector [nDimSys = 5]
            *                   Elastic strain evolution equation
            *                   Set to zero, then populated
            * @param[out] Jacobian Jacobian matrix [nDimSys × nDimSys = 25]
            *                      Row-major storage: J[i*nDimSys + j] = ∂R_i/∂x_j
            *                      Set to zero, then populated
            *                      If nullptr, Jacobian computation skipped
            * @param[in] x Solution vector [nDimSys = 5]
            *              Scaled elastic strain increment
            *              Solver iterates to find ε_e satisfying R(ε_e) = 0
            * 
            * @return true Always returns true (elastic strain update typically converges)
            * 
            * **Performance notes**:
            * - Jacobian computation optional (nullptr check)
            * - Material tangent computation optional (m_mtan_sI check)
            * - Most expensive: slip rate derivatives and chain rule
            * - Fixed rotation avoids rotation derivative computation
            * - No dynamic allocation
            * 
            * @note resid and Jacobian zeroed at start
            * @note Jacobian only computed if Jacobian != nullptr
            * @note Material tangent only computed if m_mtan_sI != nullptr
            * @note Uses fixed orientation m_xtal_ori_quat throughout
            * @note Thread-safe (no shared state modified)
            * 
            * @see EvptnLatticeStrainProblem::get_elast_strain_residual for R_ε
            * @see get_slip_rate_terms() for plastic rate computation
            * @see get_slip_rate_deriv_terms() for Jacobian chain rule
            * @see get_material_tangent_stiffness() for tangent computation
            */
            __ecmech_hdev__
            bool computeRJ(double* const resid,
                           double* const Jacobian,
                           const double* const x) {

               bool doComputeJ = (Jacobian != nullptr);

               if (doComputeJ) {
                  // zero the Jacobian so that do not need to worry about zero
                  // entries in the midst of other things later
                  //
                  for (size_t ijJ = 0; ijJ<m_nXnDim; ++ijJ) {
                     Jacobian[ijJ] = 0.0;
                  }
               }
               //
               for (int iR = 0; iR<nDimSys; ++iR) {
                  resid[iR] = 0.0;
               }

               //////////////////////////////
               // PULL VALUES out of x, with scalings
               //
               double elast_dt_d5[ecmech::ntvec];
               vecsVxa<ntvec>(elast_dt_d5, ecmech::e_scale, &(x[m_i_sub_e]) ); // elast_dt_d5 is now the delta, _not_ yet elast_dt_d5
               // elast_d5_f is end-of-step
               double elast_d5_f[ntvec];
               vecsVapb<ntvec>(elast_d5_f, elast_dt_d5, m_lattice_strain_prob.m_elast_d5_n);
               vecsVsa<ntvec>(elast_dt_d5, m_lattice_strain_prob.m_inv_dt); // _now_ elast_dt_d5 has elast_dt_d5
               //
               double xtal_rmat[ecmech::ndim * ecmech::ndim];
               quat_to_tensor(xtal_rmat, m_xtal_ori_quat);
               //
               double rmat_5x5_sample2xtal[ecmech::ntvec * ecmech::ntvec];
               get_rot_mat_vecd(rmat_5x5_sample2xtal, xtal_rmat);
               double def_rate_d5_xtal[ecmech::ntvec];
               double spin_vec_xtal[ecmech::nwvec]; // assumes nwvec = ndim

               get_xtal_frame_vel_grad_terms(def_rate_d5_xtal, spin_vec_xtal, m_def_rate_d5_sample,  m_spin_vec_sample, xtal_rmat, rmat_5x5_sample2xtal);

               //////////////////////////////
               // CALCULATIONS

               double kirchoff[ecmech::nsvec];
               m_lattice_strain_prob.elast_strain_to_kirchoff_stress(kirchoff, elast_d5_f);

               double dgdot_dtau[SlipGeom::nslip] = { 0.0 }; // crys%tmp2_slp
               double plastic_def_rate_d5[ecmech::ntvec] = { 0.0 };
               double plastic_spin_vec[ecmech::nwvec] = { 0.0 }; // \pcDhat

               get_slip_rate_terms(dgdot_dtau, plastic_def_rate_d5, plastic_spin_vec, kirchoff, m_kin_vals, m_slipGeom, m_kinetics);

               // Residual Calculations
               m_lattice_strain_prob.get_elast_strain_residual(resid, m_epsdot_scale_inv, elast_dt_d5, plastic_def_rate_d5, def_rate_d5_xtal);

               //////////////////////////////////////////////////////////////////////
               // JACOBIAN, fixed hardness and temperature
               //
               if (doComputeJ) {
                  //
                  // preliminaries
                  //
                  double dpl_deps_symm[ ecmech::ntvec * ecmech::ntvec ] = { 0.0 };
                  double dpl_deps_skew[ ecmech::nwvec * ecmech::ntvec ] = { 0.0 };
                  get_slip_rate_deriv_terms(dpl_deps_symm, dpl_deps_skew, dgdot_dtau, m_lattice_strain_prob.m_inv_a_vol, m_slipGeom, m_lattice_strain_prob.m_thermo_elast_n);

                  // d(B_S)/d(elast_d5_f)
                  //
                  m_lattice_strain_prob.template get_deriv_elast_strain_wrt_elast_strain<nDimSys>(Jacobian, dpl_deps_symm);

                  if (m_mtan_sI) {

                     // Higher-order terms related to the elasticity stuff that's used in the residuals and jacobian calculation
                     double A_e_M35[ecmech::nwvec * ecmech::ntvec];
                     double ee_wvec[ecmech::nwvec];
                     double ee_fac;
                     double xi_f[ecmech::nwvec] = {};

                     m_lattice_rot_prob.deltaOmegaFromState(xi_f, m_xtal_ori_quat, m_lattice_rot_prob.m_xtal_ori_quat_n);

                     elasticity_higher_order_terms(A_e_M35, ee_wvec, ee_fac, m_lattice_strain_prob.m_inv_a_vol, elast_d5_f, elast_dt_d5);

                     // derivatives with respect to lattice orientation changes
                     double dxtal_ori_quat_dxi_T[ ecmech::nwvec * ecmech::qdim ];
                     double dDsm_dxi[ ecmech::ntvec * ecmech::nwvec ];
                     double dWsm_dxi[ ecmech::nwvec * ecmech::nwvec ];
                     eval_d_dxi_impl_quat(dxtal_ori_quat_dxi_T, dDsm_dxi, dWsm_dxi,
                                          m_def_rate_d5_sample, m_spin_vec_sample,
                                          xi_f,
                                          m_lattice_rot_prob.m_xtal_ori_quat_n,
                                          xtal_rmat, m_xtal_ori_quat);

                     static constexpr int nDimSolve = nDimSys + ecmech::nwvec;
                     static constexpr int nDimSolve2 = nDimSolve * nDimSolve;

                     double Jacobian2[nDimSolve2] = {};

                     RAJA::View<double, RAJA::Layout<2> > pfrac_ee(Jacobian2, nDimSolve, nDimSolve);
                     RAJA::View<double, RAJA::Layout<2> > jacob_ee(Jacobian, nDimSys, nDimSys);
                     for (int i_jac = 0; i_jac < nDimSys; i_jac++) {
                        for (int j_jac = 0; j_jac < nDimSys; j_jac++)
                        pfrac_ee(i_jac, j_jac) = jacob_ee(i_jac, j_jac);
                     }

                     // d(B_S)/d(xi_f)
                     //
                     // jacob_er = -dDsm_dxi(:,:)
                     m_lattice_rot_prob.template get_deriv_elast_strain_wrt_omega<nDimSolve>(Jacobian2, dDsm_dxi);
                     // d(B_xi)/d(elast_dev_press_vecs_f)
                     //
                     m_lattice_strain_prob.template get_deriv_omega_wrt_elast_strain<nDimSolve, m_i_sub_r>(Jacobian2, elast_dt_d5, ee_fac, dpl_deps_skew, A_e_M35);
                     // d(B_xi)/d(xi_f)
                     //
                     m_lattice_rot_prob.template get_deriv_omega_wrt_omega<nDimSolve>(Jacobian2, dWsm_dxi);

                     double cauchy_stress_lattice[ ecmech::nsvec ];
                     m_lattice_strain_prob.m_thermo_elast_n.getCauchy(cauchy_stress_lattice, kirchoff, m_lattice_strain_prob.m_inv_det_v_e);
                     get_material_tangent_stiffness<ThermoElastN, nDimSolve, m_i_sub_r>(m_mtan_sI, Jacobian2,
                     dxtal_ori_quat_dxi_T, rmat_5x5_sample2xtal,
                     m_xtal_ori_quat, xtal_rmat,
                     cauchy_stress_lattice,
                     m_lattice_strain_prob.m_inv_det_v_e,
                     m_lattice_strain_prob.m_inv_a_vol,
                     m_lattice_strain_prob.m_thermo_elast_n);
                  }

                  // SCALING
                  {
                     double scaleFactorJ;
                     for (size_t iJ = 0; iJ<m_i_sub_r; ++iJ) {
                        // Jacobian(i_sub_e:i_sup_e,i_sub_e:i_sup_e) = jacob_ee * epsdot_scale_inv  * e_scale ! resid, x
                        scaleFactorJ = m_epsdot_scale_inv * ecmech::e_scale;
                        for (size_t jJ = 0; jJ<m_i_sub_r; ++jJ) { // <=_i_sup_e
                           int ijJ = ECMECH_NN_INDX(iJ, jJ, nDimSys);
                           Jacobian[ ijJ ] *= scaleFactorJ;
                        }
                     }
                  } // SCALING
               }

               return true;
            } // computeRJ

            __ecmech_hdev__
            inline
            void get_slip_contribution(double& pl_disipation_rate,
                                       double& effective_shear_rate,
                                       double* const gdot,
                                       const double* const elast_strain
                                      )
            {
               get_slip_contributions(pl_disipation_rate, effective_shear_rate, gdot,
                                      m_lattice_strain_prob.m_inv_det_v_e, elast_strain, m_kin_vals,
                                      m_slipGeom, m_kinetics, m_lattice_strain_prob);
            }

         private:
            /**
             * @brief Slip geometry reference.
             */
            const SlipGeom &m_slipGeom;
            /**
             * @brief Hardening kinetics reference.
             */
            const Kinetics &m_kinetics;
            /**
             * @brief Thermoelastic model reference.
             */
            const ThermoElastN &m_thermoElastN;
            /**
             * @brief Lattice strain subproblem (active unknowns).
             */
            const EvptnLatticeStrainProblem<ThermoElastN> m_lattice_strain_prob;
            /**
             * @brief Lattice rotation subproblem (fixed, for reference).
             * 
             * Template parameter ind_sub_r = ntvec since rotation would follow
             * elastic strain in a coupled system (but here rotation is fixed).
             */
            const EvptnLatticeRotationProblem<ecmech::ntvec> m_lattice_rot_prob;
            /**
             * @brief Hardening scale factor [dimensionless].
             */
            double m_hdn_scale;
            /**
             * @brief Inverse elastic strain rate scaling [time].
             */
            double m_epsdot_scale_inv;
            /**
             * @brief Inverse rotation increment scaling [time].
             * 
             * Computed but not used (no rotation updates in this problem).
             */
            double m_rotincr_scale_inv;
            /**
             * @brief Kinetic values array [Kinetics::nVals].
             * 
             * Contains slip resistance, reference rates, etc. from kinetics model.
             */
            double m_kin_vals[Kinetics::nVals];
            /**
             * @brief Beginning-of-step elastic strain [ntvec = 5].
             */
            const double* const m_elast_d5_n;
            /**
             * @brief Fixed orientation quaternion [qdim = 4].
             * 
             * Taken from prob_state.quat_u, remains constant during iteration.
             */
            const double* const m_xtal_ori_quat;
            /**
             * @brief Prescribed deformation rate in sample frame [ntvec = 5].
             */
            const double* const m_def_rate_d5_sample;
            /**
             * @brief Prescribed spin in sample frame [nwvec = 3].
             * 
             * Stored but not used in elastic-only formulation.
             */
            const double* const m_spin_vec_sample;
            /**
             * @brief Jacobian array size (nDimSys × nDimSys = 25).
             */
            static constexpr int m_nXnDim = nDimSys * nDimSys;
            /**
             * @brief Starting index for elastic strain DOFs (always 0).
             */
            static constexpr int m_i_sub_e = 0;
            /**
             * @brief Starting index for rotation DOFs (ntvec = 5).
             * 
             * Would be used in coupled system; here for consistency.
             */
            static constexpr int m_i_sub_r = ecmech::ntvec;
            /**
             * @brief Material tangent stiffness pointer [nsvec² = 36] or nullptr.
             * 
             * If non-null, material tangent computed during computeRJ().
             */
            double* m_mtan_sI;
      }; // class EvptnNRUpdstProblem

#endif

   } // namespace evptn
} // namespace ecmech

#endif // ECMECH_EVPTN_H
