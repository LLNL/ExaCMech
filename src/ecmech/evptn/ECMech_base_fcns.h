/**
 * @file ECMech_base_fcns.h
 * @brief Utility functions for crystal plasticity kinematics and slip system calculations.
 * 
 * This header provides a collection of template functions that implement core
 * kinematic operations for crystal plasticity finite element models:
 * 
 * **Frame transformations**:
 * - Rotate velocity gradient components from sample to crystal frame
 * - Apply 3×3 rotation matrices to vectors (spin)
 * - Apply 5×5 rotation matrices to deviatoric tensors (deformation rate)
 * 
 * **Slip system mechanics**:
 * - Resolve shear stress on slip systems (Schmid law)
 * - Evaluate slip rates from resolved stresses (kinetic laws)
 * - Compute plastic deformation and spin from slip rates
 * - Calculate slip rate derivatives for Jacobians
 * 
 * **Post-processing utilities**:
 * - Plastic dissipation rate
 * - Effective shear rate
 * - Slip contributions to work
 * 
 * **Higher-order coupling**:
 * - Elastic spin contributions (finite strain effects)
 * - Frame transformation derivatives
 * - Material tangent stiffness helpers
 * 
 * **Design philosophy**:
 * - Header-only implementation (all functions inline/template)
 * - GPU-compatible (__ecmech_hdev__ annotation)
 * - Zero dynamic allocation (stack arrays only)
 * - Compile-time dispatch via template parameters (static vs. dynamic slip systems)
 * 
 * **Usage pattern**:
 * ```cpp
 * // Example: Compute plastic deformation from slip
 * double gdot[nslip];  // Slip rates
 * double Dp[ntvec];    // Output: plastic def rate
 * double Wp[nwvec];    // Output: plastic spin
 * get_slip_rate_terms(dgdot_dtau, Dp, Wp, kirchoff, kin_vals, slip_geom, kinetics);
 * ```
 * 
 * @see ECMech_base_classes.h for problem formulation classes using these functions
 * @see ECMech_evptn.h for integration of functions into residual/Jacobian evaluations
 * @see ECMech_util.h for lower-level tensor operation primitives
 */

#pragma once

#include "ECMech_core.h"
#include "ECMech_util.h"

#include "SNLS_lup_solve.h"

namespace ecmech {
namespace evptn {

    /**
     * @brief Transform velocity gradient components from sample to crystal frame.
     * 
     * This function rotates the prescribed velocity gradient (deformation rate + spin)
     * from the sample/laboratory frame to the crystal/lattice frame using the current
     * crystal orientation. This is essential for crystal plasticity because:
     * - Slip systems are defined in the crystal frame
     * - Resolved shear stresses need crystal frame stress
     * - Plastic rates are computed in crystal frame
     * 
     * **Transformation formulas**:
     * ```
     * D_crystal = R^T * D_sample * R  (for symmetric deformation rate)
     * W_crystal = R^T * W_sample * R  (for skew spin tensor)
     * ```
     * where R is the rotation matrix from sample to crystal (from quaternion).
     * 
     * **Tensor representations**:
     * - Deformation rate: 5-vector deviatoric form, needs 5×5 rotation matrix
     * - Spin: 3-vector axial form, needs 3×3 rotation matrix
     * 
     * @param[out] def_rate_d5_xtal Deformation rate in crystal frame [1/time, 5 components]
     * @param[out] spin_vec_xtal Spin vector in crystal frame [1/time, 3 components]
     * @param[in] def_rate_d5_sample Deformation rate in sample frame [1/time, 5 components]
     * @param[in] spin_vec_sample Spin vector in sample frame [1/time, 3 components]
     * @param[in] xtal_rmat Rotation matrix sample→crystal [3×3, row-major]
     * @param[in] xtal_rot_mat5 Rotation matrix for 5-vectors [5×5, row-major]
     * 
     * **Input rotation matrices**:
     * - xtal_rmat: Standard 3×3 rotation from quat_to_tensor()
     * - xtal_rot_mat5: 5×5 deviatoric rotation from get_rot_mat_vecd()
     * 
     * **Usage in integration**:
     * Called at each iteration to rotate prescribed kinematics to crystal frame
     * before computing slip system activity.
     * 
     * @note Both output arrays must be pre-allocated (no bounds checking)
     * @note Function assumes rotation matrices are orthogonal (no validation)
     * 
     * @see quat_to_tensor() for generating xtal_rmat from quaternion
     * @see get_rot_mat_vecd() for generating xtal_rot_mat5 from 3×3 rotation
     */
    __ecmech_hdev__
    inline
    void get_xtal_frame_vel_grad_terms(double* const def_rate_d5_xtal,
                                        double* const spin_vec_xtal,
                                        const double* const def_rate_d5_sample,
                                        const double* const spin_vec_sample,
                                        const double* const xtal_rmat,
                                        const double* const xtal_rot_mat5)
    {
        vecsVMTa<ecmech::ntvec>(def_rate_d5_xtal, xtal_rot_mat5, def_rate_d5_sample);
        vecsVMTa<ecmech::ndim>(spin_vec_xtal, xtal_rmat, spin_vec_sample);
    }

    /**
     * @brief Compute slip rates and plastic deformation from resolved stresses.
     * 
     * This template function implements the core crystal plasticity kinematic calculation:
     * 1. Resolve Kirchhoff stress onto slip systems → resolved shear stresses τ^α
     * 2. Evaluate slip rates via kinetic law: γ̇^α = f(τ^α, g^α, T)
     * 3. Sum slip contributions to plastic deformation: D_p = Σ γ̇^α P^α
     * 4. Sum slip contributions to plastic spin: W_p = Σ γ̇^α Q^α
     * 
     * **Constitutive sequence**:
     * ```
     * τ^α = P^α : τ                    (Schmid's law: resolve stress)
     * γ̇^α = f(τ^α)                     (Slip kinetics laws)
     * D_p = Σ_α γ̇^α * P^α              (Plastic deformation rate)
     * W_p = Σ_α γ̇^α * Q^α              (Plastic spin)
     * ```
     * 
     * **Template parameters**:
     * @tparam SlipGeom Slip geometry class (e.g., SlipGeomFCC, SlipGeomBCC)
     *                  Must provide: nslip, P (Schmid tensors), Q (spin tensors), evalRSS(), getExtras()
     * @tparam SlipKinetics Kinetics class (e.g., KineticsKMBalD, KineticsVocePL)
     *                      Must provide: evalGdots()
     * 
     * @param[out] dgdot_dtau Slip rate derivatives ∂γ̇^α/∂τ^α [1/stress, nslip components]
     *                        Used for Jacobian evaluation (zero if kinetics is rate-independent)
     * @param[out] plastic_def_rate_d5 Plastic deformation rate D_p [1/time, 5 components]
     * @param[out] plastic_spin_vec Plastic spin W_p [1/time, 3 components]
     * @param[in] kirchoff Kirchhoff stress tensor [stress units, 6 components Voigt]
     * @param[in] kinetic_values Material state values for kinetics [various units]
     *                           Typically: [slip_resistances, dislocation_densities, ...]
     * @param[in] slip_geom Reference to slip geometry object
     * @param[in] slip_kinetics Reference to kinetics object
     * 
     * **Dynamic slip systems**:
     * If SlipGeom::dynamic is true, additional "extra" quantities (e.g., chi angles for
     * non-Schmid effects) are evaluated via getExtras() and appended to RSS array.
     * 
     * **Compile-time optimization**:
     * - if constexpr (SlipGeom::nslip > 0) allows zero-slip geometries (elastic only)
     * - Entire function body skipped at compile time if nslip=0, avoiding unused variables
     * 
     * @note All output arrays must be pre-allocated
     * @note Function does not handle temperature explicitly (passed via kinetic_values)
     * @note Kirchhoff stress used (not Cauchy) for consistent finite strain formulation
     * 
     * @see slip_geom.evalRSS() for resolved shear stress calculation
     * @see slip_kinetics.evalGdots() for slip rate evaluation
     */
    template<class SlipGeom, class SlipKinetics>
    __ecmech_hdev__
    inline
    void get_slip_rate_terms(double* const dgdot_dtau,
                            double* const plastic_def_rate_d5,
                            double* const plastic_spin_vec,
                            const double* const kirchoff,
                            const double* const kinetic_values,
                            const SlipGeom& slip_geom,
                            const SlipKinetics& slip_kinetics
                            )
    {        
        if constexpr (SlipGeom::nslip > 0) {
        // default initialize everything to 0.0
        constexpr size_t nslip_dyn = (SlipGeom::dynamic) ? (SlipGeom::nslip + SlipGeom::nSlipExtra) : SlipGeom::nslip;
        double abs_resolved_shear_stress[nslip_dyn] = {};
        double gdot[SlipGeom::nslip] = {};
        // resolve stress onto slip systems
        // CALL resolve_tau_a_n(crys%tmp4_slp, s_meas%kirchoff, crys)
        //vecsVaTM<ntvec, SlipGeom::nslip>(taua, kirchoff, slipP);
        slip_geom.evalRSS(abs_resolved_shear_stress, kirchoff, slip_geom.getP());
        if constexpr (SlipGeom::dynamic) {
            slip_geom.getExtras(&abs_resolved_shear_stress[SlipGeom::nslip]);
        }
        //
        // CALL plaw_eval(plastic_def_rate_d5, plastic_spin_vec, gss, crys, tkelv, ierr)
        // chi values are passed within extended taua array
        slip_kinetics.evalGdots(gdot, dgdot_dtau, abs_resolved_shear_stress, kinetic_values);
        
        //
        // CALL sum_slip_def(plastic_def_rate_d5, plastic_spin_vec, crys%tmp1_slp, crys) ;
        vecsVMa<ntvec, SlipGeom::nslip>(plastic_def_rate_d5, slip_geom.getP(), gdot);
        vecsVMa<nwvec, SlipGeom::nslip>(plastic_spin_vec, slip_geom.getQ(), gdot);
        }
    }

    /**
     * @brief Compute derivatives of plastic rates with respect to elastic strain.
     * 
     * This function performs the chain rule to obtain sensitivities needed for
     * implicit Jacobian evaluation:
     * ```
     * ∂D_p/∂ε_e = (∂D_p/∂γ̇)(∂γ̇/∂τ)(∂τ/∂ε_e)
     * ∂W_p/∂ε_e = (∂W_p/∂γ̇)(∂γ̇/∂τ)(∂τ/∂ε_e)
     * ```
     * 
     * **Derivative chain breakdown**:
     * 1. **∂τ/∂ε_e**: Elastic stiffness (from ThermoElastN::multDTDepsT)
     * 2. **∂γ̇/∂τ**: Kinetic rate sensitivity (from kinetics.evalGdots, dgdot_dtau)
     * 3. **∂D_p/∂γ̇**: Schmid tensor P^α (slip_geom.getP)
     * 4. **∂W_p/∂γ̇**: Spin tensor Q^α (slip_geom.getQ)
     * 
     * **Mathematical formulation**:
     * ```
     * dtaua_deps = (∂τ^α/∂ε_e) = P^α : K                 [nslip × ntvec]
     * dgdot_deps = dgdot_dtau * dtaua_deps               [nslip × ntvec]
     * dDp_deps = P^α * dgdot_deps^T                      [ntvec × ntvec]
     * dWp_deps = Q^α * dgdot_deps^T                      [nwvec × ntvec]
     * ```
     * 
     * **Template parameters**:
     * @tparam SlipGeom Slip geometry class providing Schmid and spin tensors
     * @tparam ThermoElastN Thermoelastic model providing stress derivatives
     * 
     * @param[out] dDp_hat_delast_strain Derivative ∂D_p/∂ε_e [ntvec × ntvec matrix]
     * @param[out] dWp_hat_delast_strain Derivative ∂W_p/∂ε_e [nwvec × ntvec matrix]
     * @param[in] dgdot_dtau Slip rate derivatives ∂γ̇^α/∂τ^α [1/stress, nslip]
     * @param[in] inv_a_vol Inverse volume scaling J^(-1/3) [dimensionless]
     * @param[in] slip_geom Reference to slip geometry object
     * @param[in] thermoElastN Reference to thermoelastic model
     * 
     * **Volume scaling**:
     * - inv_a_vol accounts for finite strain effects in stress-strain relation
     * - Ensures derivatives are consistent with volumetric changes
     * 
     * **Usage in Jacobian**:
     * These derivatives appear in off-diagonal blocks coupling elastic strain
     * to other equations (rotation, hardening).
     * 
     * **Compile-time optimization**:
     * - if constexpr (SlipGeom::nslip > 0) allows zero-slip case
     * - Avoids unused variable warnings for elastic-only materials
     * 
     * @note Output arrays must be pre-allocated [ntvec*ntvec] and [nwvec*ntvec]
     * @note Row-major storage: element (i,j) at index [i*ncols + j]
     * 
     * @see ThermoElastN::multDTDepsT() for elastic stiffness multiplication
     * @see vecsMABT() for matrix-matrix transpose products
     */
    template<class SlipGeom, class ThermoElastN>
    __ecmech_hdev__
    inline
    void get_slip_rate_deriv_terms(double* const dDp_hat_delast_strain,
                                    double* const dWp_hat_delast_strain,
                                    const double* const dgdot_dtau,
                                    const double inv_a_vol,
                                    const SlipGeom& slip_geom,
                                    const ThermoElastN& thermoElastN
                                )
    {
        if constexpr (SlipGeom::nslip > 0) {
        double dtaua_deps[ ecmech::ntvec * SlipGeom::nslip ];
        thermoElastN.multDTDepsT(dtaua_deps, slip_geom.getP(), inv_a_vol, SlipGeom::nslip);

        double dgdot_deps[ ecmech::ntvec * SlipGeom::nslip ];
        for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
            for (int iSlip = 0; iSlip < SlipGeom::nslip; ++iSlip) {
                int ijThis = ECMECH_NM_INDX(iTvec, iSlip, ecmech::ntvec, SlipGeom::nslip);
                dgdot_deps[ijThis] = dgdot_dtau[iSlip] * dtaua_deps[ijThis];
            }
        }
        vecsMABT<ntvec, SlipGeom::nslip>(dDp_hat_delast_strain, slip_geom.getP(), dgdot_deps);
        vecsMABT<nwvec, ntvec, SlipGeom::nslip>(dWp_hat_delast_strain, slip_geom.getQ(), dgdot_deps);
        }
    }

    /**
     * @brief Derivatives of plastic rates w.r.t. the spherical elastic strain.
     *
     * The plastic rates respond to the spherical elastic strain through the
     * resolved shear stresses:
     * ```
     * ∂γ̇^α/∂ε_s = (∂γ̇/∂τ)^α · P^α : (∂τ'/∂ε_s)
     * ∂D_p/∂ε_s = Σ_α P^α ∂γ̇^α/∂ε_s ,   ∂W_p/∂ε_s = Σ_α Q^α ∂γ̇^α/∂ε_s
     * ```
     * with ∂τ'/∂ε_s from ThermoElastN::getDTDepsSph(), which carries both the
     * symmetry coupling (K_sdax3, hexagonal and lower symmetry — elastic scale)
     * and the strain-scaling geometric term (all symmetries — stress scale).
     * These feed the volumetric right-hand side of the consistent-tangent solve
     * in get_material_tangent_stiffness(); without them the material tangent
     * misses the plastic feedback of volumetric loading, degrading host Newton
     * convergence.
     *
     * @param[out] dDp_hat_deps_sph Derivative ∂D_p/∂ε_s [ntvec]
     * @param[out] dWp_hat_deps_sph Derivative ∂W_p/∂ε_s [nwvec]
     * @param[in] dgdot_dtau Slip rate derivatives ∂γ̇^α/∂τ^α [nslip]
     * @param[in] elast_dev_press_vec Scaled strain vector [nsvec], as eval() takes
     * @param[in] energy_vol_ref Internal energy (Grüneisen sensitivity)
     * @param[in] slip_geom Slip geometry (Schmid/spin tensors)
     * @param[in] thermoElastN Thermoelastic model (provides getDTDepsSph)
     *
     * @see get_slip_rate_deriv_terms() for the deviatoric analogue
     * @see ThermoElastNHexag::getDTDepsSph() for the sensitivity vector
     */
    template<class SlipGeom, class ThermoElastN>
    __ecmech_hdev__
    inline
    void get_slip_rate_deriv_sph_terms(double* const dDp_hat_deps_sph,
                                       double* const dWp_hat_deps_sph,
                                       const double* const dgdot_dtau,
                                       const double* const elast_dev_press_vec,
                                       const double energy_vol_ref,
                                       const SlipGeom& slip_geom,
                                       const ThermoElastN& thermoElastN
                                    )
    {
        for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
            dDp_hat_deps_sph[iTvec] = 0.0;
        }
        for (int iWvec = 0; iWvec < ecmech::nwvec; ++iWvec) {
            dWp_hat_deps_sph[iWvec] = 0.0;
        }
        double dT_deps_sph[ecmech::ntvec];
        if (!thermoElastN.getDTDepsSph(dT_deps_sph, elast_dev_press_vec, energy_vol_ref)) {
            return;
        }
        if constexpr (SlipGeom::nslip > 0) {
        double dgdot_deps_sph[SlipGeom::nslip];
        const double* const slipP = slip_geom.getP();
        for (int iSlip = 0; iSlip < SlipGeom::nslip; ++iSlip) {
            double dtaua_deps_sph = 0.0;
            for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                dtaua_deps_sph += slipP[ECMECH_NM_INDX(iTvec, iSlip, ecmech::ntvec, SlipGeom::nslip)] * dT_deps_sph[iTvec];
            }
            dgdot_deps_sph[iSlip] = dgdot_dtau[iSlip] * dtaua_deps_sph;
        }
        vecsVMa<ntvec, SlipGeom::nslip>(dDp_hat_deps_sph, slip_geom.getP(), dgdot_deps_sph);
        vecsVMa<nwvec, SlipGeom::nslip>(dWp_hat_deps_sph, slip_geom.getQ(), dgdot_deps_sph);
        }
    }

    /**
     * @brief Calculate slip-related post-processing quantities.
     * 
     * This function computes diagnostic/output quantities derived from the current
     * slip state but not strictly part of the state vector:
     * - Plastic dissipation rate: Mechanical work rate from plastic deformation
     * - Effective shear rate: Scalar measure of plastic deformation intensity
     * - Slip rates on each system: Individual γ̇^α values
     * 
     * **Physical meanings**:
     * 
     * **Plastic dissipation rate** [power/volume]:
     * ```
     * Ḋ_p = (1/J) * Σ_α |τ^α| * |γ̇^α|
     * ```
     * Represents irreversible energy dissipation from plastic work.
     * Used for:
     * - Energy balance in thermomechanical coupling
     * - Damage models (plastic work as nucleation criterion)
     * - Validation (should match stress:plastic_rate)
     * 
     * **Effective shear rate** [1/time]:
     * ```
     * γ̇_eff = √(2/3 * D_p : D_p)     (if ECMECH_USE_DPEFF defined)
     * γ̇_eff = Σ_α |γ̇^α|              (otherwise)
     * ```
     * Scalar measure of plastic deformation intensity.
     * Used for:
     * - Flow stress calculation
     * - Adaptive time stepping
     * - Output visualization
     * 
     * **Template parameters**:
     * @tparam SlipGeom Slip geometry class
     * @tparam SlipKinetics Kinetics class
     * @tparam Elasticty Thermoelastic model class (typically ThermoElastN)
     * 
     * @param[out] pl_disipation_rate Plastic dissipation Ḋ_p [power/volume]
     * @param[out] effective_shear_rate Effective shear rate γ̇_eff [1/time]
     * @param[out] gdot Slip rates γ̇^α [1/time, nslip components]
     * @param[in] inv_det_v_e Inverse elastic volume J_e^(-1) [dimensionless]
     * @param[in] elast_strain Elastic strain tensor ε_e [dimensionless, ntvec]
     * @param[in] kinetic_values Material state for kinetics [various units]
     * @param[in] slip_geom Reference to slip geometry
     * @param[in] slip_kinetics Reference to kinetics
     * @param[in] elasticity Reference to elastic model
     * 
     * **Computation sequence**:
     * 1. Compute Kirchhoff stress from elastic strain
     * 2. Resolve stress onto slip systems
     * 3. Evaluate slip rates from resolved stresses
     * 4. Sum contributions to get effective rate
     * 5. Compute dissipation as stress·rate inner product
     * 
     * **Conditional compilation**:
     * - ECMECH_USE_DPEFF: Use tensor-based effective rate (more accurate)
     * - Otherwise: Use sum of absolute slip rates (faster, approximate)
     * 
     * @note This function is typically called post-convergence for output
     * @note Not used in residual/Jacobian (hence "contributions" not "terms")
     * @note Outputs initialized to zero before adding slip contributions
     * 
     * @see vecd_Deff() for tensor-based effective rate calculation
     * @see vecsssumabs() for sum of absolute values
     */
    template<class SlipGeom, class SlipKinetics, class Elasticty>
    __ecmech_hdev__
    inline
    void get_slip_contributions(double& pl_disipation_rate,
                                double& effective_shear_rate,
                                double* const gdot,
                                const double inv_det_v_e,
                                const double* const elast_strain,
                                const double* const kinetic_values,
                                const SlipGeom& slip_geom,
                                const SlipKinetics& slip_kinetics,
                                const Elasticty& elasticity
                                )
    {
        pl_disipation_rate = 0.0;
        effective_shear_rate = 0.0;
        if constexpr (SlipGeom::nslip > 0) {
        // default initialize everything to 0.0 
        constexpr size_t nslip_dyn = (SlipGeom::dynamic) ? (SlipGeom::nslip + SlipGeom::nSlipExtra) : SlipGeom::nslip;
        double abs_resolved_shear_stress[nslip_dyn] = {};
        double junk[SlipGeom::nslip] = {};
        double kirchoff[ecmech::nsvec] = {};
        elasticity.elast_strain_to_kirchoff_stress(kirchoff, elast_strain);
        // resolve stress onto slip systems
        if constexpr (SlipGeom::dynamic) {
            slip_geom.getExtras(&abs_resolved_shear_stress[SlipGeom::nslip]);
        }
        slip_geom.evalRSS(abs_resolved_shear_stress, kirchoff, slip_geom.getP());
        slip_kinetics.evalGdots(gdot, junk, abs_resolved_shear_stress, kinetic_values);
#if defined(ECMECH_USE_DPEFF)
        double plastic_def_rate_d5[ntvec] = {};
        vecsVMa<ntvec, SlipGeom::nslip>(plastic_def_rate_d5, slip_geom.getP(), gdot);
        effective_shear_rate = vecd_Deff(plastic_def_rate_d5);
#else
        effective_shear_rate = vecsssumabs<SlipGeom::nslip>(gdot);
#endif
        pl_disipation_rate = inv_det_v_e * vecsyadotb<SlipGeom::nslip>(abs_resolved_shear_stress, gdot);
        }
    }

    /**
     * @brief Compute higher-order elastic rotation coupling terms for finite strain.
     * 
     * This function calculates geometric correction terms that arise from the coupling
     * between elastic strain and lattice rotation in finite strain kinematics. These
     * terms are essential for:
     * - Elastic spin contribution to lattice rotation evolution
     * - Jacobian coupling between elastic strain and rotation unknowns
     * - Accurate representation of large elastic deformations
     * 
     * **Physical interpretation**:
     * In finite strain theory, the lattice rotation rate is not simply the difference
     * between total and plastic spin. Additional "elastic spin" terms arise from:
     * - Rotation of principal elastic strain directions
     * - Corotational derivative corrections
     * - Geometric nonlinearity in strain-rotation coupling
     * 
     * **Mathematical background**:
     * The elastic spin contribution can be expressed as:
     * ```
     * W_e = ee_fac * A_e · ε̇_e
     * ```
     * where:
     * - W_e: Elastic spin vector [nwvec = 3]
     * - A_e: Coupling matrix relating elastic strain to spin [nwvec × ntvec = 3×5]
     * - ε̇_e: Elastic strain rate [ntvec = 5]
     * - ee_fac: Scaling factor = 0.5 * (1/a_vol)²
     * 
     * **Computation sequence**:
     * 1. Compute A_e matrix via M35_d_AAoB_dA(A_e_M35, elast_d5)
     *    - This computes ∂(A ⊗ B)/∂A where A = B = elastic strain
     *    - Results in [nwvec × ntvec] matrix
     * 
     * 2. Compute elastic spin: ee_spin_vec = A_e_M35 · elast_dt_d5
     *    - Matrix-vector product: [3×5] · [5] = [3]
     *    - Represents W_e before scaling by ee_fac
     * 
     * 3. Compute scaling factor: ee_fac = 0.5 * inv_a_vol²
     *    - a_vol = J^(1/3) is the volume scaling factor
     *    - inv_a_vol = 1/a_vol
     *    - Factor of 0.5 from corotational derivative formulation
     * 
     * **Output usage**:
     * - **A_e_M35**: Used in Jacobian computation
     *   - Enters ∂R_ω/∂ε_e coupling block
     *   - Time derivative A_edot computed separately from this
     * 
     * - **ee_spin_vec**: Used in rotation residual
     *   - Added to rotation evolution equation (before scaling by ee_fac)
     *   - Typically small for moderate elastic strains
     * 
     * - **ee_fac**: Scaling factor applied to elastic spin terms
     *   - Multiplies ee_spin_vec in residual evaluation
     *   - Appears in Jacobian derivative terms
     *   - Magnitude: O(1) for typical elastic strains
     * 
     * **When these terms matter**:
     * - Large elastic strains (> 1% deviatoric strain)
     * - Stiff materials with significant elastic anisotropy
     * - High strain rate loading
     * - Accurate lattice rotation tracking
     * 
     * **When these terms can be neglected**:
     * - Small strain formulations (< 0.1% elastic strain)
     * - Isotropic elastic response
     * - Quasi-static loading
     * - Set ee_fac = 0 to disable
     * 
     * @param[out] A_e_M35 Elastic strain-spin coupling matrix [nwvec × ntvec = 15 components]
     *                     Layout: row-major, access via ECMECH_NM_INDX(i, j, nwvec, ntvec)
     *                     Physical meaning: ∂W_e/∂ε_e (before ee_fac scaling)
     * 
     * @param[out] ee_spin_vec Elastic spin vector [nwvec = 3 components]
     *                         Axial vector form of elastic spin tensor
     *                         Must be multiplied by ee_fac to get actual W_e
     *                         Units: [1/time] when multiplied by ee_fac
     * 
     * @param[out] ee_fac Elastic-elastic coupling factor [dimensionless]
     *                    Value: 0.5 / a_vol²
     *                    Scaling applied to elastic spin contributions
     *                    Typical magnitude: O(1) for a_vol ≈ 1
     * 
     * @param[in] inv_a_vol Inverse volume scaling J^(-1/3) [dimensionless]
     *                      Where J = det(F_e) is elastic deformation determinant
     *                      a_vol = J^(1/3) relates to elastic volume change
     * 
     * @param[in] elast_d5 Elastic deviatoric strain [ntvec = 5 components, dimensionless]
     *                     Beginning-of-step elastic strain in lattice frame
     *                     Used to compute A_e_M35 coupling matrix
     * 
     * @param[in] elast_dt_d5 Elastic strain rate ε̇_e [1/time, ntvec = 5 components]
     *                        Typically: ε̇_e = (ε_e^{n+1} - ε_e^n) / Δt
     *                        Used to compute ee_spin_vec = A_e · ε̇_e
     * 
     * @see M35_d_AAoB_dA() for A_e_M35 computation details
     * @see EvptnLatticeRotationProblem::get_omega_residual() for residual usage
     * @see EvptnLatticeStrainProblem::get_deriv_omega_wrt_elast_strain() for Jacobian usage
     */
    __ecmech_hdev__
    inline
    void elasticity_higher_order_terms(double* const A_e_M35,
                                        double* const ee_spin_vec,
                                        double& ee_fac,
                                        const double inv_a_vol,
                                        const double* const elast_d5,
                                        const double* const elast_dt_d5
                                    )
    {
        // from e edot product term in spin (formerly neglected)
        M35_d_AAoB_dA(A_e_M35, elast_d5);
        vecsVMa<nwvec, ntvec>(ee_spin_vec, A_e_M35, elast_dt_d5);
        ee_fac = onehalf * inv_a_vol * inv_a_vol;
    }

    /**
     * @brief Compute material tangent stiffness dσ/dD for finite element codes.
     * 
     * This function evaluates the consistent algorithmic tangent operator relating
     * Cauchy stress increments to deformation rate increments in the sample/global frame.
     * The tangent is computed via the implicit function theorem applied to the converged
     * implicit integration equations.
     * 
     * **Physical interpretation**:
     * The material tangent answers: "How does stress σ change if I perturb the
     * prescribed deformation rate D while maintaining equilibrium?"
     * 
     * **Mathematical formulation**:
     * ```
     * dσ/dD = ∂σ/∂ε_e · dε_e/dD + ∂σ/∂ω · dω/dD
     * ```
     * where:
     * - σ: Cauchy stress in sample frame [nsvec = 6]
     * - D: Deformation rate in sample frame [ntvec = 5 deviatoric]
     * - ε_e: Elastic strain in crystal frame [ntvec = 5]
     * - ω: Lattice rotation [nwvec = 3]
     * 
     * **Note**: The output is a 6×6 matrix, but this function only populates the
     * 5×5 deviatoric block (upper-left). The 6th row and column correspond to
     * volumetric/pressure contributions, which are computed separately from the
     * equation of state (EOS) and added in computeTangentStiffness().
     * 
     * **Computation approach** (implicit function theorem):
     * 1. From converged residual R(ε_e, ω, D) = 0, we have:
     *    ```
     *    dR/dε_e · dε_e/dD + dR/dω · dω/dD + dR/dD = 0
     *    ```
     * 
     * 2. Solve for sensitivities:
     *    ```
     *    [dε_e/dD] = -J^(-1) · [dR/dD]
     *    [dω/dD  ]             [  0   ]
     *    ```
     *    where J is the converged Jacobian matrix [JAC_SIZE × JAC_SIZE]
     * 
     * 3. Apply chain rule through elasticity and rotation:
     *    ```
     *    dσ/dD = (∂σ/∂ε_e · dε_e/dD) + (∂σ/∂Q · ∂Q/∂ω · dω/dD)
     *    ```
     * 
     * **Decomposition into contributions**:
     * 
     * **A. Elastic tangent** (∂σ/∂ε_e · dε_e/dD):
     * - dε_e/dD obtained by solving linear system with Jacobian
     * - Apply elastic stiffness C to get dσ/dε_e
     * - Rotate result to sample frame
     * 
     * **B. Rotation tangent** (∂σ/∂Q · ∂Q/∂ω · dω/dD):
     * - dω/dD obtained from same linear solve
     * - ∂Q/∂ω via quaternion exponential map derivative
     * - ∂σ/∂Q from stress rotation derivative
     * - Chain rule composition
     * 
     * **Implementation steps**:
     * 
     * 1. **Set up RHS**: dR/dD = rotation matrix (frame transform derivative)
     *    - nRHS = nsvec = 6: five deviatoric columns plus the spherical
     *      (volumetric) column, whose RHS is the plastic-rate sensitivity to the
     *      spherical elastic strain (nonzero only with deviatoric-volumetric
     *      elastic coupling, e.g. hexagonal symmetry)
     * 2. **Solve linear system**: J · [dε_e/dD; dω/dD] = -[dR/dD; ...]
     *    - Uses SNLS_LUP_SolveX for multiple RHS (nsvec = 6 columns)
     *    - Jacobian already factored during nonlinear solve
     * 3. **Apply elastic stiffness**: temp = C · (dε_e/dD)
     *    - Via ThermoElastN::multCauchyDif
     *    - Produces 6×6 intermediate result including the pressure row and the
     *      crystal-frame direct deviatoric-from-volumetric coupling term
     * 4. **Rotate to sample frame**: result = Q · temp · Q^T
     *    - Via qr6x6_pre_mul
     *    - Extends to 6×6 (with 6th row/column handled appropriately)
     * 5. **Add rotation contribution**: result += ∂(Qσ)/∂ω · (dω/dD)
     *    - Derivative of rotated stress w.r.t. rotation
     *    - Affects rows 0-4 across all 6 columns
     * 
     * **Template parameters**:
     * @tparam ThermoElastN Thermoelastic model class (e.g., ThermoElastNCubic)
     *                      Must provide multCauchyDif() method
     * @tparam JAC_SIZE Total Jacobian dimension from implicit solve
     *                  Typically 8 (ntvec + nwvec) for coupled elastic-rotation problem
     * @tparam ind_sub_omega Starting index for rotation DOFs in Jacobian
     *                       Typically ntvec = 5 (rotation DOFs follow elastic strain DOFs)
     * 
     * @param[out] material_tangent Material tangent stiffness [nsvec × nsvec = 6×6 = 36 components]
     *                              Row-major storage: element (i,j) at index i*nsvec + j
     *                              Units: [stress/strain_rate] or [stress·time]
     *                              **Note**: All entries populated except the EOS part of the
     *                              (S,S) pressure-volume stiffness, which the caller adds
     *                              (computeTangentStiffness). For cubic symmetry the 6th row
     *                              and column are zero, as before.
     * 
     * @param[in] jacobian Converged Jacobian matrix from implicit solve [JAC_SIZE × JAC_SIZE]
     *                     Must be the Jacobian at the converged solution
     *                     Will be used for linear solve (non-const for LU factorization)
     * 
     * @param[in] dquat_domega_t Transposed quaternion derivative [nwvec × qdim = 3×4]
     *                           ∂Q/∂ω evaluated at converged rotation
     *                           From exponential map: ∂(exp(ω))/∂ω
     * 
     * @param[in] rmat_5x5_sample2xtal Rotation matrix crystal→sample [ntvec × ntvec = 5×5]
     *                                  For rotating deviatoric tensors
     *                                  Transpose of crystal→sample rotation for 5-vectors
     * 
     * @param[in] quat Converged orientation quaternion [qdim = 4]
     *                 Represents sample→crystal frame rotation
     *                 Used for rotation derivative computations
     * 
     * @param[in] rmat Converged rotation matrix [ndim × ndim = 3×3]
     *                 Standard rotation matrix from quat
     *                 Used for vector/tensor rotations
     * 
     * @param[in] cauchy_stress Converged Cauchy stress in crystal frame [nsvec = 6]
     *                          Used in rotation derivative ∂(Qσ)/∂Q
     *
     * @param[in] dDp_hat_deps_sph Derivative ∂D_p/∂ε_s [ntvec], from
     *                             get_slip_rate_deriv_sph_terms()
     *
     * @param[in] dWp_hat_deps_sph Derivative ∂W_p/∂ε_s [nwvec], likewise
     *
     * @param[in] elast_dev_press_vec Scaled strain vector [nsvec] (as eval() takes),
     *                                for the direct J-scaling column terms
     *
     * @param[in] elast_dt_d5 Elastic strain rate [ntvec], for the 1/a_vol
     *                        strain-rate sensitivity in the volumetric RHS
     *
     * @param[in] dt Time step Δt; sets ∂ε_s/∂D_s = Δt for the volumetric column
     *
     * @param[in] inv_det_v_e Inverse elastic deformation determinant J_e^(-1) [dimensionless]
     *                        Scaling factor for Kirchhoff→Cauchy conversion
     * 
     * @param[in] inv_a_vol Inverse volume scaling a_vol^(-1) = J_e^(-1/3) [dimensionless]
     *                      Normalization for deviatoric elastic response
     * 
     * @param[in] thermo_elast_n Thermoelastic model reference
     *                           Used for multCauchyDif() elastic tangent operation
     * 
     * **Accuracy considerations**:
     * - Tangent is "consistent" with implicit integration (quadratic convergence in Newton)
     * - Assumes converged state (residual ≈ 0)
     * - Neglects EOS tangent contributions (handled separately in computeTangentStiffness)
     * - Valid for small perturbations around converged state
     * 
     * **Alternative approaches**:
     * - Finite difference: dσ/dD ≈ (σ(D+δD) - σ(D)) / δD
     *   Pros: Simple, no derivative computation
     *   Cons: Expensive (requires full nonlinear solve for each column), numerical errors
     * 
     * - This method (implicit function theorem):
     *   Pros: One linear solve for all columns, exact to machine precision, consistent
     *   Cons: Requires Jacobian, more complex implementation
     * 
     * **Limitations**:
     * - Tangent w.r.t. deformation rate D, not strain increment ε
     * - Scaling by time step Δt needed for dσ/dε (done in computeTangentStiffness wrapper)
     * - EOS part of the (S,S) pressure-volume stiffness added separately by the caller
     * - Volumetric column neglects derivatives of the J-dependent scaling factors
     *   (inv_det_v_e, inv_a_vol) and of pressure_EOS/Grüneisen terms, O(σ/B) relative
     * - Assumes small perturbations (linear approximation)
     * 
     * @warning Jacobian must be from converged solution for accurate tangent
     * @warning Assumes Jacobian is non-singular (converged state should guarantee this)
     * 
     * @see computeTangentStiffness() in ECMech_evptnSngl.h for full tangent assembly
     * @see SNLS_LUP_SolveX() for linear solver details
     * @see ThermoElastN::multCauchyDif() for elastic tangent computation
     * @see qr6x6_pre_mul() for frame rotation operations
     * @see eval_d_dxi_impl_quat() for rotation derivative computations
     */
    template<class ThermoElastN, size_t JAC_SIZE, size_t ind_sub_omega>
    __ecmech_hdev__
    inline
    void get_material_tangent_stiffness(double* const material_tangent,
                                        const double* const jacobian,
                                        const double* const dquat_domega_t,
                                        const double* const rmat_5x5_sample2xtal,
                                        const double* const quat,
                                        const double* const rmat,
                                        const double* const cauchy_stress,
                                        const double* const dDp_hat_deps_sph,
                                        const double* const dWp_hat_deps_sph,
                                        const double* const elast_dev_press_vec,
                                        const double* const elast_dt_d5,
                                        const double dt,
                                        const double inv_det_v_e,
                                        const double inv_a_vol,
                                        const ThermoElastN& thermo_elast_n
                                        )
    {
        // dCauchy/dDefRate over all nsvec vecds components of the deformation
        // rate: five deviatoric columns plus the spherical (volumetric) column.
        // The spherical column is nonzero whenever the elasticity model couples
        // deviatoric stress to volumetric strain (hexagonal K_sdax3): it carries
        // both the direct elastic term and the plastic feedback obtained from the
        // extra right-hand side below.
        constexpr int nRHS = ecmech::nsvec;

        // mtan_sample_frame
        // = d/dDefRate_sample(Cauchy) = d/dDefRate(Q * alpha * (C_{elas} : lattice_strain))
        // Q is the 5x5 rotation oper from crystal to sample
        // alpha is necessary scaling from Kirchoff to Cauchy
        // C_{elas} is the elasticity tensor (deviatoric contributions)
        // = d/dDefRate_s (Cauchy) = d(Q Cauchy) / dOmega_c * dOmega / dDefRate_s
        //   + Q * d(Cauchy)/delast_strain * delast_strain / dDefRate_s
        // From our Jacobian we can calculate the dOmega / dDefRate_s and delast_strain / dDefRate_s terms
        // The other ones are either simple to calculate or require some math...
        //
        // dstrainomega_ddef_rate_t => [dlat_strain_ddef_rate_sample; domega_ddef_rate_sample];
        // Initially we set it to be our RHS
        double dstrainomega_ddef_rate_t[ nRHS * JAC_SIZE ] = {}; // transpose for use in SNLS_LUP_SolveX !
        {
        // RHS calculations
        //
        // dstrainomega_ddef_rate_t => dResidual_ddef_rate_sample
        {
            // negatives cancel
            // dstrainomega_ddef_rate_t[0:ind_omega_vec,:] = qr5x5_c2s
            for (int jE = 0; jE < ecmech::ntvec; ++jE) {
                for (int iE = 0; iE < ntvec; ++iE) { // ntvec, _not_ nDimSys // iE is same as index
                    dstrainomega_ddef_rate_t[ECMECH_NM_INDX(jE, iE, nRHS, JAC_SIZE)] = rmat_5x5_sample2xtal[ECMECH_NN_INDX(jE, iE, ntvec)];
                }
            }
            // Spherical (volumetric) RHS: the residuals depend on the spherical
            // deformation rate through the spherical elastic strain,
            // d(eps_sph)/d(D_s) = dt, along two paths: the stress sensitivity
            // dT/deps_sph feeding the slip kinetics (dDp/dWp terms), and the
            // 1/a_vol = exp(-eps_sph/√3) scaling on the elastic strain-rate term
            // of the strain residual,
            //   d(inv_a_vol · elast_dt)/deps_sph = -(1/√3) inv_a_vol elast_dt .
            // In the unscaled Jacobian's units:
            //   strain rows: dR_eps/dD_s = dt * ( dDp/deps_sph
            //                                     - (1/√3) inv_a_vol elast_dt )
            //   omega rows:  dR_omega/dD_s = dt * ( dt * dWp/deps_sph )
            // and the tangent solve needs -dR/dD_s (no sign cancellation here,
            // unlike the deviatoric columns).
            for (int iE = 0; iE < ecmech::ntvec; ++iE) {
                dstrainomega_ddef_rate_t[ECMECH_NM_INDX(iSvecS, iE, nRHS, JAC_SIZE)] =
                    -dt * (dDp_hat_deps_sph[iE] - ecmech::sqr3i * inv_a_vol * elast_dt_d5[iE]);
            }
            for (int iW = 0; iW < ecmech::nwvec; ++iW) {
                dstrainomega_ddef_rate_t[ECMECH_NM_INDX(iSvecS, ind_sub_omega + iW, nRHS, JAC_SIZE)] = -dt * dt * dWp_hat_deps_sph[iW];
            }
            // If we had hardening terms then we'd add those here as well...
            // dRhard_ddef_rate but typically we can but for most problems safe to treat that
            // set of terms as being = 0
        }
        // Now solve for our dstrain_ddefrate and domega_ddefrate terms
        int err = SNLS_LUP_SolveX<JAC_SIZE>(const_cast<double*>(jacobian), dstrainomega_ddef_rate_t, nRHS);
        if (err != 0) {
            ECMECH_FAIL(__func__, "error from SNLS_LUP_SolveX");
        }
        }

        {
        // Gather the strain sensitivities into [ntvec x nsvec] layout,
        // dstrain_ddef_rate[i][j] = d(elast_d5)_i / d(def_rate_vecds)_j, leaving
        // the solution array untouched (its omega entries are read below).
        double dstrain_ddef_rate[ ecmech::ntvec * ecmech::nsvec ];
        for (int iTvec = 0; iTvec<ecmech::ntvec; ++iTvec) {
            for (int jSvec = 0; jSvec<ecmech::nsvec; ++jSvec) {
                dstrain_ddef_rate[ECMECH_NM_INDX(iTvec, jSvec, ecmech::ntvec, ecmech::nsvec)] =
                    dstrainomega_ddef_rate_t[ECMECH_NM_INDX(jSvec, iTvec, nRHS, JAC_SIZE)];
            }
        }
        double temp_M6[ ecmech::nsvec2 ];
        thermo_elast_n.template multCauchyDif<ecmech::ntvec, ecmech::nsvec>(temp_M6, dstrain_ddef_rate, elast_dev_press_vec, dt, inv_det_v_e, inv_a_vol);
        // Apply final rotation; the crystal-frame direct coupling term placed by
        // multCauchyDif in the spherical column rotates with the other rows here
        qr6x6_pre_mul<ecmech::nsvec, false>(material_tangent, temp_M6, rmat_5x5_sample2xtal);
        }

        // Calculate the d(QCauchy) / dDefRate_s term now and add that to
        {
        double dcauchy_dquat[ ecmech::ntvec * ecmech::qdim ];
        {
            double drmat_dquat[ ecmech::ndim * ecmech::ndim * ecmech::qdim ];
            d_quat_to_tensor(drmat_dquat, quat);
            double dcauchy_drmat[ ecmech::ntvec * ecmech::ndim * ecmech::ndim ];
            d_rot_mat_vecd_smop(dcauchy_drmat, rmat, cauchy_stress);
            vecsMAB<ntvec, qdim, ndim*ndim>(dcauchy_dquat, dcauchy_drmat, drmat_dquat);
        }

        // We now need to be able to go from our domega_ddef_rate to dquat_ddef_rate
        // dquat_ddef_rate = dquat_domega_t * domega_ddef_rate_t
        double dquat_ddef_rate[ ecmech::qdim * nRHS ];
        for (int ii_I = 0; ii_I < nRHS; ++ii_I) {
            for (int ii_Q = 0; ii_Q < ecmech::qdim; ++ii_Q) {
                int iiQI = ECMECH_NM_INDX(ii_Q, ii_I, ecmech::qdim, nRHS);
                dquat_ddef_rate[iiQI] = 0.0;
                for (int ii_W = 0; ii_W < ecmech::nwvec; ++ii_W) {
                    dquat_ddef_rate[iiQI] +=
                    dquat_domega_t[ECMECH_NM_INDX(ii_W, ii_Q, ecmech::nwvec, ecmech::qdim)] *
                    dstrainomega_ddef_rate_t[ECMECH_NM_INDX(ii_I, ind_sub_omega + ii_W, nRHS, JAC_SIZE)];
                }
            }
        }

        // Now get the dcauchy_lattice_dI terms by doing ->
        // dcauchy_lattice_dI = d_cauchy_lattice_dquat * dRmat_quat_dI
        double dqcauchy_ddefrate[ ecmech::ntvec * nRHS ];
        vecsMAB<ecmech::ntvec, nRHS, ecmech::qdim>(dqcauchy_ddefrate, dcauchy_dquat, dquat_ddef_rate);
        for (int ii_T = 0; ii_T < ecmech::ntvec; ++ii_T) {
            for (int ii_I = 0; ii_I < nRHS; ++ii_I) {
                // NOTE : only looping over ntvec, but mtan_sI is nsvec in the first dimension
                material_tangent[ECMECH_NN_INDX(ii_T, ii_I, ecmech::nsvec)] += dqcauchy_ddefrate[ECMECH_NM_INDX(ii_T, ii_I, ecmech::ntvec, nRHS)];
            }
        }
        }
    }
}
}