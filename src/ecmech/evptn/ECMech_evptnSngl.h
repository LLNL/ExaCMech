/**
 * @file ECMech_evptnSngl.h
 * @brief Single material point integration functions for crystal plasticity.
 * 
 * This file provides the high-level interface for crystal plasticity time integration
 * at a single material point (Gauss point). It orchestrates:
 * - State initialization and preparation (preprocess)
 * - Nonlinear system solution (main_problem)
 * - State extraction and post-processing (postprocess)
 * - Optional material tangent computation (computeTangentStiffness)
 * 
 * **Main entry point**:
 * - `getResponseSngl()`: Complete time step integration for one material point
 * - `getResponseNRSngl()`: Alternative with simplified rotation treatment (optional)
 * 
 * **Integration workflow**:
 * ```
 * 1. preprocess():
 *    - EOS update (pressure, temperature, bulk modulus)
 *    - Hardening state update (may involve sub-iterations)
 *    - Compute initial stress and slip rates
 *    - Accumulate deviatoric strain energy (trapezoidal rule)
 * 
 * 2. main_problem():
 *    - Set up SNLS trust region solver
 *    - Initialize guess (typically zeros)
 *    - Solve nonlinear system via Newton iterations
 *    - Check convergence
 * 
 * 3. postprocess_prob():
 *    - Extract converged elastic strain and orientation
 *    - Update hardening state with converged slip rates
 *    - Compute stress in crystal frame
 *    - Calculate plastic dissipation and flow strength
 * 
 * 4. postprocess():
 *    - Rotate stress back to sample frame
 *    - Update internal energy with plastic work
 *    - Finalize strain energy integration
 *    - Compute material properties (bulk, shear modulus)
 *    - Store final stress state
 * 
 * 5. computeTangentStiffness() [optional]:
 *    - Compute material tangent via implicit function theorem
 *    - Scale from rate to increment form
 *    - Add EOS volumetric contribution
 *    - Convert from deviatoric to full Voigt notation
 * ```
 * 
 * **Function organization**:
 * - **Template functions**: Templated by constitutive model types for flexibility
 * - **Inline functions**: Header-only for GPU compatibility and optimization
 * - **__ecmech_hdev__**: All functions work on both CPU and GPU
 * 
 * **Error handling**:
 * - Return bool: true = success, false = convergence failure
 * - Failed solves return false, allowing host code to handle errors
 * - Debug builds may print diagnostic info on failure
 * 
 * **Solver configuration**:
 * - Uses SNLS trust region dense Gauss-Newton solver
 * - Default tolerance: passed as parameter (typically 1e-8 to 1e-12)
 * - Max iterations: 200 (hard-coded in main_problem)
 * - Initial trust region: 1.0
 * 
 * **Advanced features**:
 * - Material tangent computation for consistent FEM tangent operator
 * - Multiple solver variants (full coupled, strain-only, rotation-only)
 * - RStarSolve: Optional secondary rotation solve for improved accuracy
 * 
 * @see ECMech_evptn.h for problem formulation classes
 * @see ECMech_base_classes.h for state container (ProblemState)
 * @see SNLS_TrDLDenseG.h for nonlinear solver details
 */

#pragma once

#include "ECMech_core.h"
#include "ECMech_util.h"
#include "ECMech_elastic.h"
#include "ECMech_eosSimple.h"
#include "evptn/ECMech_base_classes.h"
#include "evptn/ECMech_base_fcns.h"
#include "evptn/ECMech_evptn.h"

#include "SNLS_TrDLDenseG.h"
#include "SNLS_HybrdTrDLDenseG.h"

namespace ecmech {
namespace evptn {

/**
 * @brief Prepare state for crystal plasticity time integration.
 * 
 * This function performs all computations needed before the main nonlinear solve:
 * 1. **Thermodynamic update**: EOS evaluation for pressure, temperature, bulk modulus
 * 2. **Energy integration**: Start deviatoric strain energy accumulation (trapezoidal rule)
 * 3. **Hardening update**: Evolve hardening state using beginning-of-step slip rates
 * 
 * **Thermodynamic calculations**:
 * - Evaluate temperature at beginning of step from EOS
 * - Update pressure and temperature to end of step (simple EOS model)
 * - Compute bulk modulus at new state
 * - Store thermodynamic derivatives (dp/de, dp/dv, dT/de)
 * 
 * **Energy integration**:
 * - Uses trapezoidal rule for deviatoric part: ∫σ:D dt ≈ (σ_n + σ_{n+1})/2 * D * Δt
 * - Begins accumulation with n-state contribution
 * - Final contribution added in postprocess after stress converges
 * 
 * **Hardening update**:
 * - Uses beginning-of-step slip rates (from history or initial guess)
 * - Integrates hardening ODEs over time step
 * - May involve inner Newton iterations for implicit hardening laws
 * - Updates prob_state.h_state_u with new hardening state
 * - For dynamic slip systems, evaluates extra quantities (e.g., chi angles)
 * 
 * **Template parameters**:
 * @tparam SlipGeom Slip geometry class
 * @tparam Kinetics Kinetics class
 * @tparam EosModel Equation of state class
 * @tparam ThermoElastN Thermoelastic class
 * @tparam ProbState Problem state container
 * @tparam RStarSolve If true, enables secondary rotation solve (default: false)
 * 
 * @param[in] slipGeom Slip geometry object
 * @param[in] kinetics Kinetics object
 * @param[in] eos Equation of state object
 * @param[in] thermoElastN Thermoelastic object
 * @param[in] rel_vol_ratios Volume ratio array [nvr]:
 *                           [0] = J_n, [1] = J_{n+1}, [2] = (J_{n+1}-J_n)/dt, [3] = J_{n+1}-J_n
 * @param[in] internal_energy Energy array [ne]: [0] = total internal energy
 * @param[in] def_rate_d6v_sample Deformation rate [nsvp], uses first ntvec components
 * @param[in,out] prob_state State container, modified with updated thermodynamics and hardening
 * @param[out] halfVMidDt Half volume × midpoint time: 0.25*(J_n + J_{n+1})*Δt
 * @param[out] dev_strain_energy_total Accumulated deviatoric strain energy [energy/volume]
 * 
 * @return true if preprocessing successful, false if hardening update failed to converge
 * 
 * **State modifications** (via prob_state):
 * - tkelv: Updated to end-of-step temperature
 * - pressure_EOS: Updated to end-of-step pressure
 * - energy_new: End-of-step internal energy estimate
 * - bulk_modulus_new: Bulk modulus at new state
 * - h_state_u: Updated hardening state variables
 * 
 * **Failure modes**:
 * - Hardening solver fails to converge → return false
 * - Currently only hardening can fail; EOS and energy updates assumed robust
 * 
 * @note Called once per time step before nonlinear solve
 * @note Side effects on prob_state are essential (not const)
 * @note Energy contribution completed in postprocess() after stress converges
 */
template<class SlipGeom, class Kinetics, class EosModel, class ThermoElastN, class ProbState, bool RStarSolve=false>
__ecmech_hdev__
inline
bool preprocess(const SlipGeom& slipGeom,
                const Kinetics& kinetics,
                const EosModel& eos,
                const ThermoElastN& thermoElastN,
                const double* const rel_vol_ratios,
                const double* const internal_energy,
                const double* const def_rate_d6v_sample,
                ProbState& prob_state,
                double& halfVMidDt,
                double& dev_strain_energy_total)
{
    // total increment in the deviatoric part of the strain energy using
    // trapezoidal rule integration
    //
    // just beginning-of-step stress part so far
    //
    halfVMidDt = oneqrtr * (rel_vol_ratios[0] + rel_vol_ratios[1]) * prob_state.dt;
    dev_strain_energy_total = halfVMidDt * vecsInnerSvecDev(prob_state.cauchy_stress_d6p, def_rate_d6v_sample);

    // EOS
    //
    const double energy_old = internal_energy[ecmech::i_ne_total];
    //
    // get tkelv from beginning-of-step to avoid tangent stiffness contributions
    {
        double pressure_BOS;
        const double rel_vol_old = rel_vol_ratios[0];
        eos.evalPT(pressure_BOS, prob_state.tkelv, rel_vol_old, energy_old);
    }

    {
        const double pressure_old = prob_state.cauchy_stress_d6p[6];
        double tkelv_new, dpde, dpdv, dtde;
        updateSimple(eos, prob_state.pressure_EOS, tkelv_new, prob_state.energy_new, prob_state.bulk_modulus_new,
                     dpde, dpdv, dtde,
                     rel_vol_ratios[1], rel_vol_ratios[3],
                     energy_old, pressure_old);
    }

    // update hardness state to the end of the step
    // gdot is still at beginning-of-step
    //
    constexpr size_t nslip_dyn = (SlipGeom::dynamic) ? (SlipGeom::nslip) : 1;
    double hvals[nslip_dyn] = {}; // additional values needed to update the hardening state
    if constexpr(SlipGeom::dynamic) {

        double kirchoff[ecmech::nsvec] = {};
        double elast_d5v[ecmech::nsvec] = {};
        double stress_dev6_press[ecmech::nsvec+1] = {};

       double a_vol = pow(prob_state.rel_vol_new, onethird);
       double inv_a_vol = 1.0 / a_vol;

        vecsVxa<ntvec>(elast_d5v, inv_a_vol, prob_state.elast_d5_n);
        elast_d5v[iSvecS] = sqr3 * log(a_vol);
        thermoElastN.eval(kirchoff, elast_d5v, prob_state.tkelv, prob_state.pressure_EOS, prob_state.energy_new);
        vecdsToSvecP(stress_dev6_press, kirchoff);

        // For dynamic slip systems we need the chi angle
        double P[ecmech::ntvec * SlipGeom::nslip];
        double Q[ecmech::nwvec * SlipGeom::nslip];
        // still need to rotate stress state back to original value
        slipGeom.getPQ(hvals, P, Q, stress_dev6_press);
    }
    const int nfevals = kinetics.updateH(prob_state.h_state_u, prob_state.h_state, prob_state.dt, prob_state.gdot, hvals, prob_state.tkelv);
    if (nfevals < 0) {
        ECMECH_WARN(__func__, "Hardening failed to converge");
        return false;
    }
#if defined(ECMECH_EXTRA_SOLVERS)
    if constexpr (RStarSolve) {
        auto prob = RotUpdProblem(slipGeom, thermoElastN, prob_state);
         // update Rstar aka Rdot * dt
         // gdot is still at beginning-of-step
        snls::SNLSTrDlDenseG<decltype(prob)> solver(prob);
        const bool status = main_problem(1e-8, solver, 0);
        if (!status) {
            ECMECH_WARN(__func__, "RStar solver failed to converge");
            return false;
        }
        prob.stateFromX(prob_state.quat_u, solver._x);
    }
#endif
    return true;
}

/**
 * @brief Solve the nonlinear crystal plasticity system via trust region Newton.
 * 
 * This function configures and executes the SNLS SNLSTrDlDenseG solver to find
 * the solution of the implicit time integration equations. It handles:
 * - Solver initialization and configuration
 * - Initial guess setup
 * - Solution iteration
 * - Convergence checking
 * 
 * **Solver configuration**:
 * - Algorithm: Dogleg approximation to the trust-region sub-problem for multi-dimensional nonlinear systems of equations 
 * - Max iterations: 200
 * - Initial trust region: 1.0
 * - Tolerance: Passed as parameter (user-specified)
 * - Delta control: Default SNLS delta control region parameters
 * 
 * **Initial guess**:
 * - All unknowns initialized to zero
 * - Assumes incremental formulation (Δε_e = 0, Δω = 0 is reasonable guess)
 * - Solver typically converges in 3-8 iterations from this guess
 * 
 * **Convergence criteria**:
 * - Residual norm ‖R‖ < tolerance
 * - Trust region successfully updated
 * - Status returned by solver: converged or better
 * 
 * **Template parameter**:
 * @tparam SNLS_Solver Solver type (typically snls::SNLSTrDlDenseG<Problem>)
 * 
 * @param[in] tolerance Convergence tolerance [dimensionless]
 *                      Typical values: 1e-8 (standard), 1e-12 (tight), 1e-6 (loose)
 * @param[in,out] solver SNLS solver instance, modified during iterations
 *                       Contains problem reference and solution on exit
 * @param[in] outputLevel Diagnostic output verbosity:
 *                        0 = silent, 1 = iteration info, 2 = detailed debug
 * 
 * @return true if solver converged, false if failed to converge
 * 
 * **Convergence failure handling**:
 * - Returns false if status < converged
 * - Host code can attempt recovery (reduce time step, change strategy, etc.)
 * - Debug builds print residual norm and status code
 * 
 * **Performance notes**:
 * - Analytical Jacobian enables quadratic convergence
 * - Most time steps converge in 2-10 iterations
 * - Ill-conditioned cases may hit iteration limit (200)
 * 
 * @see SNLS documentation for trust region algorithm details
 * @see EvptnUpdstProblem::computeRJ() for residual/Jacobian evaluation
 */
template<class SNLS_Solver>
__ecmech_hdev__
inline
bool main_problem(const double tolerance,
                  SNLS_Solver& solver,
                  const int outputLevel)
{
    snls::TrDeltaControl deltaControl;
    deltaControl._deltaInit = 1e0;
    {
        static constexpr int maxIter = 200;
        solver.setupSolver(maxIter, tolerance, &deltaControl, outputLevel);
    }

    // set initial guess
    //
    for (int iX = 0; iX < solver.getNDim(); ++iX) {
        solver._x[iX] = 0e0;
    }

    snls::SNLSStatus_t status = solver.solve( );
    if (status < snls::converged ) {
#if defined(__ecmech_host_only__)
        std::cout << "trust region solver residual " << solver.getRes() << " exit status " << status << std::endl;
        ECMECH_WARN(__func__, "Solver(s) failed to converge -- will try again with implicit elastic strain solve only");
#endif
        return false;
    }
    return true;
}

/**
 * @brief Compute material tangent stiffness matrix for finite element codes.
 * 
 * This function calculates the consistent tangent operator ∂σ/∂ε needed by
 * implicit finite element codes for Newton-Raphson equilibrium iterations.
 * 
 * **Mathematical formulation**:
 * The tangent is obtained via the implicit function theorem:
 * ```
 * ∂σ/∂ε = ∂σ/∂ε_e * (∂ε_e/∂ε)
 * ```
 * where ∂ε_e/∂ε is obtained from the implicit solve by perturbing RHS.
 * 
 * **Computation approach**:
 * 1. Enable material tangent flag in problem: prob.provideMTan(mtanSD_vecds)
 * 2. Re-evaluate Jacobian: This triggers tangent accumulation
 * 3. Disable tangent flag: prob.clearMTan()
 * 4. Scale from rate to increment: multiply by Δt
 * 5. Add volumetric EOS contribution (3K on diagonal)
 * 6. Convert from deviatoric to full Voigt notation
 * 
 * **Deviatoric vs. full tangent**:
 * - Internal computation in deviatoric-pressure form (6-vector + 1)
 * - EOS provides volumetric stiffness (assumed decoupled)
 * - Final output in standard Voigt 6×6 symmetric form
 * 
 * **Template parameters**:
 * @tparam Problem Problem class (e.g., EvptnUpdstProblem)
 * @tparam Solver Solver class (e.g., snls::SNLSTrDlDenseG)
 * @tparam ProblemState State container class
 * 
 * @param[in,out] prob Problem instance, used for tangent accumulation
 * @param[in,out] solver Solver instance, Jacobian re-evaluated
 * @param[in] prob_state State container for current material state
 * @param[out] mtanSD Material tangent [nsvec2 = 36 components]
 *                    Stored in row-major Voigt notation:
 *                    dσ_11/dε_11, dσ_11/dε_22, ..., dσ_12/dε_12
 * 
 * **Coordinate frame**:
 * - Tangent returned in sample/global frame
 * - Internal computations in crystal frame, rotated back
 * 
 * **Typical values**:
 * - Diagonal: ~Elastic moduli (100-400 GPa for metals)
 * - Off-diagonal: Poisson coupling + plastic flow direction effects
 * - Symmetric: Material response is rate-independent in final form
 * 
 * **Limitations**:
 * - Assumes deviatoric-volumetric decoupling (crude for pressure dependence)
 * - Neglects thermal expansion effects on stiffness
 * - Symmetric tangent (assumes no rate effects in final form)
 * 
 * @note Tangent is with respect to strain INCREMENT Δε, not rate
 * @note EOS contribution simplified: K on (S,S) component only
 * @note Matrix symmetric by construction (no explicit symmetrization needed)
 * 
 * @see mtan_conv_sd_svec() for deviatoric-to-Voigt conversion
 * @see Problem::provideMTan() for tangent accumulation mechanism
 */
template<class Problem, class Solver, class ProblemState>
__ecmech_hdev__
inline
void computeTangentStiffness(Problem& prob,
                             Solver& solver,
                             ProblemState& prob_state,
                             double* const mtanSD)
{
    double mtanSD_vecds[ ecmech::nsvec2 ] = {};
    prob.provideMTan(mtanSD_vecds);
    {
        double residual[Problem::nDimSys] = {};
        double Jacobian[Problem::nDimSys * Problem::nDimSys] = {};
        solver.computeRJ(&residual[0], &Jacobian[0]);
    }
    prob.clearMTan();
    // currently have derivative with-respsect-to deformation rate;
    // to get derivative with-respsect-to strain increment,
    // multiply by 1/dt
    //
    double dt_ri = prob.getDtRi();
    for (int i = 0; i < ecmech::nsvec2; ++i) {
        mtanSD_vecds[i] = mtanSD_vecds[i] * dt_ri;
    }

    // contribution to stiffness from EOS
    // this is a bit crude, but should do the trick for now;
    // neglects effect of pressure_EOS and rel_vol_new on workings of evptn
    //
    mtanSD_vecds[ECMECH_NN_INDX(iSvecS, iSvecS, ecmech::nsvec)] = three * prob_state.bulk_modulus_new;

    // convert from vecds notation to svec notation
    //
    mtan_conv_sd_svec<true>(mtanSD, mtanSD_vecds);
}

/**
 * @brief Extract solution and update state after nonlinear solve convergence.
 * 
 * This function is called immediately after the main nonlinear solve succeeds.
 * It extracts physical state from the solution vector and computes derived
 * quantities needed for output and the next time step.
 * 
 * **Operations performed**:
 * 1. Copy updated hardening state to primary storage
 * 2. Compute stress in crystal frame from converged elastic strain
 * 3. Calculate plastic dissipation rate and effective shear rate
 * 4. Update accumulated equivalent plastic strain
 * 5. Compute flow strength for current state
 * 
 * **Template parameters**:
 * @tparam kinNH Number of hardening state variables (Kinetics::nH)
 * @tparam Problem Problem class
 * @tparam ProblemState State container class
 * 
 * @param[in,out] prob Problem instance, provides slip contribution calculation
 * @param[in,out] prob_state State container:
 *                           IN: elast_d5_u (from stateFromX)
 *                           OUT: h_state, eps, eps_dot, flow_strength
 * @param[out] cauchy_stress_d5p_xtal Cauchy stress in crystal frame [nsvec]
 *                                     Deviatoric components + pressure
 * 
 * **Flow strength calculation**:
 * ```
 * if (D_eff > tiny):
 *     flow_strength = plastic_dissipation_rate / D_eff
 * else:
 *     flow_strength = hardening_scale_factor
 * ```
 * This gives a generalized "yield stress" measure.
 * 
 * **Equivalent plastic strain**:
 * - Integrated as: ε_p^{n+1} = ε_p^n + γ̇_eff * Δt
 * - Provides scalar measure of accumulated plastic deformation
 * - Used for hardening laws, damage models, output
 * 
 * **State updates** (via prob_state):
 * - h_state: Hardening variables copied from h_state_u
 * - eps_dot: Effective plastic strain rate at end of step
 * - eps: Accumulated equivalent plastic strain
 * - flow_strength: Current resistance to plastic flow
 * 
 * @note Called before postprocess() which handles frame transformations
 * @note Stress returned in CRYSTAL frame; rotated to sample frame later
 * 
 * @see prob.get_slip_contribution() for plastic work calculations
 * @see postprocess() for final stress rotation and energy update
 */
template<int kinNH, class Problem, class ProblemState>
__ecmech_hdev__
inline
void postprocess_prob(Problem& prob,
                      ProblemState& prob_state,
                      double* const cauchy_stress_d5p_xtal
                     )
{
    for (int i_hstate = 0; i_hstate < kinNH; i_hstate++) {
        prob_state.h_state[i_hstate] = prob_state.h_state_u[i_hstate];
    }

    double pl_disipation_rate = 0.0;
    double effective_shear_rate = 0.0;

    prob.get_slip_contribution(pl_disipation_rate, effective_shear_rate,
                               prob_state.gdot, prob_state.elast_d5_u);

    prob_state.eps_dot = effective_shear_rate;
    prob_state.eps += prob_state.eps_dot * prob_state.dt;
    //
    {
        double dEff = vecd_Deff(prob_state.def_rate_d5_sample);
        double flow_strength = prob.getHdnScale();
        if (dEff > idp_tiny_sqrt) {
            flow_strength = pl_disipation_rate / dEff;
        }
        prob_state.flow_strength = flow_strength;
    }
        // get Cauchy stress
        //
        prob.elastNEtoC(cauchy_stress_d5p_xtal, prob_state.elast_d5_u);
}

/**
 * @brief Finalize state updates after convergence and post-processing.
 * 
 * This function completes the time step integration by:
 * 1. Rotating stress from crystal to sample frame
 * 2. Finalizing strain energy integration
 * 3. Updating internal energy with plastic work
 * 4. Computing material properties (bulk and shear moduli)
 * 5. Storing final stress state
 * 
 * **Frame transformation**:
 * - Input: Cauchy stress in crystal frame (from postprocess_prob)
 * - Rotation matrix: From updated crystal orientation
 * - Output: Cauchy stress in sample/global frame
 * - Includes both deviatoric and pressure components
 * 
 * **Energy integration completion**:
 * - Add end-of-step contribution to trapezoidal rule
 * - Total: E_dev = (V_{n+1/2} * Δt / 2) * (σ_n : D + σ_{n+1} : D)
 * - Update total internal energy: E += E_dev
 * 
 * **Material properties**:
 * - Bulk modulus: From EOS at current state
 * - Shear modulus: From elasticity at current state  
 * - Stored in sdd array for output
 * 
 * **Quaternion flip correction**:
 * - Check if updated quaternion closer to -Q_n than Q_n
 * - If so, flip sign (antipodal symmetry: Q ≡ -Q)
 * - Keeps orientation clustered, improves post-processing
 * 
 * **Template parameters**:
 * @tparam ProblemState State container class
 * @tparam ThermoElastN Thermoelastic model class
 * 
 * @param[in,out] prob_state State container:
 *                           IN: quat_u, cauchy_stress_d6p (sample frame stress storage)
 *                           OUT: quat_u (possibly flipped), energy_new
 * @param[in] elastN Thermoelastic model for property evaluation
 * @param[in] def_rate_d6v_sample Deformation rate for energy integration
 * @param[out] sdd Auxiliary output array [nsdd = 2]:
 *                 [0] = bulk modulus, [1] = shear modulus
 * @param[in,out] internal_energy Energy array [ne = 1]:
 *                                [0] = total internal energy (updated)
 * @param[in] cauchy_stress_d5p_xtal Stress in crystal frame [nsvec]
 * @param[in] dev_strain_energy_total Accumulated deviatoric strain energy
 * @param[in] halfVMidDt Half volume × midpoint time factor
 * 
 * **Final stress storage**:
 * - Converted from deviatoric-pressure to standard Voigt form
 * - Stored in prob_state.cauchy_stress_d6p for return to host code
 * - Layout: [σ11, σ22, σ33, σ23, σ13, σ12, -p]
 * 
 * **Energy balance**:
 * - Internal energy updated to include plastic work
 * - EOS could be re-evaluated for consistency (currently not done)
 * - Small discrepancy acceptable for typical loading rates
 * 
 * @note Called after postprocess_prob() completes
 * @note Final operation before returning to host code
 * @note No convergence checks (assumed converged at this point)
 * 
 * @see postprocess_prob() for crystal-frame calculations
 * @see quat_to_tensor() for rotation matrix generation
 * @see get_rot_mat_vecd() for 5-vector rotation matrix
 */
template<class ProblemState, class ThermoElastN>
__ecmech_hdev__
inline
void postprocess(ProblemState& prob_state,
                 const ThermoElastN& elastN,
                 const double* const def_rate_d6v_sample,
                 double* const sdd,
                 double* const internal_energy,
                 double* const cauchy_stress_d5p_xtal,
                 double dev_strain_energy_total,
                 double halfVMidDt
                )
{
    double xtal_rmat[ecmech::ndim * ecmech::ndim];
    quat_to_tensor(xtal_rmat, prob_state.quat_u);
    //
    double rmat_5x5_sample2xtal[ecmech::ntvec * ecmech::ntvec];
    get_rot_mat_vecd(rmat_5x5_sample2xtal, xtal_rmat);
    //
    double cauchy_stress_dev_press_sample[ecmech::nsvec];
    vecsVMa<ntvec>(cauchy_stress_dev_press_sample, rmat_5x5_sample2xtal, cauchy_stress_d5p_xtal);
    cauchy_stress_dev_press_sample[iSvecS] = cauchy_stress_d5p_xtal[iSvecS];
    //
    // put end-of-step stress in cauchy_stress_d6p
    vecdsToSvecP(prob_state.cauchy_stress_d6p, cauchy_stress_dev_press_sample);
    //
    // and now the second half of the trapezoidal integration
    //
    dev_strain_energy_total += halfVMidDt * vecsInnerSvecDev(prob_state.cauchy_stress_d6p, def_rate_d6v_sample);

    // adjust sign on quat so that as close as possible to quat_o;
    // more likely to keep orientations clustered this way;
    // this flip through the origin is equivalent under antipodal symmetry
    //
    if (vecsyadotb<qdim>(prob_state.quat_u, prob_state.quat_n) < zero) {
        for (int iQ = 0; iQ < ecmech::qdim; ++iQ) {
            prob_state.quat_u[iQ] = -prob_state.quat_u[iQ];
        }
    }

    {
        double shear_modulus = elastN.getGmod(prob_state.tkelv, prob_state.pressure_EOS, prob_state.energy_new);
        sdd[i_sdd_bulk] = prob_state.bulk_modulus_new;
        sdd[i_sdd_gmod] = shear_modulus;
    }
#ifdef ECMECH_DEBUG
    assert(ecmech::nsdd == 2);
#endif

    prob_state.energy_new = prob_state.energy_new + dev_strain_energy_total;
    //
    // could update pressure and temperature again, but do not bother

    internal_energy[ecmech::i_ne_total] = prob_state.energy_new;
#ifdef ECMECH_DEBUG
    assert(ecmech::ne == 1);
#endif
}

/**
 * @brief Complete time step integration for single material point (main interface).
 * 
 * This is the primary entry point for crystal plasticity integration. It orchestrates
 * the complete sequence: preprocessing → nonlinear solve → post-processing.
 * 
 * **Function signature summary**:
 * - Input: Beginning-of-step state, prescribed kinematics, material models, time step
 * - Output: End-of-step stress, updated history variables, optional tangent
 * - Return: Success/failure flag
 * 
 * **Integration sequence**:
 * ```
 * 1. Create ProblemState from input arrays
 * 2. preprocess(): EOS, hardening, energy initialization
 * 3. Create EvptnUpdstProblem for coupled elastic-rotation solve
 * 4. Create SNLS solver and configure
 * 5. main_problem(): Solve nonlinear system
 * 6. If tangent requested: computeTangentStiffness()
 * 7. stateFromX(): Extract converged elastic strain and orientation
 * 8. Store number of function evaluations in history
 * 9. postprocess_prob(): Compute plastic work, flow strength
 * 10. postprocess(): Rotate stress, finalize energy
 * ```
 * 
 * **Template parameters**:
 * @tparam SlipGeom Slip geometry class (e.g., SlipGeomFCC, SlipGeomBCC)
 * @tparam Kinetics Kinetics class (e.g., KineticsKMBalD, KineticsVocePL)
 * @tparam ThermoElastN Thermoelastic class (e.g., ThermoElastNCubic)
 * @tparam EosModel Equation of state class (e.g., EosModelConst)
 * 
 * @param[in] slipGeom Slip geometry object defining slip systems
 * @param[in] kinetics Kinetics object for rate-dependent strength
 * @param[in] elastN Thermoelastic object for stress-strain relation
 * @param[in] eos Equation of state object for pressure-volume-energy
 * @param[in] dt Time step size [time units]
 * @param[in] tolerance Convergence tolerance for nonlinear solver [dimensionless]
 * @param[in] def_rate_d6v_sample Prescribed deformation rate [1/time, nsvp = 7 components]
 *                                 Only first ntvec = 5 used (deviatoric)
 * @param[in] spin_vec_sample Prescribed spin [1/time, ndim = 3 components]
 * @param[in] rel_vol_ratios Volume ratio array [nvr = 4 components]:
 *                           [0] = J_n, [1] = J_{n+1}, [2] = (J_{n+1}-J_n)/dt, [3] = J_{n+1}-J_n
 * @param[in,out] internal_energy Internal energy [energy/volume, ne = 1 component]
 *                                 Updated with plastic work contribution
 * @param[in,out] cauchy_stress_d6p Cauchy stress [stress units, nsvp = 7 components]
 *                                   Input: σ_n, Output: σ_{n+1}
 *                                   Layout: [σ11, σ22, σ33, σ23, σ13, σ12, -p]
 * @param[in,out] hist History variable array [numHist components]
 *                     Contains and updates: orientation, elastic strain, hardening, etc.
 * @param[in,out] tkelv Temperature [Kelvin]
 *                      May be updated by EOS evaluation
 * @param[out] sdd Auxiliary outputs [nsdd = 2 components]:
 *                 [0] = bulk modulus, [1] = shear modulus
 * @param[out] mtanSD Material tangent stiffness [nsvec2 = 36 components or nullptr]
 *                    If nullptr: Skip tangent computation
 *                    If non-null: Compute and store ∂σ/∂Δε
 * @param[in] outputLevel Diagnostic output verbosity (0 = silent)
 * 
 * @return true if integration successful, false if solver failed to converge
 * 
 * **History array layout**:
 * Defined by case-specific indices (iHistLbQ, iHistLbGdot, etc.):
 * - Orientation quaternion [4 components]
 * - Slip rates [nslip components]
 * - Elastic deviatoric strain [ntvec = 5 components]
 * - Hardening state [nH components]
 * - Equivalent plastic strain [1 component]
 * - Equivalent plastic strain rate [1 component]
 * - Flow strength [1 component]
 * - Number of function evaluations [1 component]
 * 
 * **Convergence failure**:
 * - Returns false immediately if preprocessing fails (hardening divergence)
 * - Returns false if main nonlinear solve fails to converge
 * - Host code should handle failure (reduce time step, flag element, etc.)
 * 
 * @note Thread-safe: No static storage, can be called in parallel for different points
 * @note GPU-compatible: All functions marked __ecmech_hdev__
 * @note Const correctness: Material models passed by const reference
 * 
 * @see getResponseECM() in ECMech_evptnWrap.h for vectorized batch interface
 * @see ProblemState constructor for history array interpretation
 */
template<class SlipGeom, class Kinetics, class ThermoElastN, class EosModel>
__ecmech_hdev__
inline
bool getResponseSngl(const SlipGeom& slipGeom,
                     const Kinetics& kinetics,
                     const ThermoElastN& elastN,
                     const EosModel& eos,
                     const double dt,
                     const double tolerance,
                     const double* const def_rate_d6v_sample, // defRate,
                     const double* const spin_vec_sample, // spin
                     const double* const rel_vol_ratios,
                     double* const internal_energy,
                     double* const cauchy_stress_d6p,
                     double* const hist,
                     double& tkelv,
                     double* const sdd,
                     double* const mtanSD,
                     int outputLevel = 0)
{
    auto prob_state = ProblemState<SlipGeom, Kinetics, ThermoElastN, EosModel>(hist, cauchy_stress_d6p, tkelv, def_rate_d6v_sample, spin_vec_sample, rel_vol_ratios, dt);

    double halfVMidDt, dev_strain_energy_total;
    const bool pre_status = preprocess(slipGeom, kinetics, eos, elastN, rel_vol_ratios, internal_energy, def_rate_d6v_sample, prob_state, halfVMidDt, dev_strain_energy_total);

    if (!pre_status) { return false; }

    double cauchy_stress_d5p_xtal[ecmech::nsvec];
    {
        EvptnUpdstProblem prob(slipGeom, kinetics, elastN, prob_state);

        // Solver update of things
        {
            snls::SNLSTrDlDenseG<decltype(prob)> solver(prob);
            bool status = main_problem(tolerance, solver, outputLevel);

            if (!status) {
                return false;
            }

            if (mtanSD != nullptr) {
                computeTangentStiffness(prob, solver, prob_state, mtanSD);
            }
            // store updated state
            //
            prob.stateFromX(prob_state.elast_d5_u, prob_state.quat_u, solver._x);
            //
            hist[iHistA_nFEval] = solver.getNFEvals(); // does _not_ include updateH iterations
        }
        postprocess_prob<Kinetics::nH>(prob, prob_state, cauchy_stress_d5p_xtal);
    }
    postprocess(prob_state, elastN, def_rate_d6v_sample, sdd, internal_energy, cauchy_stress_d5p_xtal, dev_strain_energy_total, halfVMidDt);
    return true;
} // getResponseSngl

#if defined(ECMECH_EXTRA_SOLVERS)
/**
 * @brief Alternative integration with simplified rotation treatment.
 * 
 * Similar to getResponseSngl() but uses EvptnNRUpdstProblem which solves only
 * for elastic strain (5 DOFs instead of 8). Lattice rotation handled separately
 * or extrapolated.
 * 
 * **When to use**:
 * - Small rotation increments (quasi-static, small strain)
 * - Faster convergence more important than rotation accuracy
 * - Debugging: isolate elastic strain solver issues
 * - Operator-split schemes with separate rotation update
 * 
 * **Differences from getResponseSngl**:
 * - Uses EvptnNRUpdstProblem instead of EvptnUpdstProblem
 * - Solves smaller system (5 unknowns vs 8)
 * - Potentially faster convergence (fewer coupling terms)
 * - Less accurate for large rotations
 * 
 * **Compile flag**:
 * Only available if ECMECH_EXTRA_SOLVERS is defined.
 * 
 * **Optional RStarSolve**:
 * If template parameter RStarSolve=true, performs secondary rotation solve
 * after elastic strain converges for improved accuracy.
 * 
 * **Template parameters**:
 * @tparam SlipGeom Slip geometry class (e.g., SlipGeomFCC, SlipGeomBCC)
 * @tparam Kinetics Kinetics class (e.g., KineticsKMBalD, KineticsVocePL)
 * @tparam ThermoElastN Thermoelastic class (e.g., ThermoElastNCubic)
 * @tparam EosModel Equation of state class (e.g., EosModelConst)
 * 
 * @param[in] slipGeom Slip geometry object defining slip systems
 * @param[in] kinetics Kinetics object for rate-dependent strength
 * @param[in] elastN Thermoelastic object for stress-strain relation
 * @param[in] eos Equation of state object for pressure-volume-energy
 * @param[in] dt Time step size [time units]
 * @param[in] tolerance Convergence tolerance for nonlinear solver [dimensionless]
 * @param[in] def_rate_d6v_sample Prescribed deformation rate [1/time, nsvp = 7 components]
 *                                 Only first ntvec = 5 used (deviatoric)
 * @param[in] spin_vec_sample Prescribed spin [1/time, ndim = 3 components]
 * @param[in] rel_vol_ratios Volume ratio array [nvr = 4 components]:
 *                           [0] = J_n, [1] = J_{n+1}, [2] = (J_{n+1}-J_n)/dt, [3] = J_{n+1}-J_n
 * @param[in,out] internal_energy Internal energy [energy/volume, ne = 1 component]
 *                                 Updated with plastic work contribution
 * @param[in,out] cauchy_stress_d6p Cauchy stress [stress units, nsvp = 7 components]
 *                                   Input: σ_n, Output: σ_{n+1}
 *                                   Layout: [σ11, σ22, σ33, σ23, σ13, σ12, -p]
 * @param[in,out] hist History variable array [numHist components]
 *                     Contains and updates: orientation, elastic strain, hardening, etc.
 * @param[in,out] tkelv Temperature [Kelvin]
 *                      May be updated by EOS evaluation
 * @param[out] sdd Auxiliary outputs [nsdd = 2 components]:
 *                 [0] = bulk modulus, [1] = shear modulus
 * @param[out] mtanSD Material tangent stiffness [nsvec2 = 36 components or nullptr]
 *                    If nullptr: Skip tangent computation
 *                    If non-null: Compute and store ∂σ/∂Δε
 * @param[in] outputLevel Diagnostic output verbosity (0 = silent)
 * 
 * @return true if integration successful, false if solver failed to converge
 * 
 * @see getResponseSngl() for full coupled integration
 * @see EvptnNRUpdstProblem for problem formulation
 */
template<class SlipGeom, class Kinetics, class ThermoElastN, class EosModel>
__ecmech_hdev__
inline
bool getResponseNRSngl(
                     const SlipGeom& slipGeom,
                     const Kinetics& kinetics,
                     const ThermoElastN& elastN,
                     const EosModel& eos,
                     const double dt,
                     const double tolerance,
                     const double* const def_rate_d6v_sample, // defRate,
                     const double* const spin_vec_sample, // spin
                     const double* const rel_vol_ratios,
                     double* const internal_energy,
                     double* const cauchy_stress_d6p,
                     double* const hist,
                     double& tkelv,
                     double* const sdd,
                     double* const mtanSD,
                     int outputLevel = 0)
{
    auto prob_state = ProblemState<SlipGeom, Kinetics, ThermoElastN, EosModel>(hist, cauchy_stress_d6p, tkelv, def_rate_d6v_sample, spin_vec_sample, rel_vol_ratios, dt);

    double halfVMidDt, dev_strain_energy_total;
    preprocess<SlipGeom, Kinetics, EosModel, ThermoElastN, decltype(prob_state), true>(slipGeom, kinetics, eos, elastN, rel_vol_ratios, internal_energy, def_rate_d6v_sample, prob_state, halfVMidDt, dev_strain_energy_total);

    double cauchy_stress_d5p_xtal[ecmech::nsvec];
    {
        EvptnNRUpdstProblem prob(slipGeom, kinetics, elastN, prob_state);

        // Solver update of things
        {
            snls::SNLSTrDlDenseG<decltype(prob)> solver(prob);
            bool status = main_problem(tolerance, solver, outputLevel);

            if (!status) {
                return false;
            }

            if (mtanSD != nullptr) {
                computeTangentStiffness(prob, solver, prob_state, mtanSD);
            }
            // store updated state
            //
            prob.stateFromX(prob_state.elast_d5_u, solver._x);
            //
            hist[iHistA_nFEval] = solver.getNFEvals(); // does _not_ include updateH iterations
        }
        postprocess_prob<Kinetics::nH>(prob, prob_state, cauchy_stress_d5p_xtal);
    }
    postprocess(prob_state, elastN, def_rate_d6v_sample, sdd, internal_energy, cauchy_stress_d5p_xtal, dev_strain_energy_total, halfVMidDt);
    return true;
} // getResponseSngl
#endif

}
}