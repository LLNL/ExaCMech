/**
 * @file ECMech_base_classes.h
 * @brief Core problem formulation classes for elasto-viscoplastic crystal plasticity integration.
 * 
 * This file defines the fundamental building blocks for implicit time integration of
 * crystal plasticity constitutive models. It provides classes that encapsulate:
 * - Material state evolution (ProblemState)
 * - Elastic lattice strain updates (EvptnLatticeStrainProblem)
 * - Lattice rotation updates (EvptnLatticeRotationProblem)
 * 
 * These classes serve as subproblems within the larger crystal plasticity integration
 * framework, separating concerns between:
 * 1. **Elastic deformation**: Evolution of elastic strain in the lattice frame
 * 2. **Lattice rotation**: Evolution of crystal orientation quaternions
 * 3. **State management**: Tracking all material state variables
 * 
 * **Integration approach**:
 * - Implicit backward Euler time integration
 * - Coupled nonlinear system for elastic strain and rotation
 * - Residual and Jacobian formulations for SNLS solvers
 * - Scaling for numerical conditioning
 * 
 * **Key features**:
 * - Template-based design for flexibility with different constitutive models
 * - GPU-compatible (__ecmech_hdev__ functions)
 * - Efficient tensor operations via ECMech_util.h functions
 * - Support for dynamic slip systems (crystal plasticity variants)
 * 
 * @see ECMech_evptn.h for higher-level problem formulations
 * @see ECMech_evptnSngl.h for single material point integration routines
 * @see ECMech_elastic.h for thermoelastic constitutive models
 */

#pragma once

#include <cassert>

#include "ECMech_core.h"
#include "ECMech_util.h"

#include "RAJA/RAJA.hpp"

namespace ecmech {
namespace evptn {

    /**
     * @brief Number of auxiliary history variables.
     * 
     * These are derived/diagnostic quantities computed after convergence:
     * - Effective shear rate (eps_dot or dplas_eff)
     * - Accumulated effective shear (eps or shr_eff)
     * - Flow strength
     * - Number of function evaluations
     * 
     * Value: 4
     */
    const int numHistAux = 4; // effective shearing rate, accumulated shear, flow strength, nFEval
    /**
     * @brief Starting index for auxiliary history variables.
     * 
     * Auxiliary variables are stored first in the history array.
     * Value: 0
     */
    const int iHistLbA = 0;
    /**
     * @brief Index for effective shearing rate in history array.
     * 
     * Storage for equivalent plastic strain rate γ̇_eff or D_p_eff depending on
     * ECMECH_USE_DPEFF compile flag.
     * 
     * **Conditional naming**:
     * - If ECMECH_USE_DPEFF: "dplas_eff" (tensor-based effective rate)
     * - Otherwise: "shrate_eff" (sum of absolute slip rates)
     * @see ProblemState::eps_dot for access via state container
     */
    const int iHistA_shrateEff = iHistLbA + 0;
    /**
     * @brief Index for accumulated effective shear in history array.
     * 
     * Storage for integrated equivalent plastic strain ε_p = ∫ γ̇_eff dt.
     * 
     * **Conditional naming**:
     * - If ECMECH_USE_DPEFF: "eps" (effective plastic strain from tensor norm)
     * - Otherwise: "shr_eff" (effective shear from slip rate sum)
     * @see ProblemState::eps for access via state container
     */
    const int iHistA_shrEff = iHistLbA + 1;
    /**
     * @brief Index for flow strength in history array.
     * 
     * Storage for generalized "yield stress" computed as:
     * @see ProblemState::flow_strength for access via state container
     */
    const int iHistA_flowStr = iHistLbA + 2;
    /**
     * @brief Index for number of function evaluations in history array.
     * 
     * Storage for iteration count from nonlinear solver (SNLS).
     * @see SNLS solver getNFEvals() for how this is populated
     */
    const int iHistA_nFEval = iHistLbA + 3;
    /**
     * @brief Starting index for elastic deviatoric strain in history array.
     * 
     * Elastic strain stored as deviatoric 5-vector in crystal/lattice frame.
     * Occupies indices [iHistLbE : iHistLbE + ntvec - 1] = [4:8].
     * @see ProblemState::elast_d5_n and ProblemState::elast_d5_u
     */
    const int iHistLbE = numHistAux;
    /**
     * @brief Starting index for crystal orientation quaternion in history array.
     * 
     * Orientation stored as unit quaternion representing rotation mapping crystal to sample.
     * Occupies indices [iHistLbQ : iHistLbQ + qdim - 1] = [9:12].
     * @see ProblemState::quat_n and ProblemState::quat_u
     */
    const int iHistLbQ = numHistAux + ecmech::ntvec;
    /**
     * @brief Starting index for hardening state variables in history array.
     * 
     * Hardening variables are material-model specific (from Kinetics class).
     * Occupies indices [iHistLbH : iHistLbH + Kinetics::nH - 1].
     * @see ProblemState::h_state and ProblemState::h_state_u
     * @see Kinetics classes for specific hardening variable definitions
     */
    const int iHistLbH = numHistAux + ecmech::ntvec + ecmech::qdim;

    /**
     * @brief Traits class for computing history array dimensions.
     * 
     * This template class provides compile-time constants for:
     * - Total number of history variables (numHist)
     * - Starting index for slip rates (iHistLbGdot)
     * 
     * **Design pattern**: Traits class
     * - No data members, no methods
     * - Only static constexpr values
     * - Template specialization based on constitutive model combination
     * 
     * **Purpose**:
     * - Centralize history array layout calculations
     * - Ensure consistency across code
     * - Enable compile-time sizing of arrays
     * 
     * **Template parameters**:
     * @tparam SlipGeom Slip geometry class (defines nslip)
     * @tparam Kinetics Kinetics class (defines nH)
     * @tparam ThermoElastN Thermoelastic model (not used in calculations)
     * @tparam EosModel Equation of state (not used in calculations)
     * 
     * **History array layout** (for reference):
     * ```
     * [0:3]              : Auxiliary (shrateEff, shrEff, flowStr, nFEval)
     * [4:8]              : Elastic deviatoric strain (ntvec = 5)
     * [9:12]             : Orientation quaternion (qdim = 4)
     * [13:13+nH-1]       : Hardening state (Kinetics::nH components)
     * [13+nH:13+nH+nslip-1] : Slip rates (SlipGeom::nslip components)
     * ```
     * Total: numHist = 13 + nH + nslip
     * 
     * @see matModel for usage in material model wrapper
     * @see ProblemState constructor for how indices are used
     */
    template<class SlipGeom, class Kinetics, class ThermoElastN, class EosModel>
    class NumHist
    {
        public:
        /**
         * @brief Starting index for slip rates in history array.
         * 
         * Slip rates stored after auxiliary variables, elastic strain, orientation,
         * and hardening state.
         * 
         * Value: iHistLbH + Kinetics::nH = 13 + nH
         * @see ProblemState::gdot for access
         */
        static constexpr int iHistLbGdot = iHistLbH + Kinetics::nH;
        /**
         * @brief Total number of history variables.
         * 
         * Includes all auxiliary, state, and diagnostic variables:
         * - 4 auxiliary (shrateEff, shrEff, flowStr, nFEval)
         * - 5 elastic strain components (ntvec)
         * - 4 quaternion components (qdim)
         * - nH hardening variables (Kinetics::nH)
         * - nslip slip rates (SlipGeom::nslip)
         * 
         * Value: iHistLbH + Kinetics::nH + SlipGeom::nslip = 13 + nH + nslip
         */
        static constexpr int numHist = iHistLbH + Kinetics::nH + SlipGeom::nslip;
    }; // NumHist

    /**
     * @brief State variable container for a single material point during time integration.
     * 
     * This template class aggregates all state variables required for crystal plasticity
     * integration at a material point, including:
     * - Crystal orientation (quaternion)
     * - Elastic strain in lattice frame
     * - Hardening state variables
     * - Slip rates
     * - Stress state
     * - Thermodynamic variables (temperature, pressure, energy)
     * - Deformation kinematics
     * 
     * **Design philosophy**:
     * - Provides a unified interface to heterogeneous state data
     * - Separates "n" (beginning-of-step) from "u" (updated/end-of-step) values
     * - Handles mapping between host code arrays and local variables
     * - Supports both static and dynamic slip system configurations
     * 
     * **Template parameters**:
     * @tparam SlipGeom Slip geometry class defining slip systems (e.g., SlipGeomFCC)
     * @tparam Kinetics Kinetics class for rate-dependent strength (e.g., KineticsKMBalD)
     * @tparam ThermoElastN Thermoelastic model (e.g., ThermoElastNCubic)
     * @tparam EosModel Equation of state model (e.g., EosModelConst)
     * 
     * **Memory layout**:
     * State variables are extracted from a flat history array (histV) according to
     * predefined index conventions (iHistLbQ, iHistLbH, etc.)
     * 
     * @see getResponseSngl() for usage in time integration
     * @see EvptnUpdstProblem for the nonlinear problem formulation using this state
     */
    template<class SlipGeom, class Kinetics, class ThermoElastN, class EosModel>
    struct ProblemState
    {
        static constexpr int iHistLbGdot = NumHist<SlipGeom, Kinetics, ThermoElastN, EosModel>::iHistLbGdot;

        /**
         * @brief Hardening state variables at beginning of step [various units].
         * Material-specific internal variables (e.g., slip resistances, dislocation densities).
         * Size: Kinetics::nH
         */
        double* const h_state;
        /**
         * @brief Slip rates on all slip systems [1/time].
         * Computed from resolved shear stresses and kinetic law.
         * Size: SlipGeom::nslip
         */
        double* const gdot;
        /**
         * @brief Updated elastic strain at end of step, deviatoric 5-vector [dimensionless].
         * Result of implicit integration, stored here after convergence.
         */
        double* const elast_d5_u;
        /**
         * @brief Updated crystal orientation quaternion at end of step [unit quaternion].
         * Result of lattice rotation integration.
         */
        double* const quat_u;
        /**
         * @brief Equivalent plastic strain rate [1/time].
         * Current rate of plastic deformation accumulation.
         */
        double& eps_dot;
        /**
         * @brief Accumulated equivalent plastic strain [dimensionless].
         * Integrated measure of total plastic deformation.
         */
        double& eps;
        /**
         * @brief Flow strength [stress units].
         * Current resistance to plastic deformation (generalized "yield stress").
         */
        double& flow_strength;
        /**
         * @brief Cauchy stress in deviatoric-pressure form [stress units, 7 components].
         * Layout: [σ₁₁, σ₂₂, σ₃₃, σ₂₃, σ₁₃, σ₁₂, -p]
         * where p = -(σ₁₁ + σ₂₂ + σ₃₃)/3 is the pressure.
         * This array is both input (n) and output (n+1).
         * Sample frame
         */
        double* const cauchy_stress_d6p;
        /**
         * @brief Spin vector in sample frame [1/time, 3 components].
         * Prescribed velocity gradient skew-symmetric part, axial vector form.
         * Components: [W₁₂, W₂₃, W₁₃] from skew tensor.
         */
        const double* const spin_vec_sample;
        /**
         * @brief Relative volume at end of step (J_{n+1} = V_{n+1} / V_ref).
         * Updated value after deformation increment.
         */
        const double rel_vol_new;
        /**
         * @brief Time step size [time units].
         * Fundamental parameter controlling integration accuracy and stability.
         */
        const double dt;
        /**
         * @brief Temperature [Kelvin].
         * Thermodynamic state variable affecting kinetics and thermoelasticity.
         */
        double& tkelv;

        /**
         * @brief Deformation rate in sample frame, deviatoric 5-vector [1/time].
         * Prescribed velocity gradient symmetric part, traceless representation.
         * Uses √2 and √3 normalization for norm preservation.
         */
        double def_rate_d5_sample[ecmech::ntvec];
        /**
         * @brief Elastic strain at beginning of step, deviatoric 5-vector [dimensionless].
         * Lattice frame elastic deformation, traceless symmetric representation.
         */
        double elast_d5_n[ecmech::ntvec];
        /**
         * @brief Crystal orientation quaternion at beginning of step [unit quaternion].
         * Represents rotation from lattice frame to sample frame.
         * Layout: [q₀, q₁, q₂, q₃] with |q| = 1.
         */
        double quat_n[ecmech::qdim];
        /**
         * @brief Updated hardening state at end of step [various units].
         * Result of hardening law integration.
         * Size: Kinetics::nH
         */
        double h_state_u[Kinetics::nH];
        /**
         * @brief EOS pressure contribution [stress units].
         * Pressure from equation of state, separate from elastic deviatoric stress.
         */
        double pressure_EOS;
        /**
         * @brief Internal energy at end of step [energy/volume].
         * Thermodynamic state variable updated from mechanical work.
         */
        double energy_new;
        /**
         * @brief Bulk modulus at current state [stress units].
         * Volumetric stiffness, may be temperature and pressure dependent.
         */
        double bulk_modulus_new;

        /**
         * @brief Constructor: Extract state from history array and initialize working variables.
         * 
         * This constructor performs the critical task of unpacking the flat history array
         * (histV) used by the host code into structured state variables. It also initializes
         * beginning-of-step values and prepares storage for updated values.
         * 
         * **Initialization sequence**:
         * 1. Extract time step and volume ratios
         * 2. Map pointers to orientation, elastic strain, and hardening state in histV
         * 3. Copy Cauchy stress for modification
         * 4. Initialize temperature reference
         * 5. Store prescribed kinematics
         * 6. Normalize crystal orientation quaternion (numerical safety)
         * 
         * @param[in,out] hist Flat history array from host code [numHist elements]
         *                     Modified during integration, stores both input and output
         * @param[in] cauchy_stress_d6p_in Stress state [7 elements], copied to member
         * @param[in] tkelv_in Temperature [K], stored by reference
         * @param[in] def_rate_d6v_sample Deformation rate [nsvp elements], only [0:4] used
         * @param[in] spin_vec_sample Spin vector [ndim elements]
         * @param[in] rel_vol_ratios Volume ratio array [nvr elements]:
         *            [0]: J_n (old), [1]: J_{n+1} (new), [2]: (J_{n+1}-J_n)/dt, [3]: J_{n+1}-J_n
         * @param[in] dt_in Time step size [time units]
         * 
         * **History array layout** (indices from ECMech_cases.h):
         * - iHistLbQ: Start of quaternion (4 components)
         * - iHistLbGdot: Start of slip rates (nslip components) 
         * - iHistLbEdev: Start of elastic deviatoric strain (ntvec components)
         * - iHistLbH: Start of hardening state (nH components)
         * - iHistA_ep: Scalar equivalent plastic strain
         * - iHistA_epdot: Scalar equivalent plastic strain rate
         * - iHistA_flow_stress: Scalar flow strength
         * 
         * @note Quaternion normalization ensures numerical stability over many increments
         * @note Const correctness: input arrays remain unchanged, outputs clearly marked
         */
        __ecmech_hdev__
        ProblemState(double* const hist, double* const cauchy_stress_d6p,
                    double& tkelv,
                    const double* const def_rate_d6v_sample,
                    const double* const spin_vec_sample,
                    const double* const rel_vol_ratios,
                    const double dt) :
        h_state(&(hist[iHistLbH])),
        gdot(&(hist[iHistLbGdot])),
        elast_d5_u(&(hist[iHistLbE])),
        quat_u(&(hist[iHistLbQ])),
        eps_dot(hist[iHistA_shrateEff]),
        eps(hist[iHistA_shrEff]),
        flow_strength(hist[iHistA_flowStr]),
        cauchy_stress_d6p(cauchy_stress_d6p),
        spin_vec_sample(spin_vec_sample),
        rel_vol_new(rel_vol_ratios[1]),
        dt(dt),
        tkelv(tkelv)
        {
            // convert deformation rate convention
            //
            // double def_rate_d5_sample[ecmech::ntvec];
            svecToVecd(def_rate_d5_sample, def_rate_d6v_sample);
            //
            // copies, to keep beginning-of-step state safe
            //
            for (int i_hist = 0; i_hist < ecmech::ntvec; i_hist++) {
                elast_d5_n[i_hist] = hist[iHistLbE + i_hist];
            }

            for (int i_hist = 0; i_hist < ecmech::qdim; i_hist++) {
                quat_n[i_hist] = hist[iHistLbQ + i_hist];
            }
            //
            // normalize quat just in case
            vecsVNormalize<qdim>(quat_n);
        }

        /**
         * @brief Destructor (default).
         */
        ~ProblemState() = default;
    };

    /**
     * @brief Subproblem formulation for elastic lattice strain evolution.
     * 
     * This template class encapsulates the implicit time integration of elastic strain
     * in the crystal lattice frame. It provides:
     * - Residual evaluation for strain evolution equation
     * - Jacobian computation for Newton-Raphson solver
     * - Stress evaluation from elastic strain
     * - Material tangent stiffness contributions
     * 
     * **Governing equation** (implicit backward Euler):
     * ```
     * Residual: R_e = Δε_e - Δt(D_sample - D_plastic)
     * ```
     * where:
     * - Δε_e: Increment in elastic strain (unknowns)
     * - D_sample: Prescribed deformation rate in sample frame (input)
     * - D_plastic: Plastic deformation rate from slip (depends on ε_e via stress)
     * 
     * **Elastic strain representation**:
     * - Deviatoric 5-vector form (traceless symmetric tensor)
     * - Lattice/crystal reference frame
     * - Volumetric logarithmic strain measure (ln V_e basis)
     * - Volumetric part handled separately via EOS
     * 
     * **Coupling to other subproblems**:
     * - Stress from elastic strain drives slip kinetics
     * - Slip rates determine plastic deformation rate
     * - Lattice rotation affects frame transformations
     * 
     * @tparam ThermoElastN Thermoelastic model class (e.g., ThermoElastNCubic, ThermoElastNHexag)
     *                      Must provide eval(), getCauchy(), multDTDepsT() methods
     * 
     * **Numerical considerations**:
     * - Scaling by e_scale for conditioning (x = Δε_e / e_scale)
     * - Inverse time step pre-computed for efficiency
     * - Volume scaling (a_vol = J^(1/3)) for finite strain kinematics
     * 
     * @see EvptnUpdstProblem for coupled elastic strain + rotation problem
     * @see ThermoElastN classes for specific anisotropic elasticity models
     */
    template<class ThermoElastN>
    class EvptnLatticeStrainProblem
    {
        public:
        /**
         * @brief System dimension: number of deviatoric elastic strain components.
         * Value: 5 (deviatoric symmetric 3×3 tensor has 5 independent components)
         */
        static constexpr int nDimSys = ecmech::ntvec;

        /**
         * @brief Constructor: Initialize elastic strain subproblem with state data.
         * 
         * Sets up all parameters needed for elastic strain integration over one time step.
         * Pre-computes inverse and volume scaling quantities for efficiency.
         * 
         * @param[in] thermoElastN Thermoelastic model providing stress-strain relation
         * @param[in] dt Time step size [time units]
         * @param[in] det_v_e Elastic volume determinant J_e = det(F_e)
         * @param[in] energy_vol_ref Internal energy at end of time step [energy]
         * @param[in] pressure_EOS Pressure from equation of state [stress units]
         * @param[in] tkelv Temperature [Kelvin]
         * @param[in] elast_d5_n Beginning-of-step elastic strain [5 components]
         * 
         * **Computed derived quantities**:
         * - m_inv_dt = 1/dt (rate conversion)
         * - m_inv_det_v_e = 1/J_e (stress push-forward)
         * - m_a_vol = J_e^(1/3) (isotropic volume scaling)
         * - m_inv_a_vol = J_e^(-1/3) (inverse scaling)
         */
        __ecmech_hdev__
        EvptnLatticeStrainProblem(const ThermoElastN& thermoElastN,
                                const double dt,
                                const double det_v_e, 
                                const double energy_vol_ref, 
                                const double pressure_EOS, 
                                const double tkelv,
                                const double* const elast_d5_n)
        : m_thermo_elast_n(thermoElastN),
        m_dt(dt), m_det_v_e(det_v_e), m_energy_vol_ref(energy_vol_ref),
        m_pressure_EOS(pressure_EOS), m_tkelv(tkelv),
        m_elast_d5_n(elast_d5_n),
        m_inv_dt(1.0 / dt),
        m_inv_det_v_e(1.0 / m_det_v_e),
        m_a_vol(pow(m_det_v_e, onethird)),
        m_inv_a_vol(1.0 / m_a_vol)
        {}

        /**
         * @brief Destructor (default).
         */
        ~EvptnLatticeStrainProblem() = default;

        /**
         * @brief Convert elastic strain to Kirchhoff stress.
         * 
         * Evaluates Kirchhoff stress τ = J σ from elastic deviatoric strain using:
         * 1. Scale strain by inverse volume: ε_scaled = (1/a_vol) * ε_dev
         * 2. Add volumetric part: ε_scaled[5] = √3 * ln(a_vol)
         * 3. Apply elastic constitutive law: τ = K(ε_scaled, T, p_EOS, E)
         * 
         * **Kirchhoff vs. Cauchy stress**:
         * - Kirchhoff: τ = J σ (weighted by volume)
         * - Cauchy: σ = (1/J) τ (true stress)
         * - Kirchhoff is natural for constitutive models in finite strain
         * 
         * @param[out] kirchoff_stress Kirchhoff stress tensor [stress units, 6 components]
         *                              Voigt notation: [τ₁₁, τ₂₂, τ₃₃, τ₂₃, τ₁₃, τ₁₂]
         * @param[in] elast_d5 Deviatoric elastic strain [dimensionless, 5 components]
         * 
         * @see elast_strain_to_cauchy_stress() for direct Cauchy stress evaluation
         */
        __ecmech_hdev__
        inline
        void elast_strain_to_kirchoff_stress(double* const kirchoff_stress, // nsvec
                                            const double* const elast_d5 // ntvec
                                        ) const
        {
        //// do not need to use elaw_T_BT here as T and BT are the same
        //
        // specialize to cem%l_lin_lnsd
        // CALL elawn_T(s_meas, elast_d5_f, crys%elas, tkelv, .TRUE., a_V, &
        // & pressure_EOS, energy_vol_ref, crys%i_eos_model, crys%eos_const &
        // &)
        double elast_d5v[ecmech::nsvec];
        vecsVxa<ntvec>(elast_d5v, m_inv_a_vol, elast_d5);
        //// tr_Ee = three * DLOG(a_V%r)
        //// CALL trace_to_vecds_s(s_meas%elast_dev_press_vec(SVEC), tr_Ee)
        elast_d5v[iSvecS] = sqr3 * log(m_a_vol); // could go into constructor
        //
        //// Kirchhoff stress from elast_d5v
        // CALL elawn_lin_op(s_meas%kirchoff, s_meas%elast_dev_press_vec, cem, tkelv, &
        // & pressure_EOS, energy_vol_ref, i_eos_model, eos_const)
        m_thermo_elast_n.eval(kirchoff_stress, elast_d5v, m_tkelv, m_pressure_EOS, m_energy_vol_ref);
        }

        /**
         * @brief Convert elastic strain to Cauchy stress.
         * 
         * Convenience wrapper that computes Cauchy stress σ directly:
         * 1. Compute Kirchhoff stress: τ = f(ε_e)
         * 2. Push forward to Cauchy: σ = (1/J) τ
         * 
         * @param[out] cauchy Cauchy stress tensor [stress units, 6 components]
         * @param[in] elast_d5_f Deviatoric elastic strain [dimensionless, 5 components]
         * 
         * @note "f" suffix typically denotes "final" (end-of-step) value
         */
        __ecmech_hdev__
        inline
        void elast_strain_to_cauchy_stress(double* const cauchy, // nsvec
                                        const double* const elast_d5_f // ntvec
                                        ) const
        {
        double kirchoff[ecmech::nsvec];
        this->elast_strain_to_kirchoff_stress(kirchoff, elast_d5_f);
        m_thermo_elast_n.getCauchy(cauchy, kirchoff, m_inv_det_v_e);
        }

        /**
         * @brief Extract elastic strain and rate from solution vector.
         * 
         * This function decodes the scaled solution vector x into physical
         * elastic strain and optionally computes the strain rate.
         * 
         * **Scaling convention**:
         * - Unknowns stored as: x = Δε_e / e_scale
         * - Physical increment: Δε_e = e_scale * x
         * - End-of-step: ε_e^{n+1} = ε_e^n + Δε_e
         * - Rate: ε̇_e = Δε_e / Δt
         * 
         * @tparam calc_strain_rate If true, also compute ε̇_e (default: false)
         * 
         * @param[out] elast_delta_d5 Elastic strain increment Δε_e [dimensionless, 5 components]
         * @param[out] elast_dt_d5 Elastic strain rate ε̇_e [1/time, 5 components]
         *                          Only computed if calc_strain_rate=true
         * @param[in] x Scaled solution vector [dimensionless, ntvec components]
         * 
         * **Usage in Jacobian evaluation**:
         * ```cpp
         * get_elast_strain_state<true>(delta, rate, x);  // Need rate for derivatives
         * ```
         */
        template<bool calc_strain_rate = false>
        __ecmech_hdev__
        inline
        void get_elast_strain_state(double* const elast_delta_d5,
                                double* const elast_dt_d5,
                                const double* const x) const
        {
        //////////////////////////////
        // PULL VALUES out of x, with scalings
        //
        // double elast_dt_d5[ecmech::ntvec];
        vecsVxa<ntvec>(elast_dt_d5, ecmech::e_scale, x); // elast_dt_d5 is now the delta, _not_ yet elast_dt_d5
        // elast_d5_f is end-of-step
        // double elast_d5_f[ntvec];
        vecsVapb<ntvec>(elast_delta_d5, elast_dt_d5, m_elast_d5_n);
        if constexpr(calc_strain_rate) {
            vecsVsa<ntvec>(elast_dt_d5, m_inv_dt); // _now_ elast_dt_d5 has dt contributions
        }
        }

        /**
         * @brief Reconstruct end-of-step elastic strain from solution vector.
         * 
         * Simpler interface to get_elast_strain_state<false> for when only
         * the final strain value is needed (not the rate).
         * 
         * @param[out] elast_d5 End-of-step elastic strain ε_e^{n+1} [dimensionless, 5 components]
         * @param[in] x Scaled solution vector [dimensionless, ntvec components]
         * 
         * @warning Not safe if elast_d5 and m_elast_d5_n point to same memory
         */
        __ecmech_hdev__
        inline
        void stateFromX(double* const elast_d5,
                    const double* const x) const
        {
        double elast_d5_delta[ecmech::ntvec] = {};
        this->get_elast_strain_state(elast_d5, elast_d5_delta, x);
        }

        /**
         * @brief Evaluate elastic strain residual equation.
         * 
         * Computes the residual vector for the elastic strain evolution equation:
         * ```
         * R_e = (ε̇_e - D_sample + D_plastic) * epsdot_scale_inv
         * ```
         * This should be zero at convergence, meaning elastic rate equals total
         * rate minus plastic rate.
         * 
         * **Physical interpretation**:
         * - D_sample: Prescribed total deformation rate (input)
         * - D_plastic: Plastic deformation from slip (computed from ε_e via stress)
         * - ε̇_e: Elastic rate (unknown, from time difference of ε_e)
         * - Balance: ε̇_e = D_sample - D_plastic (compatibility)
         * 
         * @param[out] residual Residual vector [dimensionless, ntvec components]
         *                       Scaled for numerical conditioning
         * @param[in] epsdot_scale_inv Inverse strain rate scaling (1/D_ref)
         * @param[in] elast_dt_d5 Elastic strain rate ε̇_e [1/time, ntvec components]
         * @param[in] plastic_def_rate_d5 Plastic deformation rate D_p [1/time, ntvec components]
         * @param[in] def_rate_d5_xtal Total deformation rate in crystal frame [1/time, ntvec components]
         * 
         * **Scaling rationale**:
         * - epsdot_scale_inv ≈ 1/‖D_sample‖ makes residual O(1) for better conditioning
         * - Helps Newton-Raphson convergence regardless of loading rate magnitude
         */
        __ecmech_hdev__
        inline
        void get_elast_strain_residual(double* const residual,
                                    const double epsdot_scale_inv,
                                    const double* const elast_dt_d5,
                                    const double* const plastic_def_rate_d5,
                                    const double* const def_rate_d5_xtal) const
        {
        for (size_t iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
            residual[m_ind_sub_elas + iTvec] = epsdot_scale_inv * ( // SCALING
                m_inv_a_vol * elast_dt_d5[iTvec] + plastic_def_rate_d5[iTvec] - def_rate_d5_xtal[iTvec]);
        }
        }

        /**
         * @brief Compute Jacobian contribution: ∂R_e/∂ε_e (elastic strain block).
         * 
         * Evaluates the derivative of the elastic residual with respect to elastic
         * strain unknowns. This forms the diagonal block in the full Jacobian matrix.
         * 
         * **Mathematical form**:
         * ```
         * J_ee = ∂R_e/∂ε_e = (I/Δt - ∂D_p/∂ε_e) * epsdot_scale_inv * e_scale
         * ```
         * where:
         * - I/Δt: Identity scaled by inverse time step (explicit elastic rate term)
         * - ∂D_p/∂ε_e: Derivative of plastic rate w.r.t. strain (implicit slip coupling)
         * - Scalings: epsdot_scale_inv (rate conditioning), e_scale (unknown scaling)
         * 
         * **Implementation notes**:
         * - Uses RAJA::View for 2D array access (row-major indexing)
         * - Adds identity scaled contribution first
         * - Subtracts plastic rate derivative (chain rule through stress)
         * - Applied scaling transformations for solver conditioning
         * 
         * @tparam JAC_SIZE Total Jacobian dimension (may include rotation DOFs)
         * 
         * @param[in,out] jacobian Full Jacobian matrix [JAC_SIZE × JAC_SIZE]
         *                          Modified to include elastic-elastic block
         * @param[in] dDp_hat_delast_strain Derivative ∂D_p/∂ε_e [ntvec × ntvec]
         * @param[in] epsdot_scale_inv Inverse strain rate scaling factor
         * 
         * **Memory layout**:
         * Jacobian accessed as: jacobian[i*JAC_SIZE + j] for element (i,j)
         * Elastic block: rows 0:ntvec-1, columns 0:ntvec-1
         */
        template<size_t JAC_SIZE>
        __ecmech_hdev__
        inline
        void get_deriv_elast_strain_wrt_elast_strain(double* const jacobian,
                                                    const double* const dDp_hat_delast_strain) const
        {
        RAJA::View<double, RAJA::Layout<2>> jacob_ee(jacobian, JAC_SIZE, JAC_SIZE);
        // dislocation plasticity;
        // first contribution; overwrite
        for (size_t jTvec = 0; jTvec < ecmech::ntvec; ++jTvec) {
            for (size_t iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                jacob_ee(iTvec, jTvec) = dDp_hat_delast_strain[ECMECH_NN_INDX(iTvec, jTvec, ecmech::ntvec)];
            }
        }
        // elastic rate
        //
        {
            const double adti = m_inv_a_vol * m_inv_dt;
            for (size_t iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                jacob_ee(iTvec, iTvec) += adti;
            }
        }
        } 

        /**
         * @brief Compute Jacobian: ∂R_ω/∂ε_e (rotation residual w.r.t. elastic strain).
         * 
         * This method evaluates the off-diagonal coupling between the rotation residual
         * and elastic strain unknowns. This captures how changes in elastic strain affect
         * the lattice rotation evolution through:
         * - Changes in plastic spin (stress → slip rates → plastic spin)
         * - Elastic spin higher-order corrections
         * 
         * **Mathematical formulation**:
         * ```
         * J_re = ∂R_ω/∂ε_e = -Δt * ∂W_p/∂ε_e + ∂(elastic_corrections)/∂ε_e
         * ```
         * where:
         * - ∂W_p/∂ε_e: Derivative of plastic spin w.r.t. elastic strain
         * - Elastic corrections: Higher-order terms from finite strain kinematics
         * 
         * **Physical coupling mechanisms**:
         * 1. **Plastic spin coupling**: ε_e → τ → γ̇ → W_p
         * 2. **Elastic spin coupling**: Geometric terms from corotational derivatives
         * 
         * **Higher-order elastic terms**:
         * The method includes contributions from:
         * - A_e_M35: Elastic rotation coupling matrix [nwvec × ntvec]
         * - Time derivatives of elastic rotation (A_edot_M35)
         * - Scaled by dt_ee_fac = Δt * elast_elast_factor
         * 
         * **Template parameters**:
         * @tparam JAC_SIZE Total Jacobian dimension (8 for full coupled problem)
         * @tparam ind_sub_r Starting index for rotation DOFs in Jacobian (typically ntvec = 5)
         * 
         * @param[in,out] jacobian Full Jacobian matrix [JAC_SIZE × JAC_SIZE]
         *                         Modified to include rotation-elastic strain coupling block
         * @param[in] elast_dt_d5 Elastic strain rate ε̇_e [1/time, ntvec = 5]
         *                        Used to compute A_edot_M35 time derivative
         * @param[in] elast_elast_factor Scaling factor for elastic-elastic contributions
         *                               Controls magnitude of higher-order geometric terms
         * @param[in] dWp_hat_delast_strain Derivative ∂W_p/∂ε_e [nwvec × ntvec]
         *                                   Plastic spin sensitivity to elastic strain
         *                                   Already includes full chain rule: ∂W_p/∂τ * ∂τ/∂ε_e
         * @param[in] A_e_M35 Elastic rotation coupling matrix [nwvec × ntvec = 3×5]
         *                    Relates elastic strain to elastic spin corrections
         * 
         * **Block location in Jacobian**:
         * - Rows: [ind_sub_r : ind_sub_r + nwvec - 1] (rotation equations)
         * - Columns: [0 : ntvec - 1] (elastic strain unknowns)
         * - Size: [3 × 5]
         * 
         * **Scaling notes**:
         * - dWp_hat_delast_strain already includes all necessary stress and kinetic derivatives
         * - Time step Δt multiplies plastic spin derivative (backward Euler)
         * - dt_ee_fac scales elastic corrections appropriately
         * - Additional scaling applied in problem's computeRJ for solver conditioning
         * 
         * **Higher-order terms significance**:
         * - Vanish for infinitesimal strains
         * - Important for: large elastic strains, stiff materials, high-rate loading
         * - Can be neglected in small-strain formulations (set elast_elast_factor = 0)
         * 
         * @note This is an off-diagonal block; diagonal blocks handled by separate methods
         * @note Symmetry: The transpose block ∂R_ε/∂ω computed by rotation problem's method
         * 
         * @see EvptnLatticeRotationProblem::get_deriv_elast_strain_wrt_omega for transpose block
         * @see elasticity_higher_order_terms() for A_e_M35 computation
         * @see M35_d_AAoB_dA() for A_edot_M35 calculation
         */
        template<size_t JAC_SIZE, size_t ind_sub_r>
        __ecmech_hdev__
        inline
        void get_deriv_omega_wrt_elast_strain(double* const jacobian,
                                            const double* const elast_dt_d5,
                                            const double elast_elast_factor,
                                            const double* const dWp_hat_delast_strain,
                                            const double* const A_e_M35) const
        {
        // d(B_xi)/d(elast_dev_press_vecs_f)
        //
        RAJA::View<double, RAJA::Layout<2>> jacob_re(jacobian, JAC_SIZE, JAC_SIZE);

        double A_edot_M35[ecmech::nwvec * ecmech::ntvec];
        M35_d_AAoB_dA(A_edot_M35, elast_dt_d5);

        double dt_ee_fac = m_dt * elast_elast_factor;

        for (size_t iWvec = 0; iWvec < ecmech::nwvec; ++iWvec) {
            for (size_t jTvec = 0; jTvec < ecmech::ntvec; ++jTvec) {
                int ijWT = ECMECH_NM_INDX(iWvec, jTvec, ecmech::nwvec, ecmech::ntvec);
                jacob_re(iWvec + ind_sub_r, jTvec) =
                    m_dt * dWp_hat_delast_strain[ijWT] - dt_ee_fac * (A_e_M35[ijWT] * m_inv_dt - A_edot_M35[ijWT]);
            }
        }         
        }

        /**
         * @brief Compute Jacobian contribution: ∂R_h/∂ε_e (hardening w.r.t. elastic strain).
         * 
         * Evaluates how hardening state residuals depend on elastic strain unknowns.
         * This coupling arises because hardening evolution depends on slip rates,
         * which in turn depend on resolved shear stress from elastic strain.
         * 
         * **Chain rule**:
         * ```
         * ∂R_h/∂ε_e = -Δt * (∂ḣ/∂γ̇)(∂γ̇/∂ε_e)
         * ```
         * where:
         * - ∂ḣ/∂γ̇: Hardening rate sensitivity to slip rates (from hardening law)
         * - ∂γ̇/∂ε_e: Slip rate sensitivity to elastic strain (from kinetics via stress)
         * - Minus sign: Residual form R_h = h^{n+1} - h^n - Δt*ḣ
         * 
         * @tparam JAC_SIZE Total Jacobian dimension
         * @tparam ind_sub_h Starting index for hardening unknowns
         * @tparam num_hard Number of hardening state variables
         * @tparam num_slip Number of slip systems
         * 
         * @param[in,out] jacobian Full Jacobian matrix [JAC_SIZE × JAC_SIZE]
         * @param[in] dhard_dgdot Derivative ∂ḣ/∂γ̇ [num_hard × num_slip]
         * @param[in] dgdot_delast_strain Derivative ∂γ̇/∂ε_e [num_slip × ntvec]
         * 
         * **Implementation**:
         * 1. Matrix multiply: dhdot_deps = (∂ḣ/∂γ̇)(∂γ̇/∂ε_e)ᵀ  [num_hard × ntvec]
         * 2. Scale by -Δt and insert into Jacobian
         * 3. Block location: Rows ind_sub_h:ind_sub_h+num_hard-1, Columns 0:ntvec-1
         */
        template<size_t JAC_SIZE, size_t ind_sub_h, size_t num_hard, size_t num_slip>
        __ecmech_hdev__
        inline
        void get_deriv_hardening_wrt_elast_strain(double* const jacobian,
                                                const double* const dhard_dgdot,
                                                const double* const dgdot_delast_strain) const
        {
        // d(B_h) / d(e)
        // jacob_he = dt * (dhdot/dgdot)(dgdot/de)
        // nh x ntvec matrix
        double dhdot_delast_strain[num_hard * ecmech::ntvec];
        vecsMABT<num_hard, ecmech::ntvec, num_slip>(dhdot_delast_strain, dhard_dgdot, dgdot_delast_strain);
        RAJA::View<double, RAJA::Layout<2>> jacob_he(jacobian, JAC_SIZE, JAC_SIZE);
        for (size_t iH = 0; iH < num_hard; ++iH) {
            for (size_t jE = 0; jE < ecmech::ntvec; ++jE) {
                // could also make dhdot_delast_strain into a RAJA view, but not really needed
                jacob_he(iH + ind_sub_h, jE) = -m_dt * dhdot_delast_strain[ECMECH_NM_INDX(iH, jE, num_hard, ecmech::ntvec) ];
            }
        }       
        }

        public:
        /**
         * @brief Index offset for elastic strain unknowns in combined systems.
         * Used when elastic strain is part of a larger unknown vector.
         * Value: 0 (elastic strain comes first)
         */
        static constexpr size_t m_ind_sub_elas = 0; // ntvec end_point
        /**
         * @brief Reference to thermoelastic constitutive model.
         * Provides stress evaluation and elastic stiffness.
         */
        const ThermoElastN& m_thermo_elast_n;
        /**
         * @brief Time step size [time units].
         */
        const double m_dt;
        /**
         * @brief Determinant of elastic volume ratio (J_e = det(F_e)).
         * Relates current to reference configuration volume.
         */
        const double m_det_v_e;
        /**
         * @brief Internal energy at end of time step [energy/volume].
         * Used in anisotropic Grüneisen contributions to stress.
         */
        const double m_energy_vol_ref;
        /**
         * @brief EOS pressure contribution [stress units].
         */
        const double m_pressure_EOS;
        /**
         * @brief Temperature [tempreture units].
         */
        const double m_tkelv;
        /**
         * @brief Elastic strain at beginning of step [dimensionless, 5 components].
         * Pointer to state in history array.
         */
        const double* const m_elast_d5_n;
        /**
         * @brief Inverse of time step (1/Δt) [1/time].
         * Pre-computed for efficiency in residual evaluation.
         */
        const double m_inv_dt;
        /**
         * @brief Inverse of elastic volume determinant (1/J_e).
         */
        const double m_inv_det_v_e;
        /**
         * @brief Volume scaling factor (a_vol = J_e^(1/3)) [dimensionless].
         * Converts between spatial and material tensor norms.
         */
        const double m_a_vol;
        /**
         * @brief Inverse volume scaling factor (1/a_vol).
         * Used in stress evaluation: σ = (1/a_vol) * K * ε_dev
         */
        const double m_inv_a_vol;
    };

    /**
     * @brief Subproblem formulation for crystal lattice rotation evolution.
     * 
     * This template class handles the implicit integration of crystal orientation
     * (represented as a unit quaternion) due to plastic spin from slip. It provides:
     * - Quaternion update formulation using exponential map
     * - Residual evaluation for rotation evolution equation  
     * - Jacobian computation for implicit integration
     * - Frame transformation utilities
     * 
     * **Governing equation**:
     * ```
     * R_ω = ξ - Δt * (W_sample - W_plastic)
     * ```
     * where:
     * - ξ: Incremental rotation vector (unknowns, 3 components)
     * - W_sample: Prescribed spin in sample frame (input)
     * - W_plastic: Plastic spin from slip (depends on stress via elastic strain)
     * - Quaternion updated via: Q^{n+1} = exp(ξ/2) ∘ Q^n
     * 
     * **Quaternion representation**:
     * - Unit quaternions Q = [q₀, q₁, q₂, q₃] with |Q| = 1
     * - Represents rotation from sample frame to crystal/lattice frame
     * - Exponential map: small rotation vector → quaternion increment
     * - Composition: Quaternion multiplication for sequential rotations
     * 
     * **Numerical approach**:
     * - Incremental rotation vector ξ as unknowns (not quaternion components)
     * - Avoids quaternion normalization constraints in Newton solver
     * - Exponential map provides smooth, singularity-free parametrization
     * - Scaling by r_scale for numerical conditioning
     * 
     * @tparam ind_sub_r Starting index for rotation unknowns in combined systems (default: ntvec)
     *                   When elastic strain and rotation solved together, rotation comes second
     * 
     * **Integration with larger problems**:
     * - Used standalone for rotation-only updates (rare)
     * - Typically coupled with EvptnLatticeStrainProblem in EvptnUpdstProblem
     * - Template parameter ind_sub_r enables flexible assembly in combined systems
     * 
     * @see EvptnUpdstProblem for coupled elastic strain + rotation integration
     * @see exponential map functions in ECMech_util.h (emap_to_quat, etc.)
     */
    template <size_t ind_sub_r=ecmech::ntvec>
    class EvptnLatticeRotationProblem {
        public:
        /**
         * @brief System dimension: number of rotation vector components.
         * Value: 3 (axial vector representation of rotation)
         */
        static constexpr size_t nDimSys = ecmech::nwvec;

        public:
        /**
         * @brief Constructor: Initialize lattice rotation subproblem.
         * 
         * Sets up rotation integration for one time step. Stores reference to
         * beginning-of-step orientation and time step size.
         * 
         * @param[in] dt Time step size [time units]
         * @param[in] xtal_ori_quat_n Beginning-of-step quaternion [4 components]
         *                             Represents rotation from sample to crystal frame
         * 
         * **Memory management**:
         * - xtal_ori_quat_n stored as const pointer (not copied)
         * - Caller responsible for ensuring validity during subproblem lifetime
         */
        __ecmech_hdev__
        EvptnLatticeRotationProblem(const double dt,
                                const double* const xtal_ori_quat_n)
        : m_dt(dt), m_xtal_ori_quat_n(xtal_ori_quat_n) {}
        /**
         * @brief Destructor (default).
         */
        ~EvptnLatticeRotationProblem() = default;

        /**
         * @brief Extract rotation state and compute frame transformation matrices.
         * 
         * This overload of get_rotation_state decodes the solution vector and additionally
         * computes the rotation matrices needed for frame transformations. This is more
         * efficient than separate calls when both rotation state and matrices are needed.
         * 
         * **Computation sequence**:
         * 1. Unscale solution vector: ξ = r_scale * x
         * 2. Convert to quaternion increment: A = exp(ξ/2)
         * 3. Compose quaternions: Q^{n+1} = A ∘ Q^n
         * 4. Generate 3×3 rotation matrix: R = quat_to_tensor(Q^{n+1})
         * 5. Generate 5×5 deviatoric rotation: R₅ = get_rot_mat_vecd(R)
         * 
         * **Output matrices**:
         * - xtal_rmat: 3×3 rotation matrix for vectors and skew tensors
         * - xtal_rot_mat5: 5×5 rotation matrix for deviatoric symmetric tensors
         * 
         * **Template parameter**:
         * @tparam ind_sub_r Starting index for rotation unknowns (default: ntvec = 5)
         * 
         * @param[out] delta_omega Incremental rotation vector ξ [radians, nwvec = 3]
         *                         Total rotation angle × axis over time step
         * @param[out] xtal_rmat Rotation matrix crystal→sample [3×3, row-major]
         *                       Use for rotating vectors (spin, positions)
         * @param[out] xtal_rot_mat5 Rotation matrix for 5-vectors [5×5, row-major]
         *                           Use for rotating deviatoric tensors (deformation rate)
         * @param[in] x Scaled solution vector [dimensionless, nwvec components]
         *              Offset by ind_sub_r if part of larger system
         * 
         * **Memory layout**:
         * - xtal_rmat: 9 elements, access via ECMECH_NN_INDX(i, j, 3)
         * - xtal_rot_mat5: 25 elements, access via ECMECH_NN_INDX(i, j, 5)
         * 
         * **Efficiency note**:
         * This combined operation avoids redundant quaternion operations when both
         * state extraction and matrix computation are needed in the same context.
         * 
         * @see get_rotation_state(double*, double*, const double*) for simpler overload
         * @see quat_to_tensor() for 3×3 rotation matrix generation
         * @see get_rot_mat_vecd() for 5×5 rotation matrix generation
         */
        __ecmech_hdev__
        inline
        void get_rotation_state(double* const delta_omega,
                                double* const xtal_rmat,
                                double* const xtal_rot_mat5,
                                const double* const x) const
        {
        vecsVxa<ecmech::nwvec>(delta_omega, ecmech::r_scale, &(x[ind_sub_r]));
        //
        // not done in EvpC :
        // CALL exp_map_cpvec(A, xi_f)
        // CALL get_c(c, A, C_n)
        //
        double xtal_ori_quat_delta[ecmech::qdim];
        double xtal_ori_quat_n1[ecmech::qdim];

        emap_to_quat(xtal_ori_quat_delta, delta_omega);
        get_c_quat(xtal_ori_quat_n1, xtal_ori_quat_delta, m_xtal_ori_quat_n);
        quat_to_tensor(xtal_rmat, xtal_ori_quat_n1);
        get_rot_mat_vecd(xtal_rot_mat5, xtal_rmat);
        }

        /**
         * @brief Reconstruct end-of-step orientation quaternion from solution vector.
         * 
         * Computes updated quaternion Q^{n+1} from incremental rotation:
         * 1. Unscale solution: ξ = r_scale * x
         * 2. Convert to quaternion increment: A = exp(ξ/2)
         * 3. Compose with beginning-of-step: Q^{n+1} = A ∘ Q^n
         * 
         * **Exponential map**:
         * - Maps rotation vector to quaternion via emap_to_quat()
         * - Factor of 1/2 relates quaternion space to rotation vector space
         * - Ensures Q remains unit quaternion (no normalization needed if exp accurate)
         * 
         * @param[out] xtal_ori_quat End-of-step orientation quaternion Q^{n+1} [4 components]
         * @param[in] x Scaled incremental rotation vector [dimensionless, 3 components]
         * 
         * @warning Not safe if quat and m_xtal_ori_quat_n point to same memory
         */
        __ecmech_hdev__
        inline
        void stateFromX(double* const xtal_ori_quat,
                    const double* const x) const
        {
        double delta_omega[ecmech::nwvec];
        double xtal_ori_quat_delta[ecmech::qdim];
        vecsVxa<ecmech::nwvec>(delta_omega, ecmech::r_scale, x);
        emap_to_quat(xtal_ori_quat_delta, delta_omega);
        get_c_quat(xtal_ori_quat, xtal_ori_quat_delta, m_xtal_ori_quat_n);
        }

        /**
         * @brief Compute incremental rotation from beginning and end quaternions.
         * 
         * This function extracts the rotation vector ξ that represents the incremental
         * rotation from orientation Q^n to Q^{n+1}. Useful for:
         * - Post-processing converged solutions
         * - Computing rotation residuals from state
         * - Analysis of rotation increments
         * 
         * **Mathematical operation**:
         * ```
         * 1. Compute relative quaternion: Q_Δ = Q^{n+1} ∘ (Q^n)*
         * 2. Convert to rotation vector: ξ = quat_to_emap(Q_Δ)
         * ```
         * 
         * **Quaternion conjugate**:
         * (Q^n)* represents the inverse rotation (conjugate quaternion)
         * 
         * **Rotation vector**:
         * - Magnitude: rotation angle in radians
         * - Direction: rotation axis
         * - Small rotations: ξ ≈ 2 * [q₁, q₂, q₃] for Q_Δ = [q₀, q₁, q₂, q₃]
         * 
         * @param[out] delta_omega Incremental rotation vector ξ [radians, nwvec = 3]
         * @param[in] xtal_ori_quat_n Beginning-of-step quaternion Q^n [qdim = 4]
         * @param[in] xtal_ori_quat_n1 End-of-step quaternion Q^{n+1} [qdim = 4]
         * 
         * @note Assumes quaternions are normalized (should be enforced during integration)
         * @note Result is unscaled (physical rotation vector, not scaled by r_scale)
         * 
         * @see quat_rel_rotation() for relative quaternion computation
         * @see quat_to_emap() for exponential map inversion
         * @see stateFromX() for the forward operation
         */
        __ecmech_hdev__
        inline
        void deltaOmegaFromState(double* const delta_omega,
                                 const double* const xtal_ori_quat_n,
                                 const double* const xtal_ori_quat_n1) const
        {
            double xtal_ori_quat_delta[ecmech::qdim] = {};

            quat_rel_rotation(xtal_ori_quat_delta, xtal_ori_quat_n1, xtal_ori_quat_n);
            quat_to_emap(delta_omega, xtal_ori_quat_delta);
        }

        /**
         * @brief Evaluate lattice rotation residual equation.
         * 
         * Computes the residual for the rotation evolution equation:
         * ```
         * R_ω = (ξ - Δt * (W_lat - W_p + ee_fac * W_e)) * rot_incr_scale_inv
         * ```
         * This should be zero at convergence.
         * 
         * **Physical interpretation**:
         * - ξ: Incremental lattice rotation (unknowns)
         * - W_lat: Lattice spin = prescribed spin in lattice frame
         * - W_p: Plastic spin from slip
         * - W_e: Elastic spin corrections (higher-order finite strain)
         * - Balance: ξ/Δt = W_lat - W_p + ee_fac * W_e
         * 
         * **Spin components**:
         * - W_lat (spin_vec_lat): Total prescribed spin rotated to crystal frame
         * - W_p (plastic_spin_vec): Computed from slip geometry and slip rates
         * - W_e (ee_spin_vec): Geometric corrections from elastic strain evolution
         * 
         * **Elastic spin factor**:
         * - ee_fac controls magnitude of elastic spin corrections
         * - Typically O(1) for standard formulations
         * - Can be set to 0 to neglect higher-order effects
         * 
         * **Scaling**:
         * - rot_incr_scale_inv: Conditioning factor ≈ 1/‖W_sample‖
         * - Makes residual O(1) for better Newton convergence
         * - Same scaling applied to Jacobian
         * 
         * **Template parameter**:
         * @tparam ind_sub_r Starting index for rotation DOFs in residual vector (default: ntvec = 5)
         * 
         * @param[out] residual Residual vector [dimensionless, full system size]
         *                      Only indices [ind_sub_r : ind_sub_r + nwvec - 1] modified
         * @param[in] rot_incr_scale_inv Inverse rotation scaling factor [time]
         * @param[in] ee_fac Elastic-elastic coupling factor [dimensionless]
         *                   Controls magnitude of elastic spin contribution
         * @param[in] delta_omega Incremental rotation ξ [radians, nwvec = 3]
         * @param[in] spin_vec_lat Lattice frame spin W_lat [rad/time, nwvec = 3]
         *                         Total prescribed spin rotated to crystal frame
         * @param[in] plastic_spin_vec Plastic spin W_p [rad/time, nwvec = 3]
         *                              From slip: W_p = Σ γ̇^α Q^α
         * @param[in] ee_spin_vec Elastic spin correction W_e [rad/time, nwvec = 3]
         *                        Higher-order geometric terms from elastic strain
         * 
         * **Residual indexing**:
         * - residual[ind_sub_r + 0]: W₁₂ component
         * - residual[ind_sub_r + 1]: W₂₃ component
         * - residual[ind_sub_r + 2]: W₁₃ component
         * (Indices for skew-symmetric tensor components in axial vector form)
         * 
         * **Frame considerations**:
         * All spin vectors must be in the same frame (crystal/lattice frame) for
         * consistent residual evaluation. spin_vec_lat obtained by rotating sample
         * frame spin to crystal frame.
         * 
         * **Convergence criterion**:
         * At convergence: ‖R_ω‖ < tolerance, meaning:
         * ```
         * ξ ≈ Δt * (W_lat - W_p + ee_fac * W_e)
         * ```
         * The incremental rotation balances prescribed and plastic spins.
         * 
         * @note Only modifies rotation DOF entries; elastic strain residual set elsewhere
         * @note ee_spin_vec typically small for moderate elastic strains
         * 
         * @see elasticity_higher_order_terms() for ee_spin_vec computation
         * @see get_slip_rate_terms() for plastic_spin_vec computation
         */
        __ecmech_hdev__
        inline
        void get_omega_residual(double* const residual,
                                const double rot_incr_scale_inv,
                                const double ee_fac,
                                const double* const delta_omega,
                                const double* const spin_vec_lat,
                                const double* const plastic_spin_vec,
                                const double* const ee_spin_vec) const
        {
        // RESIDUAL B_omega
        for (int iWvec = 0; iWvec < ecmech::nwvec; ++iWvec) {
            residual[ind_sub_r + iWvec] = rot_incr_scale_inv * // SCALING
                                        (delta_omega[iWvec] - m_dt * (spin_vec_lat[iWvec]
                                        - plastic_spin_vec[iWvec] + ee_fac * ee_spin_vec[iWvec]));
        }
        }

        /**
         * @brief Compute Jacobian: ∂R_ω/∂ω (rotation residual w.r.t. rotation).
         * 
         * Evaluates the diagonal block of the Jacobian corresponding to rotation DOFs.
         * This represents how the rotation residual changes with respect to the
         * rotation unknowns themselves.
         * 
         * **Mathematical formulation**:
         * ```
         * J_rr = ∂R_ω/∂ω = I - Δt * ∂W_lat/∂ω
         * ```
         * where:
         * - I: Identity matrix (from ξ term in residual)
         * - ∂W_lat/∂ω: Derivative of lattice frame spin w.r.t. rotation increment
         * 
         * **Derivative components**:
         * - dspin_samp_domega: Captures how rotation change affects frame transformation
         * - Includes geometric coupling: rotating spin vector to different frame
         * - Additional terms from rotation rate of frame itself
         * 
         * **Physical coupling**:
         * The derivative dspin_samp_domega captures:
         * 1. How changing rotation affects the rotation matrix
         * 2. How rotation matrix change affects transformed spin
         * 3. Rate of change of rotation itself
         * 
         * **Template parameter**:
         * @tparam JAC_SIZE Total Jacobian dimension (8 for full coupled problem)
         * 
         * @param[in,out] jacobian Full Jacobian matrix [JAC_SIZE × JAC_SIZE]
         *                         Modified to include rotation-rotation block
         * @param[in] dspin_samp_domega Derivative ∂W_lat/∂ω [nwvec × nwvec = 3×3]
         *                               Sensitivity of lattice frame spin to rotation
         * 
         * **Block location in Jacobian**:
         * - Rows: [ind_sub_r : ind_sub_r + nwvec - 1] (rotation equations)
         * - Columns: [ind_sub_r : ind_sub_r + nwvec - 1] (rotation unknowns)
         * - Size: [3 × 3]
         * 
         * **Matrix structure**:
         * Typically near-diagonal dominant:
         * - Diagonal: 1 - Δt * (rate coupling)
         * - Off-diagonal: -Δt * (geometric coupling)
         * 
         * **Scaling notes**:
         * - Raw derivative -Δt * dspin_samp_domega computed here
         * - Additional scaling (rotincr_scale_inv * r_scale) applied in problem's computeRJ
         * - Ensures Jacobian has O(1) entries for solver conditioning
         * 
         * @note This is the diagonal block for rotation; off-diagonal coupling to strain separate
         * @note Identity contribution ensures Jacobian is non-singular
         * 
         * @see eval_d_dxi_impl_quat() for dspin_samp_domega computation
         * @see get_deriv_elast_strain_wrt_omega() for rotation-strain coupling
         */
        template<size_t JAC_SIZE>
        __ecmech_hdev__
        inline
        void get_deriv_omega_wrt_omega(double* const jacobian,
                                    const double* const dspin_samp_domega) const
        {
        // d(B_xi)/d(xi_f)
        RAJA::View<double, RAJA::Layout<2>> jacob_rr(jacobian, JAC_SIZE, JAC_SIZE);
        for (int iWvec = 0; iWvec < ecmech::nwvec; ++iWvec) {
            for (int jWvec = 0; jWvec < ecmech::nwvec; ++jWvec) {
                int ijWW = ECMECH_NN_INDX(iWvec, jWvec, ecmech::nwvec);
                jacob_rr(iWvec + ind_sub_r, jWvec + ind_sub_r) = -m_dt * dspin_samp_domega[ijWW];
            }
            jacob_rr(iWvec + ind_sub_r, iWvec + ind_sub_r) += one;
        }
        }

        /**
         * @brief Compute Jacobian: ∂R_ε/∂ω (elastic strain residual w.r.t. rotation).
         * 
         * Evaluates the off-diagonal coupling between elastic strain residual and
         * rotation unknowns. This captures how changes in lattice rotation affect
         * the elastic strain evolution.
         * 
         * **Mathematical formulation**:
         * ```
         * J_er = ∂R_ε/∂ω = -∂D_sample_xtal/∂ω
         * ```
         * where D_sample_xtal is the deformation rate rotated to crystal frame.
         * 
         * **Physical coupling mechanism**:
         * - Rotation changes crystal frame orientation
         * - Different orientation means different D_sample in crystal coordinates
         * - Changed D_sample affects elastic strain evolution equation
         * 
         * **Derivative computation**:
         * ddef_rate_samp_domega captures:
         * ```
         * ∂(R^T * D_sample * R)/∂ω
         * ```
         * where R is the rotation matrix function of ω.
         * 
         * Note the negative sign: increased rotation typically decreases deformation
         * rate in crystal frame for most orientations.
         * 
         * **Template parameter**:
         * @tparam JAC_SIZE Total Jacobian dimension (8 for full coupled problem)
         * 
         * @param[in,out] jacobian Full Jacobian matrix [JAC_SIZE × JAC_SIZE]
         *                         Modified to include strain-rotation coupling block
         * @param[in] ddef_rate_samp_domega Derivative ∂D_xtal/∂ω [ntvec × nwvec = 5×3]
         *                                   How crystal frame deformation rate changes with rotation
         * 
         * **Block location in Jacobian**:
         * - Rows: [0 : ntvec - 1] (elastic strain equations)
         * - Columns: [ind_sub_r : ind_sub_r + nwvec - 1] (rotation unknowns)
         * - Size: [5 × 3]
         * 
         * **Coupling strength**:
         * - Magnitude depends on ‖D_sample‖ and current orientation
         * - Zero if D_sample = 0 (no prescribed deformation)
         * - Larger for orientations where rotation significantly reorients slip systems
         * 
         * **Scaling notes**:
         * - Raw derivative stored here (no time step or other scaling)
         * - Additional scaling applied in problem's computeRJ for conditioning
         * - Ensures off-diagonal coupling has appropriate magnitude relative to diagonal
         * 
         * **Symmetry consideration**:
         * This is the transpose coupling to get_deriv_omega_wrt_elast_strain:
         * - This method: How rotation affects elastic strain equation
         * - Transpose: How elastic strain affects rotation equation
         * 
         * @note Off-diagonal block; diagonal blocks handled by respective subproblems
         * @note Negative sign reflects kinematic constraint relationship
         * 
         * @see eval_d_dxi_impl_quat() for ddef_rate_samp_domega computation
         * @see EvptnLatticeStrainProblem::get_deriv_omega_wrt_elast_strain() for transpose
         */
        template<size_t JAC_SIZE>
        __ecmech_hdev__
        inline
        void get_deriv_elast_strain_wrt_omega(double* const jacobian,
                                            const double* const ddef_rate_samp_domega) const
        {
        // d(B_S)/d(xi_f)
        //
        // jacob_er = -dDsm_dxi(:,:)
        RAJA::View<double, RAJA::Layout<2> > jacob_er(jacobian, JAC_SIZE, JAC_SIZE);
        for (int jWvec = 0; jWvec<ecmech::nwvec; ++jWvec) {
            for (int iTvec = 0; iTvec<ecmech::ntvec; ++iTvec) {
                // could also make dDsm_dxi into a RAJA view, but not really needed
                const size_t ind = ECMECH_NM_INDX(iTvec, jWvec, ecmech::ntvec, ecmech::nwvec);
                jacob_er(iTvec, jWvec + ind_sub_r) = -ddef_rate_samp_domega[ind];
            }
        }
        }

        /**
         * @brief Compute Jacobian: ∂R_h/∂ω (hardening residual w.r.t. rotation).
         * 
         * This method exists for interface completeness but does nothing (empty function body).
         * 
         * @param[in,out] jacobian Jacobian matrix (not modified)
         * 
         * @note Function parameter commented out to avoid unused variable warning
         * @note Inline empty function optimized away by compiler (zero overhead)
         * 
         * @see EvptnLatticeStrainProblem::get_deriv_hardening_wrt_elast_strain() for actual coupling
         */
        __ecmech_hdev__
        inline
        void get_deriv_hardening_wrt_omega(double* const /* jacobian */){}

        public:
        /**
         * @brief Time step size [time units].
         */
        const double m_dt;
        /**
         * @brief Crystal orientation quaternion at beginning of step [unit quaternion, 4 components].
         * Pointer to state in history array.
         * Layout: [q₀, q₁, q₂, q₃] with |Q| = 1
         */
        const double* const m_xtal_ori_quat_n;
    };

}
}