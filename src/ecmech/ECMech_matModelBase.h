/**
 * @file ECMech_matModelBase.h
 * @brief Base class interface for all material constitutive models in ECMech.
 * 
 * This abstract base class defines the common interface that all ECMech material models
 * must implement. It provides a unified API for:
 * - Material parameter initialization and retrieval
 * - Constitutive response computation (stress update)
 * - History variable management
 * - Execution strategy selection (CPU/GPU)
 * - Model completion and validation
 * 
 * Derived classes implement specific constitutive models (e.g., crystal plasticity,
 * porosity models) by providing concrete implementations of the pure virtual methods.
 * 
 * Design pattern: Template Method pattern with Strategy pattern for execution
 * 
 * @see matModelEvptn for crystal plasticity implementations
 * @see EosModel for equation of state implementations
 * @see ThermoElastN for thermoelastic behavior
 */
// -*-c++-*-

#ifndef ECMech_matModelBase_include
#define ECMech_matModelBase_include

#include <string>
#include <sstream>
#include <vector>

#include "ECMech_core.h"

/**
 * @brief Utility macro for formatted vector output to ostringstream.
 * 
 * Outputs vector name followed by comma-separated values.
 * Example output: "params : 1.0, 2.5, 3.7"
 * 
 * @param aname Name of the vector (as string)
 * @param a Vector container (must have size() method and [] operator)
 */
#define DUMPVECOSS(aname, a) oss << aname << " : "; \
   for (unsigned int iThing = 0; iThing<a.size(); ++iThing) { if (iThing) { oss << ", "; } oss << a[iThing]; } oss << std::endl;

namespace ecmech {

   /**
    * @brief Abstract base class for material constitutive models.
    * 
    * matModelBase defines the interface that all material models must implement
    * for integration with finite element codes. It provides a unified interface
    * for material response evaluation, parameter management, and execution strategy
    * configuration.
    * 
    * Core responsibilities:
    * - Material response computation (stress update over time step)
    * - History variable management (internal state)
    * - Parameter initialization and retrieval
    * - Execution strategy control (CPU/GPU selection)
    * 
    * Derived classes must implement:
    * - initFromParams() - Initialize from parameter vectors
    * - getParams() - Retrieve current parameters
    * - getResponse() - Compute material response
    * - Query methods for dimensions and properties
    * 
    * @ingroup ECMech_material_models
    * 
    * @see ECMech_evptnWrap.h for template-based implementations
    * @see ECMech_cases.h for concrete model instantiations
    */
   class matModelBase
   {
      protected:
         /** @brief Flag indicating whether model initialization is complete */
         bool  m_complete;
         /** @brief Reference density (mass per unit reference volume) */
         double m_density0;
         /** @brief Average specific heat at constant volume for temperature evolution */
         double m_cvav;
         /** @brief Reference relative volume (typically 1.0 for unstressed state) */
         double m_rel_vol0;
         /** @brief Reference internal energy corresponding to reference state */
         double m_energy0;
         /** @brief Reference bulk modulus for equation of state scaling */
         double m_bulkRef;
         /** @brief Output verbosity level (0=silent, higher values=more verbose) */
         int m_outputLevel;
         /** @brief Execution strategy for parallel operations (CPU/GPU/OpenMP) */
         ecmech::ExecutionStrategy m_accel;

         /**
          * @brief Protected constructor for abstract base class.
          * 
          * Initializes all members to invalid/default values:
          * - m_complete = false (model not ready for use)
          * - Material properties = -1.0 (invalid, must be set)
          * - m_outputLevel = 0 (quiet)
          * - m_accel = CPU execution
          * 
          * Derived classes must call initFromParams() to set valid parameters.
          */
         __ecmech_host__
         matModelBase() :
            m_complete(false),
            m_density0(-1.0),
            m_cvav(-1.0),
            m_rel_vol0(-1.0),
            m_energy0(-1.0),
            m_bulkRef(-1.0),
            m_outputLevel(0),
            m_accel(ECM_EXEC_STRAT_CPU)
         {}

      public:
         /** @brief Virtual destructor for proper cleanup of derived classes */
         __ecmech_host__
         virtual ~matModelBase() {}

         /**
          * @brief Initialize material model from parameter vectors.
          * 
          * This method configures the material model from three parameter vectors containing
          * integer options, floating-point parameters, and string parameters. The exact
          * interpretation depends on the specific material model implementation.
          * 
          * Parameter organization (typical for crystal plasticity models):
          * - opts: Integer flags controlling model behavior (solver options, kinematic assumptions)
          * - pars: Material parameters in order: density, specific heat, tolerance, elastic constants,
          *         kinetic parameters, hardening parameters, EOS parameters
          * - strs: String identifiers for slip systems, kinetic laws, or model variants
          * 
          * The method performs:
          * 1. Parameter validation and bounds checking
          * 2. Derived quantity computation (e.g., elastic moduli from stiffness constants)
          * 3. Internal state initialization
          * 4. Memory allocation for model-specific data structures
          * 
          * @param[in] opts Integer options controlling model configuration
          * @param[in] pars Floating-point material parameters
          * @param[in] strs String parameters related to model instantiation
          * @param[in] call_back Optional callback pointer for sub-models (default: nullptr)
          * 
          * @note This method must be called before getResponse() or the model cannot function
          * @note Parameter order and number are model-specific; consult derived class documentation
          */
         __ecmech_host__
         virtual void initFromParams(const std::vector<int>& opts,
                                     const std::vector<double>& pars,
                                     const std::vector<std::string>& strs,
                                     void* call_back = nullptr) = 0;

         /**
          * @brief Retrieve current material model parameters.
          * 
          * Inverse operation of initFromParams(). Extracts all parameters needed to
          * reconstruct the current model state (excluding history variables).
          * 
          * @param[out] opts Integer options/flags (cleared and filled)
          * @param[out] pars Material parameters (cleared and filled)
          * @param[out] strs String parameters (cleared and filled)
          * 
          * **Purpose**:
          * - Serialization for restart files
          * - Parameter inspection/debugging
          * - Model comparison
          * 
          * **Guarantee**: After getParams(), calling initFromParams() with the same
          *                vectors should reconstruct an equivalent model
          * 
          * @note History variable state is NOT included (see getHistInfo() for initial values)
          */
         __ecmech_host__
         virtual void getParams(std::vector<int>& opts,
                                std::vector<double>& pars,
                                std::vector<std::string>& strs) const = 0;

         /**
          * @brief Log human-readable parameter information to output stream.
          * 
          * Provides formatted output of model parameters more suitable for human
          * inspection than the raw getParams() output. Default implementation
          * calls getParams() and formats using DUMPVECOSS macro.
          * 
          * @param[in,out] oss Output string stream to receive formatted text
          * 
          * **Default output format**:
          * ```
          * evptn constitutive model
          *   opts : <comma-separated integers>
          *   pars : <comma-separated doubles>
          *   strs : <comma-separated strings>
          * ```
          * 
          * **Override guidance**: Derived classes may override for more descriptive output
          *                       (e.g., label each parameter by name, group by category)
          * 
          * **Common usage**: Logging at simulation start, parameter verification
          */
         __ecmech_host__
         virtual void logParameters(std::ostringstream& oss) const {
            std::vector<int>         opts;
            std::vector<double>      pars;
            std::vector<std::string> strs;
            this->getParams(opts, pars, strs);
            oss << "evptn constitutive model" << std::endl;
            DUMPVECOSS("  opts", opts);
            DUMPVECOSS("  pars", pars);
            DUMPVECOSS("  strs", strs);
         }

         /**
          * @brief Request constitutive response for a batch of material points.
          * 
          * Central computational method: updates stress and history variables over a
          * finite time step given prescribed deformation. This is the main interface
          * between ECMech and host codes (FEM, MPM, etc.).
          * 
          * **Computational model**: Implicit time integration
          * - Beginning-of-step state provided as input
          * - End-of-step state computed and written to output
          * - Typically uses Newton-Raphson iteration internally
          * - Deviatoric stress and history variables: first-order (implicit, e.g., backward-Euler)
          * - Equation of state: may use third-order or first-order schemes
          * - Some implementations use closed-form analytic ODE solutions when available
          *   (assumes constant values of certain arguments over the time step)
          * 
          * **2D vs 3D Interface**:
          * - Interface is always for 3D deformation
          * - For 2D problems: some def_rate_d6vV and spin_vecV components will be zero
          * - Anisotropic materials: stress response can still be fully populated (non-zero in all components)
          * - State encoding: If stress is used to encode state for the material model
          *   (implementation-dependent), all stress components should be stored even in 2D
          * 
          * **Integration scheme compatibility**:
          * - Designed for finite time step integration (not instantaneous rates)
          * - NOT directly compatible with integration schemes requiring instantaneous
          *   rates of evolution for state variables and stress components
          * - Workaround: Possible to make multiple getResponse calls (e.g., at half-step
          *   and full-step dt) for multi-stage integration schemes
          * 
          * **Parallel execution**: Operates on batches of material points
          * - nPassed independent material points processed
          * - Execution strategy (CPU/GPU/OpenMP) set by setExecutionStrategy()
          * - Thread-safe: method is const, points are independent
          * 
          * @param[in] dt Time step size [s or simulation time units]
          *               - Must be positive
          *               - Convergence may require small steps for highly nonlinear behavior
          * 
          * @param[in] def_rate_d6vV Deviatoric deformation rate + volumetric rate [nsvp × nPassed]
          *            - Layout: [D_11', D_22', D_33', D_23, D_13, D_12, Dv] per point
          *            - Deviatoric: traceless (D_11' + D_22' + D_33' = 0)
          *            - Dv = volumetric rate = trace(D) / (dt × average_rel_vol)
          *            - Frame: Sample/global/laboratory
          * 
          * @param[in] spin_vecV Spin (vorticity) vector [ndim × nPassed]
          *            - Layout: [wxx, wyy, wzz] per point (axial vector form)
          *            - wxx = (L₃₂ - L₂₃)/2 = W₂₁ (rotation about x₁-axis)
          *            - wyy = (L₁₃ - L₃₁)/2 = W₀₂ (rotation about x₂-axis)
          *            - wzz = (L₂₁ - L₁₂)/2 = W₁₀ (rotation about x₃-axis)
          *            - W_ij = skew part of velocity gradient L
          *            - Frame: Sample/global
          * 
          * @param[in] rel_vol_ratiosV Relative volume ratios [nvr × nPassed]
          *            - Layout: [V/V₀, ΔV, dt, V_old/V₀] per point
          *            - V/V₀ = current relative volume (end of step)
          *            - ΔV = V_new - V_old (volume increment)
          *            - V_old/V₀ = relative volume at start of step
          * 
          * @param[in,out] internal_energyV Internal energy [ne × nPassed]
          *                - Layout: [e_total, e_cold, eQ, e_thermal, ..., e_melt]
          *                - Input: e_total (beginning-of-step), others zero
          *                - Output: e_total updated to end-of-step
          * 
          * @param[in,out] cauchy_stress_d6pV Cauchy stress [nsvp × nPassed]
          *                - Layout: [σ_11', σ_22', σ_33', σ_23, σ_13, σ_12, pressure]
          *                - Deviatoric: first 6 components traceless
          *                - Pressure: p = -⅓ trace(σ) (positive in compression)
          *                - Input: beginning-of-step state
          *                - Output: end-of-step state
          * 
          * @param[in,out] histV History (state) variables [numHist × nPassed]
          *                - Layout: Model-specific ordering (see getHistInfo())
          *                - Contains: slip resistances, dislocation densities, orientations, etc.
          *                - Input: beginning-of-step state
          *                - Output: end-of-step state
          *                - Size: numHist from getNumHist()
          * 
          * @param[out] tkelvV Temperature at end of step [nPassed]
          *             - Computed from equation of state model
          * 
          * @param[out] sddV Auxiliary output quantities [nsdd × nPassed]
          *            - Layout: Indexed by i_sdd_* constants
          *            - Example: i_sdd_gmod = shear modulus
          *            - Purpose: Derived quantities for output/analysis
          * 
          * @param[out] mtanSDV Tangent stiffness matrix [nsvec × nsvec × nPassed]
          *             - Layout: Row-major 6×6 symmetric matrix per point
          *             - Represents: ∂σ/∂ε (stress-strain tangent)
          *             - Purpose: Implicit FEM global tangent assembly
          *             - May be nullptr if not needed
          * 
          * @param[in] nPassed Number of material points to process
          *            - Must be ≥ 1
          *            - All arrays sized accordingly
          * 
          * **Indexing**: For length-X×nPassed arrays, indexing is fastest along X
          *               Example: stress[iPoint * nsvp + iComponent]
          * 
          * **Const correctness**: Method is const for thread safety
          *                       - No modification of material parameters
          *                       - Enables concurrent evaluation for different points
          * 
          * **Error handling**: Failed convergence may update histV with negative flags
          *                    (e.g., nFEval < 0 indicates failure)
          * 
          * @see getHistInfo() for history variable layout
          * @see setExecutionStrategy() for CPU/GPU selection
          */
         __ecmech_host__
         virtual void getResponseECM(const double & dt,
                                     const double * def_rate_d6vV,
                                     const double * spin_vecV,
                                     const double * rel_vol_ratiosV,
                                     double * internal_energyV,
                                     double * cauchy_stress_d6pV,
                                     double * histV,
                                     double * tkelvV,
                                     double * sddV,
                                     double * mtanSDV,
                                     const int & nPassed) const = 0;

         /**
          * @brief Get history variable metadata: names, initial values, output flags.
          * 
          * Provides complete information about the model's internal state variables
          * (history variables). Essential for:
          * - Allocating history arrays
          * - Initializing new material points
          * - Configuring output/visualization
          * - Post-processing and restart files
          * 
          * @param[out] histNames Names of history variables (cleared and filled)
          *             - Length: numHist
          *             - Example: "flowStr", "nFEval", "quat_0", "slipRes_1", etc.
          *             - Use for output file headers, visualization labels
          * 
          * @param[out] initVals Initial values for history variables (cleared and filled)
          *             - Length: numHist
          *             - Use to initialize histV array for new material points
          *             - Example: orientations=[1,0,0,0] (identity quaternion)
          *             - Example: slip resistances from initial hardening
          * 
          * @param[out] plot Flags indicating whether to output variable (cleared and filled)
          *             - Length: numHist
          *             - true: Variable should appear in standard visualization output
          *             - false: Internal/diagnostic variable (suppress from normal output)
          *             - Purpose: Reduce output file size by excluding low-level variables
          * 
          * @param[out] state Flags indicating whether variable is "state" (cleared and filled)
          *             - Length: numHist
          *             - true: True state variable (e.g., orientation, dislocation density)
          *             - false: Derived/diagnostic (e.g., number of function evals, rates)
          *             - Purpose: Distinguish essential restart data from auxiliary
          * 
          * **Typical history variable categories**:
          * 1. **Crystal orientation**: Quaternions or rotation matrix components
          * 2. **Slip resistances**: Current strength of each slip system
          * 3. **Dislocation densities**: Mobile and forest densities
          * 4. **Plastic strain**: Elastic strain in lattice frame
          * 5. **Diagnostic counters**: Function evaluations, convergence metrics
          * 
          * **Usage pattern**:
          * ```cpp
          * std::vector<std::string> names;
          * std::vector<double> initVals;
          * std::vector<bool> plot, state;
          * model->getHistInfo(names, initVals, plot, state);
          * 
          * // Allocate history array
          * double* histV = new double[model->getNumHist() * nPoints];
          * 
          * // Initialize all points
          * for (int ipt=0; ipt<nPoints; ++ipt) {
          *     std::copy(initVals.begin(), initVals.end(), 
          *               histV + ipt*model->getNumHist());
          * }
          * ```
          * 
          * @note Vector ordering matches histV array layout in getResponseECM()
          */
         virtual void getHistInfo(std::vector<std::string>& histNames,
                                  std::vector<double>& initVals,
                                  std::vector<bool>& plot,
                                  std::vector<bool>& state) const = 0;


         /**
          * @brief Query number of history variables per integration point.
          * @return Number of history variables
          */
         __ecmech_host__
         virtual int getNumHist( ) const = 0;

         /**
          * @brief Update array stride configuration for flexibility in memory layouts.
          * 
          * Allows customization of how data arrays are laid out in memory. Default
          * strides assume close-packed layout, but this method enables:
          * - Padding for alignment
          * - Interleaving with other data
          * - Compatibility with existing data structures
          * 
          * @param[in] strides Vector of stride values [nstride elements]
          *            - Must have size == ECMECH_NSTRIDE
          *            - Each element specifies stride for one array type:
          *              - strides[ISTRIDE_DEF_RATE]: def_rate_d6vV stride (≥ nsvp)
          *              - strides[ISTRIDE_SPIN_V]: spin_vecV stride (≥ ndim)
          *              - strides[ISTRIDE_VOL_RATIO]: rel_vol_ratiosV stride (≥ nvr)
          *              - strides[ISTRIDE_INT_ENG]: internal_energyV stride (≥ ne)
          *              - strides[ISTRIDE_STRESS]: cauchy_stress_d6pV stride (≥ nsvp)
          *              - strides[ISTRIDE_HISTORY]: histV stride (≥ numHist)
          *              - strides[ISTRIDE_TKELV]: tkelvV stride (≥ 1)
          *              - strides[ISTRIDE_SDD]: sddV stride (≥ nsdd)
          * 
          * **Constraints**:
          * - Must be called BEFORE complete()
          * - Each stride must be ≥ corresponding minimum size
          * - Calling after complete() triggers ECMECH_FAIL
          * 
          * **Default behavior**: If not called, strides = minimum required sizes
          * 
          * **Purpose**: Enable integration with codes using different memory layouts
          *             (e.g., AoS vs SoA, padding for vectorization)
          * 
          * **Example**:
          * ```cpp
          * std::vector<size_t> strides(ECMECH_NSTRIDE);
          * strides[ISTRIDE_STRESS] = 16;  // Pad stress to 16 for alignment
          * strides[ISTRIDE_HISTORY] = model->getNumHist() + 4;  // Extra space
          * model->updateStrides(strides);
          * ```
          * 
          * @warning Incorrect strides lead to memory corruption
          */
         __ecmech_host__
         virtual void updateStrides(std::vector<size_t> strides) = 0;

         /**
          * @brief Set execution strategy for parallel operations.
          * 
          * Configures how getResponseECM() executes batch computations:
          * - CPU: Serial execution on host
          * - GPU: CUDA/HIP execution on device
          * - OPENMP: OpenMP threading on host
          * 
          * @param[in] accel Desired execution strategy
          *            - ECM_EXEC_STRAT_CPU: Serial CPU execution
          *            - ECM_EXEC_STRAT_GPU: GPU execution (if compiled with support)
          *            - ECM_EXEC_STRAT_OPENMP: OpenMP parallel execution
          * 
          * **Default**: CPU execution if not called
          * 
          * **Runtime configuration**: Can be changed between getResponseECM() calls
          * 
          * **Availability**: GPU and OpenMP options depend on compile-time flags:
          * - GPU requires: RAJA_ENABLE_CUDA or RAJA_ENABLE_HIP
          * - OpenMP requires: RAJA_ENABLE_OPENMP
          * 
          * **Data location**: User responsible for ensuring data is accessible
          *                   (e.g., unified memory for GPU, host memory for CPU)
          * 
          * **Performance**: GPU strategy amortizes launch overhead over batch
          *                 - Efficient for nPassed ≥ 100-1000 points
          *                 - Small batches may be faster on CPU
          * 
          * @note Virtual to allow derived classes to validate or configure strategy
          */
         __ecmech_host__
         virtual void setExecutionStrategy(ecmech::ExecutionStrategy accel)  {
            m_accel = accel;
         }

         /**
          * @brief Get the reference density
          * @return the reference density unless value is < 0 incase an error is thrown
          */
         __ecmech_host__
         virtual double getRhoRef() const {
            if (m_density0 < 0.0) { // want to be able to call this before m_complete
               ECMECH_FAIL(__func__, "density0 does not appear to have been set");
            }
            return m_density0;
         }

         /**
          * @brief Set output verbosity level.
          * @param outputLevel Verbosity level (0 = silent)
          */
         __ecmech_host__
         void setOutputLevel(int outputLevel) { m_outputLevel = outputLevel; }


         /**
          * @brief Mark model as complete and ready for use.
          * 
          * Finalizes model setup after parameter initialization. Must be called
          * before getResponseECM(). Typically performs:
          * - Final validation of parameters
          * - Precomputation of derived quantities
          * - Allocation of any remaining storage
          * - Lock-down of configuration (no further parameter changes)
          * 
          * **Calling sequence**:
          * 1. Construct model
          * 2. Call initFromParams()
          * 3. (Optionally) call updateStrides()
          * 4. **Call complete()**
          * 5. Call getResponseECM() repeatedly
          * 
          * **Requirement**: Must be called exactly once before first getResponseECM()
          * 
          * **State change**: Sets m_complete = true
          * 
          * **Override pattern**: Derived classes can override to add additional finalization.
          * 
          * @note Virtual to allow derived class customization
          */
         __ecmech_host__
         virtual void complete() { m_complete = true; }

         /**
          * @brief Check if model initialization is complete.
          * @return true if initialized
          */
         __ecmech_host__
         virtual bool isComplete() { return m_complete; }
   }; // class matModelBase
} // ecmech namespace

#endif // ECMech_matModelBase_include
