/**
 * @file ECMech_elastic.h
 * @brief Thermoelastic constitutive models for cubic and hexagonal crystal symmetries.
 * 
 * This file provides thermoelastic response classes that compute:
 * - Kirchhoff stress from elastic lattice strain
 * - Cauchy stress from Kirchhoff stress and deformation
 * - Elastic stiffness contributions to tangent matrices
 * - Thermal stress effects via Grüneisen parameters
 * 
 * **Crystal symmetries supported**:
 * - **Cubic** (ThermoElastNCubic): FCC, BCC crystal structures
 *   - Parameters: c11, c12, c44 (3 independent elastic constants)
 *   - Symmetry: Isotropic in (111) family directions
 * 
 * - **Hexagonal** (ThermoElastNHexag): HCP crystal structures
 *   - Parameters: c11, c12, c13, c33, c44, g_vecd2 (5 elastic + 1 thermal)
 *   - Symmetry: Transverse isotropy about c-axis
 * 
 * **Physical formulation**:
 * - Elastic law: τ = C : ε^e (Kirchhoff stress from elastic strain)
 * - Deviatoric formulation: Uses 5-component deviatoric representation
 * - Thermal coupling: Grüneisen stress terms for pressure-energy coupling
 * - Frame: Crystal lattice frame (not sample frame)
 * 
 * **Key operations**:
 * 1. **eval()**: Compute Kirchhoff stress from elastic strain
 * 2. **getCauchy()**: Transform Kirchhoff to Cauchy stress
 * 3. **multDTDepsT()**: Apply elastic stiffness to matrix (tangent operations)
 * 4. **multCauchyDif()**: Compute Cauchy tangent stiffness contributions
 * 5. **getBulkMod()**: Extract bulk modulus
 * 6. **getGmod()**: Extract shear modulus (effective for anisotropic)
 * 
 * **Usage in crystal plasticity**:
 * - Compute stress from elastic lattice strain (after subtracting plastic strain)
 * - Provide elastic tangent for implicit Newton-Raphson solver
 * - Couple with slip geometry to get sample-frame response
 * 
 * **Naming convention**:
 * - "ThermoElastN" = Thermoelastic in Natural (crystal lattice) frame
 * - Contrast with sample-frame formulations
 * 
 * @see ThermoElastNCubic for cubic crystal elasticity
 * @see ThermoElastNHexag for hexagonal crystal elasticity
 * @see matModel for integration into crystal plasticity
 */

// -*-c++-*-
#pragma once

#include "ECMech_core.h"
#include "ECMech_util.h"
#include "ECMech_eosSimple.h"


namespace ecmech {
namespace evptn {

    /**
     * @brief Anisotropic thermoelastic constitutive model for cubic crystal symmetry.
     * 
     * ThermoElastNCubic implements the stress-strain relationship for materials with
     * cubic crystal symmetry (e.g., FCC, BCC metals) using a linear anisotropic elastic
     * formulation in the crystal reference frame. It computes lattice elastic stress
     * from logarithmic elastic strain tensors typical in finite deformation crystal
     * plasticity formulations.
     * 
     * Crystal symmetry and stiffness representation:
     * Cubic crystals have three independent elastic constants:
     * - C₁₁: Extensional stiffness along <100> directions
     * - C₁₂: Coupling between normal strains (Poisson-like effect)
     * - C₄₄: Shear stiffness along <110> type directions
     * 
     * The fourth-order stiffness tensor in Voigt notation:
     * ┌─                            ─┐
     * │ C₁₁  C₁₂  C₁₂  0    0    0   │
     * │ C₁₂  C₁₁  C₁₂  0    0    0   │
     * │ C₁₂  C₁₂  C₁₁  0    0    0   │
     * │ 0    0    0   C₄₄   0    0   │
     * │ 0    0    0    0   C₄₄   0   │
     * │ 0    0    0    0    0   C₄₄  │
     * └─                            ─┘
     * 
     * Strain measure and constitutive relation:
     * The model operates on deviatoric logarithmic elastic strain (e_dev) defined such that:
     *   tr(e) = 0 (traceless deviatoric part)
     *   σ = C : e_dev (linear elastic relation in crystal frame)
     * 
     * The stress is computed as Kirchhoff stress τ and then transformed to Cauchy stress σ
     * via appropriate volume scaling for the finite deformation framework.
     * 
     * Deviatoric representation (5-component vectors):
     * To enforce incompressibility and computational efficiency, the formulation uses
     * 5-component deviatoric vectors instead of full 6-component symmetric tensors:
     *   e_dev5 = [(e₁₁-e₂₂)/√2, (2e₃₃-e₁₁-e₂₂)/√6, √2·e₂₃, √2·e₁₃, √2·e₁₂]ᵀ
     * 
     * This representation automatically satisfies tr(e) = 0 and reduces degrees of freedom.
     * 
     * Material anisotropy and equivalent moduli:
     * For cubic materials, isotropic elastic moduli can be computed:
     * - Bulk modulus: K = (C₁₁ + 2C₁₂)/3
     * - Shear modulus: G = (2C₄₄ + C₁₁ - C₁₂)/5 (Voigt average)
     * - Anisotropy ratio: A = 2C₄₄/(C₁₁ - C₁₂)
     *   - A = 1: Isotropic material
     *   - A ≠ 1: Anisotropic (typical for single crystals)
     * 
     * Thermomechanical coupling:
     * While the class name includes "Thermo", temperature dependence of elastic constants
     * is not directly implemented in this version. Temperature effects enter through:
     * - Equation of state (pressure-volume-temperature relation)
     * - Reference configuration definition
     * Extensions could include temperature-dependent C_ij parameters.
     * 
     * Computational aspects:
     * - Diagonal stiffness representation (K_diag) enables efficient matrix-vector products
     * - Formulation supports both CPU and GPU execution (device-compatible functions)
     * - Tangent stiffness computation for Newton-Raphson implicit integration
     * - Transformation operators between crystal and sample reference frames
     * 
     * @ingroup ECMech_elasticity
     * 
     * @see ThermoElastNHexag for hexagonal symmetry materials
     * @see ECMech_evptn.h for crystal plasticity integration using this elasticity
     */
    class ThermoElastNCubic
    {
        public:
        /** @brief Number of parameters (elastic constants) required for cubic symmetry */
        static constexpr int nParams = 3;

        /**
         * @brief Default constructor creates uninitialized elasticity model.
         * 
         * Moduli are set to sentinel negative values indicating unconfigured state.
         * Must call setParams() before use.
         */
        __ecmech_hdev__
        inline ThermoElastNCubic() : m_bulk_modulus(-1.0), m_shear_modulus(-1.0) {}

        /** @brief Destructor */
        ~ThermoElastNCubic() = default;

        /**
         * @brief Constructor with parameter array initialization.
         * @param params Array of 3 elastic constants [C₁₁, C₁₂, C₄₄]
         */
         __ecmech_hdev__
         ThermoElastNCubic(const double* const params) {
            setParams(params);
         }

        /**
         * @brief Set parameters from vector (host-side convenience).
         * @param params Vector containing [C₁₁, C₁₂, C₄₄] in consistent units
         */
        __ecmech_host__
        inline void setParams(const std::vector<double> & params) {
            setParams(params.data());
        }

         /**
          * @brief Set cubic elastic constants and compute derived quantities.
          * 
          * Configures the anisotropic elastic stiffness for cubic symmetry and precomputes
          * quantities needed for efficient stress evaluation and tangent calculations.
          * 
          * Parameter requirements and physical constraints:
          * All elastic constants must be positive and satisfy thermodynamic stability:
          * - C₁₁ > 0 (extensional stability)
          * - C₄₄ > 0 (shear stability)
          * - C₁₁ > |C₁₂| (bulk stability)
          * - C₁₁ + 2C₁₂ > 0 (3D hydrostatic stability)
          * 
          * Parameter order:
          * 1. m_c11: Extensional stiffness C₁₁ [pressure units, typically GPa]
          *    - Controls resistance to extension along crystal axes
          *    - Largest of the three constants for most materials
          * 
          * 2. m_c12: Cross-coupling stiffness C₁₂ [pressure units]
          *    - Controls Poisson-like coupling between normal strains
          *    - Related to lateral contraction under uniaxial loading
          * 
          * 3. m_c44: Shear stiffness C₄₄ [pressure units]
          *    - Controls resistance to shear deformation
          *    - Independent from C₁₁ and C₁₂ due to crystal symmetry
          * 
          * Derived quantities computed:
          * - K_diag[5]: Diagonal of stiffness in 5-vector deviatoric representation
          *   - K_diag[0] = K_diag[1] = C₁₁ - C₁₂ (in-plane deviatoric stiffness)
          *   - K_diag[2] = K_diag[3] = K_diag[4] = 2C₄₄ (shear stiffnesses)
          * - Bulk modulus: K = (C₁₁ + 2C₁₂)/3
          * - Shear modulus: G = (2C₄₄ + C₁₁ - C₁₂)/5 (Voigt average)
          * 
          * Units consistency:
          * All three constants must use the same pressure units.
          * The code does not perform unit conversion; consistency is user's responsibility.
          * 
          * @param params Pointer to array of 3 doubles: [C₁₁, C₁₂, C₄₄]
          * 
          * @note In debug builds, validates that exactly nParams=3 values are read
          * @note Moduli validation (positive, thermodynamically stable) not enforced here
          */
         __ecmech_hdev__
         inline
         void setParams(const double* const params) {
            const double* parsIt = params;

            m_c11 = *parsIt; ++parsIt;
            m_c12 = *parsIt; ++parsIt;
            m_c44 = *parsIt; ++parsIt;
            //
#if defined(ECMECH_DEBUG)
            int iParam = parsIt - params;
            if (iParam != nParams) {
               ECMECH_FAIL(__func__, "iParam != nParams");
            }
#endif

            // Compute diagonal entries of stiffness in deviatoric 5-vector representation
            // For cubic symmetry, the deviatoric response is characterized by:
            // - In-plane deviatoric modes: (C₁₁ - C₁₂)
            // - Shear modes: 2C₄₄ (factor of 2 from tensor-to-vector conversion)
            m_K_diag[0] = m_c11 - m_c12;
            m_K_diag[1] = m_c11 - m_c12;
            m_K_diag[2] = two * m_c44;
            m_K_diag[3] = two * m_c44;
            m_K_diag[4] = two * m_c44;

            // Compute volumetric stiffness (trace of first 3x3 block)
            double K_vecds_s = m_c11 + two * m_c12;
            m_bulk_modulus = onethird * K_vecds_s;

            // Compute effective isotropic shear modulus using Voigt averaging
            // This provides a scalar approximation to the anisotropic shear response
            m_shear_modulus = (two * m_c11 - two * m_c12 + six * m_c44) * 0.2;
        }

         /**
          * @brief Retrieve elastic constants for serialization.
          * 
          * Extracts the three independent elastic constants in the order expected by
          * setParams(), enabling parameter inspection and model serialization.
          * 
          * @param[out] params Vector to receive [C₁₁, C₁₂, C₄₄] (not cleared; appended to)
          */
        __ecmech_host__
        inline void getParams(std::vector<double> & params
                                ) const {
    #ifdef ECMECH_DEBUG
            // do not clear params in case adding to an existing set
            int paramsStart = params.size();
    #endif

            params.push_back(m_c11);
            params.push_back(m_c12);
            params.push_back(m_c44);

    #ifdef ECMECH_DEBUG
            assert((params.size() - paramsStart) == nParams);
    #endif
        }

        /**
         * @brief Evaluate Kirchhoff stress from elastic strain in crystal frame.
         * 
         * Computes the Kirchhoff stress (volume-weighted Cauchy stress) from the deviatoric
         * elastic strain tensor using the cubic anisotropic elastic law.
         * The computation is performed in the crystal reference frame where the stiffness
         * tensor has its simple diagonal form.
         * 
         * Constitutive relation:
         *   τ = C : ε_dev
         * where:
         * - τ: Kirchhoff stress (6-component symmetric tensor)
         * - C: Fourth-order stiffness tensor (cubic symmetry)
         * - ε_dev: Deviatoric elastic strain (5-component traceless)
         * 
         * Volumetric contributions:
         * This method computes only deviatoric stress. The volumetric (pressure) part
         * is handled separately by the equation of state based on volume change and
         * thermal energy evolution.
         * 
         * Scaling factors:
         * - tkelv: Temperature - currently unused but reserved for temperature-dependent moduli
         * - pressure_EOS: Pressure from EOS - may affect elastic moduli in extended formulations
         * - energy_vol_ref: internal energy - provides thermodynamic consistency
         * 
         * The method applies volume scaling inv_a_vol to account for the relationship between
         * elastic stretch V^e and the spatial configuration. This ensures proper stress measures
         * for finite deformation kinematics.
         * 
         * @param[out] kirchoff_stress Computed Kirchhoff stress (6-component symmetric Voigt)
         * @param[in] elast_d5v Deviatoric elastic strain (5-component crystal frame)
         * @param[in] tkelv Temperature [K] (currently unused, reserved for future extensions)
         * @param[in] pressure_EOS Pressure from equation of state (currently unused)
         * @param[in] energy_vol_ref internal energy (currently unused)
         * @param[in] inv_a_vol Inverse of elastic volume scaling factor
         * 
         * @note Output kirchoff_stress has 6 independent components (symmetric Voigt notation)
         * @note The 7th component (pressure) is not set here; managed by EOS separately
         */
        __ecmech_hdev__
        inline
        void eval(double* const kirchoff,
                    const double* const elast_dev_press_vec,
                    double, // tkelv
                    double pressure_EOS,
                    double // energy_vol_ref
                    ) const {
            double ln_J = sqr3 * elast_dev_press_vec[iSvecS]; // vecds_s_to_trace
            double J = exp(ln_J);
            double kirchoff_pressure = -sqr3 * J * pressure_EOS;

            vecsVAdiagB<ntvec>(kirchoff, m_K_diag, elast_dev_press_vec);
            kirchoff[iSvecS] = kirchoff_pressure; // _K_vecds_s * elast_dev_press_vec(SVEC)
        }

        /**
         * @brief Multiply elastic stiffness diagonal by matrix (element-wise).
         * 
         * Applies the diagonal elastic stiffness tensor to a matrix A in a
         * computationally efficient manner. Used for:
         * - Computing elastic contributions to Jacobian matrices
         * - Applying elastic tangent in implicit solve
         * - Building tangent stiffness for FEM global assembly
         * 
         * **Operation**: P = D × A (element-wise along first index)
         * 
         * where:
         * - D = diag(K_diag × inv_a_vol) is the diagonal stiffness
         * - K_diag[i] are the elastic moduli for deviatoric components
         * - inv_a_vol = 1 / V^e (inverse elastic volume)
         * 
         * **Mathematical detail**:
         * ```
         * For cubic symmetry, the elastic stiffness dT/dε is diagonal:
         * 
         * dT/dε = diag([c11-c12, c11-c12, 2c44, 2c44, 2c44]) / V^e
         * 
         * This function computes:
         * P[i,j] = (K_diag[i] / V^e) × A[i,j]  for i ∈ [0,ntvec), j ∈ [0,p)
         * ```
         * 
         * **Indexing**:
         * - Deviatoric components indexed 0 to ntvec-1 (5 components)
         * - A and P are stored in row-major order
         * - Output P has same dimensions as input A [ntvec × p]
         * 
         * **Why diagonal?**: For cubic symmetry, elastic response in deviatoric
         * principal axes is uncoupled, so stiffness is diagonal.
         * 
         * **Usage context**: Called during Jacobian assembly in implicit solver
         * to add elastic stiffness contributions to tangent matrix.
         * 
         * @param[out] P Output matrix [ntvec × p]
         *             - P[i,j] = (K_diag[i] / V^e) × A[i,j]
         *             - Row-major storage: P[i*p + j]
         * 
         * @param[in] A Input matrix [ntvec × p]
         *            - Typically derivative matrix from implicit solve
         *            - Row-major storage: A[i*p + j]
         * 
         * @param[in] inv_a_vol Inverse elastic volume (1 / V^e)
         *            - V^e = det(F^e) where F^e is elastic deformation gradient
         *            - Accounts for volume change in elastic configuration
         * 
         * @param[in] p Number of columns in A and P
         *            - Typically ntvec (5) for deviatoric tangent
         *            - Can be other sizes for different derivative operations
         * 
         * **Performance**: O(ntvec × p) with diagonal stiffness optimization
         * 
         * @note For cubic symmetry, dT/dε[iSvecS,:] = 0 (no deviatoric-pressure coupling)
         * @note Method signature identical for cubic and hexagonal (same interface)
         */
        __ecmech_hdev__
        inline
        void multDTDepsT(double* const P, // ntvec*p
                            const double* const A, // ntvec*p
                            double inv_a_vol,
                            int p) const {
            for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                double dTdepsThis = m_K_diag[iTvec] * inv_a_vol;
                for (int iP = 0; iP < p; ++iP) {
                    int ii = ECMECH_NM_INDX(iTvec, iP, ecmech::ntvec, p);
                    P[ii] = dTdepsThis * A[ii];
                }
            }
        }

         /**
          * @brief Convert Kirchhoff stress to Cauchy stress with volume scaling.
          * 
          * Transforms the Kirchhoff stress τ to Cauchy stress σ using the determinant
          * of the elastic deformation:
          *   σ = τ / det(F^e) = τ / det(V^e)
          * 
          * This scaling accounts for the difference between stress measures:
          * - Kirchhoff stress: Work-conjugate to logarithmic strain rate
          * - Cauchy stress: True physical stress (force per current area)
          * 
          * @param[out] cauchy_stress Cauchy stress (6-component symmetric)
          * @param[in] kirchoff_stress Kirchhoff stress (6-component symmetric)
          * @param[in] inv_det_v_e Inverse of elastic deformation determinant 1/det(V^e)
          */
        __ecmech_hdev__
        inline
        void getCauchy(double* const cauchy_xtal,
                        const double* const kirchoff,
                        double inv_det_v_e) const
        {
            for (int iSvec = 0; iSvec < ecmech::nsvec; ++iSvec) {
                cauchy_xtal[iSvec] = inv_det_v_e * kirchoff[iSvec];
            }
        }

        /**
         * @brief Derivative of deviatoric Kirchhoff stress w.r.t. spherical elastic strain.
         *
         * Fills dT_deps_sph[0:ntvec) with ∂τ'_i/∂ε_s, where ε_s is the spherical
         * (vecds) component of the elastic strain passed to eval(). For cubic
         * symmetry there is no deviatoric-volumetric elastic coupling, so the
         * vector is identically zero.
         *
         * @param[out] dT_deps_sph Coupling vector [ntvec], zero-filled here
         * @return false, indicating no coupling (callers may skip dependent work)
         */
        __ecmech_hdev__
        inline
        bool getDTDepsSph(double* const dT_deps_sph) const {
            for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                dT_deps_sph[iTvec] = 0.0;
            }
            return false;
        }

        /**
         * @brief Compute tangent stiffness contribution for implicit integration.
         *
         * Calculates the material tangent operator needed for Newton-Raphson iteration
         * in implicit finite element analysis. For the elastic part, this is:
         *   ∂σ/∂ε_rate = (1/det(V^e)) · (1/a_vol) · C
         *
         * The tangent includes scaling factors accounting for:
         * - Finite deformation kinematics (inv_det_v_e)
         * - Volume normalization (inv_a_vol)
         * - Transformation to appropriate rate form
         *
         * A holds the converged sensitivities of the deviatoric lattice strain with
         * respect to all nsvec vecds components of the sample-frame deformation
         * rate: A[i,j] = ∂(elast_d5)_i / ∂(def_rate_vecds)_j, i < ntvec, j < nsvec.
         * The final (spherical) column of A comes from the volumetric right-hand
         * side of the tangent solve; for cubic symmetry that column is zero because
         * there is no deviatoric-volumetric elastic coupling.
         *
         * @tparam N Row dimension of input matrix A (ntvec)
         * @tparam M Column dimension of input matrix A (nsvec)
         *
         * @param[out] M6 Tangent stiffness matrix (nsvec×nsvec), crystal-frame rows
         * @param[in] A Strain sensitivity matrix (N×M), see above
         * @param[in] deps_sph_ddef_rate_sph ∂ε_s/∂(def_rate_vecds)_s over the step
         *            (Δt at the current level of EOS decoupling); unused for cubic
         * @param[in] inv_det_v_e Inverse elastic deformation determinant
         * @param[in] inv_a_vol Inverse volume scaling factor
         *
         * @note The volumetric (S,S) pressure stiffness comes from the equation of
         *       state and is added separately by the caller
         */
        template<size_t N=ecmech::ntvec, size_t M=ecmech::nsvec>
        __ecmech_hdev__
        inline
        void multCauchyDif(double* const M6,
                            const double* const A,
                            double, // deps_sph_ddef_rate_sph -- no coupling for cubic
                            double inv_det_v_e,
                            double inv_a_vol
                            ) const {
            // CALL vecds_s_to_trace(tr_ln_V, s_meas%elast_dev_press_vec(SVEC))
            // det_v_e = DEXP(tr_ln_V)
            // inv_det_v_e = one / det_v_e

            // dsigC_de(:,:) = inv_det_v_e * s_meas%dT_deps(:,:)
            // for cubic, dT_deps is diag(K_diag * a_V%ri) (symmetric) ; dT_deps[iSvecS,:] = 0
            // M65_ij = dd_ii A_ij
            // Compute scaled tangent: M6 = (inv_det_v_e * inv_a_vol * K_diag) ⊙ A
            // where ⊙ represents element-wise multiplication in the tangent structure
            for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                double vFact = inv_det_v_e * inv_a_vol * m_K_diag[iTvec];
                for (int jSvec = 0; jSvec < ecmech::nsvec; ++jSvec) {
                    M6[ECMECH_NN_INDX(iTvec, jSvec, ecmech::nsvec)] = vFact * A[ECMECH_NM_INDX(iTvec, jSvec, N, M)];
                }
            }

            // no deviatoric-volumetric elastic coupling for cubic: the pressure row
            // is entirely the equation-of-state's business
            for (int jSvec = 0; jSvec < ecmech::nsvec; ++jSvec) {
                M6[ECMECH_NN_INDX(iSvecS, jSvec, ecmech::nsvec)] = 0.0;
            }
        }

        /**
         * @brief Query bulk modulus with initialization check.
         * @return Bulk modulus K = (C₁₁ + 2C₁₂)/3
         * @throws Error if bulk modulus is negative (indicates uninitialized model)
         */
        __ecmech_hdev__
        inline
        double getBulkMod( ) const {
            if (m_bulk_modulus <= 0.0) {
                ECMECH_FAIL(__func__, "bulk modulus negative -- not initialized?");
            }
            return m_bulk_modulus;
        }

        /**
         * @brief Query effective shear modulus with initialization check.
         * 
         * Returns the Voigt-averaged isotropic shear modulus approximation for the
         * anisotropic cubic material. Unused parameters are reserved for potential
         * temperature or pressure dependence in future extensions.
         * 
         * @param tkelv Temperature (currently unused)
         * @param pressure_EOS Pressure (currently unused)
         * @param energy_vol_ref Volumetric energy (currently unused)
         * @return Shear modulus G = (2C₄₄ + C₁₁ - C₁₂)/5
         * @throws Error if shear modulus is negative (indicates uninitialized model)
         */
        __ecmech_hdev__
        inline
        double getGmod(double, // tkelv
                        double, // pressure_EOS
                        double // energy_vol_ref
                        ) const {
            if (m_shear_modulus <= 0.0) {
                ECMECH_FAIL(__func__, "effective shear modulus negative -- not initialized?");
            }
            return m_shear_modulus;
        }

        private:
        /** @brief Elastic constants for cubic symmetry */
        double m_c11, m_c12, m_c44;
        /** @brief Diagonal stiffness entries in deviatoric 5-vector representation */
        double m_K_diag[ecmech::ntvec];
        /** @brief Bulk modulus (volumetric stiffness) K = (C₁₁ + 2C₁₂)/3 */
        double m_bulk_modulus;
        /** @brief Effective isotropic shear modulus (Voigt average) */
        double m_shear_modulus;
    };

    /**
     * @brief Anisotropic thermoelastic constitutive model for hexagonal crystal symmetry.
     * 
     * ThermoElastNHexag implements stress-strain relations for hexagonal close-packed (HCP)
     * materials like titanium, magnesium, zinc, and zirconium. HCP crystals exhibit
     * transverse isotropy with the unique c-axis perpendicular to the basal plane.
     * 
     * Crystal symmetry and stiffness representation:
     * Hexagonal crystals have five independent elastic constants:
     * - C₁₁: In-plane extensional stiffness (basal plane)
     * - C₁₂: In-plane Poisson coupling (basal plane)
     * - C₁₃: Out-of-plane Poisson coupling (basal-prism interaction)
     * - C₃₃: Out-of-plane extensional stiffness (c-axis direction)
     * - C₄₄: Out-of-plane shear stiffness (prismatic shear)
     * 
     * Note: In-plane shear stiffness is constrained by: C₆₆ = (C₁₁ - C₁₂)/2
     * 
     * The stiffness tensor in Voigt notation (crystal coordinates):
     * ┌─                                ─┐
     * │ C₁₁   C₁₂   C₁₃    0     0     0 │
     * │ C₁₂   C₁₁   C₁₃    0     0     0 │
     * │ C₁₃   C₁₃   C₃₃    0     0     0 │
     * │  0     0     0    C₄₄    0     0 │
     * │  0     0     0     0    C₄₄    0 │
     * │  0     0     0     0     0   C₆₆ │ where C₆₆ = (C₁₁-C₁₂)/2
     * └─                                ─┘
     * 
     * Coordinate system convention:
     * - x₁, x₂: Basal plane directions (equivalent by 6-fold symmetry)
     * - x₃: c-axis direction (unique axis of hexagonal symmetry)
     * - The c-axis typically aligns with [0001] crystallographic direction
     * 
     * Material anisotropy characteristics:
     * - c/a ratio: Atomic structure parameter affecting elastic anisotropy
     *   - Ideal HCP: c/a = √(8/3) ≈ 1.633
     * - Basal-to-prism stiffness ratio: C₃₃/C₁₁
     *   - > 1: Stiffer along c-axis
     *   - < 1: Softer along c-axis
     * - Zener anisotropy analog: 2C₄₄/(C₁₁ + C₃₃ - 2C₁₃)
     * 
     * Gruneisen parameter anisotropy:
     * HCP materials exhibit anisotropic thermal expansion, characterized by:
     * - Γ_a: Gruneisen parameter in basal plane directions
     * - Γ_c: Gruneisen parameter along c-axis
     * Stored efficiently as single parameter m_g_vecd2 representing the difference.
     * 
     * The Gruneisen tensor is diagonal: Γ = diag(Γ_a, Γ_a, Γ_c)
     * In deviatoric 5-vector form, only one component is non-zero:
     *   g_vecd2 = 2(Γ_c - Γ_a)/√6
     * 
     * @ingroup ECMech_elasticity
     * 
     * @see ThermoElastNCubic for cubic symmetry materials
     * @see ECMech_slipgeom_hcp.h for hexagonal slip geometry definitions
     */
    class ThermoElastNHexag
    {
        public:
        /** @brief Number of parameters: 5 elastic constants + 1 Gruneisen parameter */
        static constexpr int nParams = 6;

        /**
         * @brief Default constructor creates uninitialized elasticity model.
         */
        __ecmech_hdev__
        inline ThermoElastNHexag() : m_bulk_modulus(-1.0), m_shear_modulus(-1.0) {}

        /** @brief Destructor */
        ~ThermoElastNHexag() = default;

        /**
         * @brief Constructor with parameter array initialization.
         * @param params Array of 6 values: [C₁₁, C₁₂, C₁₃, C₃₃, C₄₄, g_vecd2]
         */
         __ecmech_hdev__
         ThermoElastNHexag(const double* const params) {
            setParams(params);
         }

        /**
         * @brief Set parameters from vector (host-side convenience).
         * @param params Vector of 6 parameters
         */
        __ecmech_host__
        inline void setParams(const std::vector<double> & params) {
            setParams(params.data());
        }

         /**
          * @brief Set hexagonal elastic constants and compute derived quantities.
          * 
          * Parameters and physical meaning:
          * 1. C₁₁: In-plane (basal) extensional stiffness
          * 2. C₁₂: In-plane Poisson coupling  
          * 3. C₁₃: Basal-prism coupling
          * 4. C₃₃: Out-of-plane (c-axis) extensional stiffness
          * 5. C₄₄: Out-of-plane shear stiffness
          * 6. g_vecd2: Anisotropic Gruneisen parameter (thermal property)
          * 
          * @param params Array of 6 doubles in the order listed above
          */
         __ecmech_hdev__
         inline
         void setParams(const double* const params) {
            const double* parsIt = params;

            m_c11 = *parsIt; ++parsIt;
            m_c12 = *parsIt; ++parsIt;
            m_c13 = *parsIt; ++parsIt;
            m_c33 = *parsIt; ++parsIt;
            m_c44 = *parsIt; ++parsIt;
            //
            m_g_vecd2 = *parsIt; ++parsIt;
            //
#if defined(ECMECH_DEBUG)
            int iParam = parsIt - params;
            if (iParam != nParams) {
               ECMECH_FAIL(__func__, "iParam != nParams");
            }
#endif
            // Diagonal stiffness entries for deviatoric 5-vector representation
            // These account for hexagonal symmetry and deviatoric decomposition
            m_K_diag[0] = m_c11 - m_c12;  // In-plane deviatoric mode 1
            m_K_diag[1] = m_c11 * onethird + m_c12 * onethird - fourthirds * m_c13 + twothird * m_c33;  // Hexagonal deviatoric mode
            m_K_diag[2] = m_c11 - m_c12;  // In-plane deviatoric mode 2
            m_K_diag[3] = two * m_c44;    // Out-of-plane shear mode 1
            m_K_diag[4] = two * m_c44;    // Out-of-plane shear mode 2
            // Volumetric stiffness
            double K_vecds_s = twothird * m_c11 + twothird * m_c12 + fourthirds * m_c13 + m_c33 * onethird;
            // Off-diagonal coupling term for hexagonal symmetry
            m_K_sdax3 = sqr2 * (-m_c11 - m_c12 + m_c13 + m_c33) * onethird;
            m_bulk_modulus = onethird * K_vecds_s;
            //
            // m_shear_modulus below ignores the m_K_sdax3 contribution, but it is just meant to be approximate anyway
            m_shear_modulus = 0.5 * 0.2 * vecsssum<ecmech::ntvec>(m_K_diag); // 0.5 * (average of m_K_diag entries)
        }

         /**
          * @brief Retrieve elastic constants for serialization.
          * @param[out] params Vector to receive [C₁₁, C₁₂, C₁₃, C₃₃, C₄₄, g_vecd2]
          */
        __ecmech_host__
        inline void getParams(std::vector<double> & params
                                ) const {
    #ifdef ECMECH_DEBUG
            // do not clear params in case adding to an existing set
            int paramsStart = params.size();
    #endif

            params.push_back(m_c11);
            params.push_back(m_c12);
            params.push_back(m_c13);
            params.push_back(m_c33);
            params.push_back(m_c44);
            //
            params.push_back(m_g_vecd2);

    #ifdef ECMECH_DEBUG
            assert((params.size() - paramsStart) == nParams);
    #endif
        }

        /**
         * @brief Evaluate Kirchhoff stress from elastic strain (hexagonal symmetry).
         * 
         * Similar to cubic case but accounting for hexagonal anisotropy with its
         * additional independent elastic constants and off-diagonal coupling.
         * 
         * @param[out] kirchoff_stress Kirchhoff stress (6-component)
         * @param[in] elast_d5v Deviatoric elastic strain (5-component)
         * @param[in] tkelv Temperature [K]
         * @param[in] pressure_EOS Pressure from EOS
         * @param[in] energy_vol_ref internal energy
         * @param[in] inv_a_vol Inverse volume scaling
         */
        __ecmech_hdev__
        inline
        void eval(double* const kirchoff,
                    const double* const elast_dev_press_vec,
                    double, // tkelv
                    double pressure_EOS,
                    double energy_vol_ref
                    ) const {
            double ln_J = sqr3 * elast_dev_press_vec[iSvecS]; // vecds_s_to_trace
            double J = exp(ln_J);
            double kirchoff_pressure = -sqr3 * J * pressure_EOS;

            vecsVAdiagB<ntvec>(kirchoff, m_K_diag, elast_dev_press_vec);
            kirchoff[iSvecS] = kirchoff_pressure; // _K_vecds_s * elast_dev_press_vec(SVEC)

            kirchoff[iTvecHex] += m_K_sdax3 * elast_dev_press_vec[iSvecS];
            kirchoff[iSvecS] += m_K_sdax3 * elast_dev_press_vec[iTvecHex];

            // anisotropic Gruneisen contribution; pressure part of Gruneisen tensor contribution should already be in pressure_EOS
            // CALL eos_eval_e_Csdev(Cauchy_eos_vecd, energy_vol_ref, J, &
            // & i_eos_model, eos_const)
            // -(Gamma' + a' * mu) * energy_vol_ref // but do not do a'*mu part
            // Cauchy_eos_vecd(:) = -eos_const(4:8) * energy_vol_ref
            // kirchoff(1:TVEC) = kirchoff(1:TVEC) + J * Cauchy_eos_vecd(:)
            kirchoff[iTvecHex] += J * (-m_g_vecd2 * energy_vol_ref);
        }

        /**
         * @brief Multiply elastic stiffness diagonal by matrix (element-wise).
         * 
         * Applies the diagonal elastic stiffness tensor to a matrix A. For hexagonal
         * symmetry, this is IDENTICAL to the cubic case because the additional coupling
         * term m_K_sdax3 does not enter this diagonal multiplication.
         * 
         * **Operation**: P = D × A (element-wise along first index)
         * 
         * where:
         * - D = diag(K_diag × inv_a_vol) is the diagonal stiffness
         * - K_diag[i] are the elastic moduli for deviatoric components
         * - inv_a_vol = 1 / V^e (inverse elastic volume)
         * 
         * **Hexagonal elasticity note**:
         * The hexagonal stiffness has an additional off-diagonal term m_K_sdax3
         * that couples component iTvecHex (index 1) with pressure. However, this
         * coupling does NOT appear in the diagonal operation multDTDepsT - it only
         * appears in multCauchyDif where the full stiffness tensor is needed.
         * 
         * **Mathematical detail**:
         * ```
         * For hexagonal symmetry, the diagonal part of dT/dε is:
         * 
         * K_diag[0] = c11 - c12
         * K_diag[1] = (c11 + c12)/3 - 4c13/3 + 2c33/3  (basal plane contribution)
         * K_diag[2] = c11 - c12
         * K_diag[3] = 2c44
         * K_diag[4] = 2c44
         * 
         * Off-diagonal coupling (not used here):
         * K_sdax3 = √2 × (-c11 - c12 + c13 + c33) / 3
         * 
         * This function applies ONLY the diagonal:
         * P[i,j] = (K_diag[i] / V^e) × A[i,j]  for i ∈ [0,ntvec), j ∈ [0,p)
         * ```
         * 
         * **Hexagonal vs Cubic**:
         * - Cubic: 3 independent constants (c11, c12, c44)
         * - Hexagonal: 5 independent constants (c11, c12, c13, c33, c44)
         * - Both have diagonal K_diag, but hexagonal has additional off-diagonal terms
         * - This function uses ONLY K_diag, so implementation is identical
         * 
         * **Usage context**: Same as cubic - Jacobian assembly for implicit solver.
         * The off-diagonal m_K_sdax3 coupling is handled separately in multCauchyDif.
         * 
         * @param[out] P Output matrix [ntvec × p]
         *             - P[i,j] = (K_diag[i] / V^e) × A[i,j]
         *             - Row-major storage
         * 
         * @param[in] A Input matrix [ntvec × p]
         *            - Input derivative matrix
         *            - Row-major storage
         * 
         * @param[in] inv_a_vol Inverse elastic volume (1 / V^e)
         *            - V^e = det(F^e)
         * 
         * @param[in] p Number of columns in A and P
         *            - Typically ntvec (5) for deviatoric operations
         * 
         * **Implementation note**: The comment in the original code states:
         * "multDTDepsT ends up looking the same as in the cubic case because
         * m_K_sdax3 does not enter" - this function applies only the diagonal
         * stiffness, not the full coupling tensor.
         * 
         * @note Full hexagonal stiffness (with m_K_sdax3) used in multCauchyDif
         * @note iTvecHex = 1 (basal plane component) has special anisotropy
         */
        __ecmech_hdev__
        inline
        void multDTDepsT(double* const P, // ntvec*p
                            const double* const A, // ntvec*p
                            double inv_a_vol,
                            int p) const {
            for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                double dTdepsThis = m_K_diag[iTvec] * inv_a_vol;
                for (int iP = 0; iP < p; ++iP) {
                    int ii = ECMECH_NM_INDX(iTvec, iP, ecmech::ntvec, p);
                    P[ii] = dTdepsThis * A[ii];
                }
            }
        }

         /**
          * @brief Convert Kirchhoff to Cauchy stress (hexagonal).
          * @param[out] cauchy_stress Cauchy stress
          * @param[in] kirchoff_stress Kirchhoff stress
          * @param[in] inv_det_v_e Inverse elastic deformation determinant
          */
        __ecmech_hdev__
        inline
        void getCauchy(double* const cauchy_xtal,
                        const double* const kirchoff,
                        double inv_det_v_e) const
        {
            for (int iSvec = 0; iSvec < ecmech::nsvec; ++iSvec) {
                cauchy_xtal[iSvec] = inv_det_v_e * kirchoff[iSvec];
            }
        }

        /**
         * @brief Derivative of deviatoric Kirchhoff stress w.r.t. spherical elastic strain.
         *
         * Fills dT_deps_sph[0:ntvec) with ∂τ'_i/∂ε_s, where ε_s is the spherical
         * (vecds) component of the elastic strain passed to eval(). For hexagonal
         * symmetry the coupling is a single entry:
         *   ∂τ'_{iTvecHex}/∂ε_s = K_sdax3 = √2(-c₁₁ - c₁₂ + c₁₃ + c₃₃)/3 .
         * Note eval() applies no a_vol scaling to the spherical strain slot, so no
         * inv_a_vol factor belongs here either.
         *
         * @param[out] dT_deps_sph Coupling vector [ntvec]
         * @return true, indicating coupling is present
         */
        __ecmech_hdev__
        inline
        bool getDTDepsSph(double* const dT_deps_sph) const {
            for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                dT_deps_sph[iTvec] = 0.0;
            }
            dT_deps_sph[iTvecHex] = m_K_sdax3;
            return true;
        }

         /**
          * @brief Compute tangent stiffness for hexagonal symmetry.
          *
          * Similar to cubic case but includes the off-diagonal deviatoric-volumetric
          * coupling K_sdax3 characteristic of hexagonal symmetry:
          * - the pressure row: ∂σ_s/∂D_j = K_sdax3 · ∂ε'_{iTvecHex}/∂D_j for all
          *   nsvec columns of A (including the spherical column, which carries the
          *   plastic feedback of pressure-induced deviatoric stress);
          * - the direct deviatoric-from-volumetric term
          *   ∂σ'_{iTvecHex}/∂D_s = K_sdax3 · ∂ε_s/∂D_s, placed in crystal-frame rows
          *   BEFORE the caller's crystal-to-sample row rotation so that it transforms
          *   like every other deviatoric row. (An earlier version filled the whole
          *   spherical column with the pressure-row values pre-rotation, which
          *   after the row rotation produced the coupling along the doubly-rotated
          *   axis Q²e₁ instead of Qe₁ — wrong at any general orientation.)
          *
          * @tparam N Row dimension of A (ntvec)
          * @tparam M Column dimension of A (nsvec)
          * @param[out] M6 Tangent stiffness (nsvec×nsvec), crystal-frame rows
          * @param[in] A Strain sensitivities ∂(elast_d5)_i/∂(def_rate_vecds)_j (N×M)
          * @param[in] deps_sph_ddef_rate_sph ∂ε_s/∂(def_rate_vecds)_s over the step
          *            (Δt at the current level of EOS decoupling)
          * @param[in] inv_det_v_e Inverse elastic deformation determinant
          * @param[in] inv_a_vol Inverse volume factor
          */
        template<size_t N=ecmech::ntvec, size_t M=ecmech::nsvec>
        __ecmech_hdev__
        inline
        void multCauchyDif(double* const M6,
                            const double* const A,
                            double deps_sph_ddef_rate_sph,
                            double inv_det_v_e,
                            double inv_a_vol
                            ) const {
            // Diagonal contributions, all nsvec columns
            for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                double vFact = inv_det_v_e * inv_a_vol * m_K_diag[iTvec];
                for (int jSvec = 0; jSvec < ecmech::nsvec; ++jSvec) {
                    M6[ECMECH_NN_INDX(iTvec, jSvec, ecmech::nsvec)] = vFact * A[ECMECH_NM_INDX(iTvec, jSvec, N, M)];
                }
            }

            // Pressure row: dsigC_de[iSvecS, iTvecHex] * A[iTvecHex, :]
            // (the EOS part of the pressure stiffness is added by the caller)
            {
                double vFact = inv_det_v_e * inv_a_vol * m_K_sdax3;
                for (int jSvec = 0; jSvec < ecmech::nsvec; ++jSvec) {
                    M6[ECMECH_NN_INDX(iSvecS, jSvec, ecmech::nsvec)] = vFact * A[ECMECH_NM_INDX(iTvecHex, jSvec, N, M)];
                }
            }

            // Direct deviatoric response to volumetric loading, crystal frame:
            // dsigC_de[iTvecHex, iSvecS] * deps_sph/dD_s. No inv_a_vol: eval()
            // applies none to the spherical strain slot.
            M6[ECMECH_NN_INDX(iTvecHex, iSvecS, ecmech::nsvec)] +=
                inv_det_v_e * m_K_sdax3 * deps_sph_ddef_rate_sph;
        }

        /**
         * @brief Query bulk modulus.
         * @return Bulk modulus for hexagonal material
         */
        __ecmech_hdev__
        inline
        double getBulkMod( ) const {
            if (m_bulk_modulus <= 0.0) {
                ECMECH_FAIL(__func__, "bulk modulus negative -- not initialized?");
            }
            return m_bulk_modulus;
        }

        /**
         * @brief Query effective shear modulus with initialization check.
         * 
         * @param tkelv Temperature (currently unused)
         * @param pressure_EOS Pressure (currently unused)
         * @param energy_vol_ref Volumetric energy (currently unused)
         * @return Shear modulus
         * @throws Error if shear modulus is negative (indicates uninitialized model)
         */
        __ecmech_hdev__
        inline
        double getGmod(double, // tkelv
                        double, // pressure_EOS
                        double // energy_vol_ref
                        ) const {
            if (m_shear_modulus <= 0.0) {
                ECMECH_FAIL(__func__, "effective shear modulus negative -- not initialized?");
            }
            return m_shear_modulus;
        }

        private:
        /** @brief Elastic constants for hexagonal symmetry */
        double m_c11, m_c12, m_c13, m_c33, m_c44;
        /** @brief Off-diagonal coupling term specific to hexagonal symmetry */
        double m_K_sdax3;
        /** @brief Anisotropic Gruneisen parameter (thermal expansion anisotropy) */
        double m_g_vecd2;
        /** @brief Diagonal stiffness entries in deviatoric representation */
        double m_K_diag[ecmech::ntvec];
        /** @brief Bulk modulus */
        double m_bulk_modulus;
        /** @brief Effective shear modulus */
        double m_shear_modulus;
        /** @brief Index for hexagonal mode in 5-vector representation */
        static constexpr int iTvecHex = 1;
    };

}
}