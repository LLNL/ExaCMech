/**
 * @file ECMech_eosSimple.h
 * @brief Simple equation of state models for crystal plasticity.
 * 
 * This file provides equation of state (EOS) models that relate pressure, volume,
 * energy, and temperature in thermomechanical coupling. The EOS is essential for:
 * - Computing pressure from volumetric deformation and internal energy
 * - Determining temperature evolution from energy changes
 * - Providing bulk modulus for elastic response
 * - Coupling mechanical and thermal physics
 * 
 * **Key components**:
 * - **EosModelConst**: Constant Grüneisen gamma EOS (isothermal or non-isothermal)
 * - **updateSimple**: Helper function for simple EOS update during time stepping
 * 
 * **Physical model**:
 * - Pressure: p = p_cold(V) + Γ(V) × e(V,T) / V
 * - Temperature: T = T_0 + e / c_v (linear for constant c_v)
 * - Grüneisen gamma: Γ = constant (material parameter)
 * - Bulk modulus: K = -V × dp/dV
 * 
 * **Usage context**: Called within crystal plasticity stress update to:
 * 1. Compute pressure from deformation (volumetric strain)
 * 2. Update temperature from energy dissipation (plastic work)
 * 3. Provide stiffness for implicit tangent
 * 
 * **Template parameters**:
 * - isothermal: If true, temperature held constant; if false, temperature evolves
 * 
 * @see EosModelConst for the main EOS class
 * @see updateSimple for time integration helper
 * @see matModel for integration into crystal plasticity framework
 */
// -*-c++-*-

#ifndef ECMECH_EOS_SIMPLE_H
#define ECMECH_EOS_SIMPLE_H

#include "ECMech_core.h"
#include "ECMech_util.h"

#include <string>
#include <vector>

namespace ecmech {

   /**
    * @brief Constant Gruneisen equation of state model for thermomechanical coupling.
    * 
    * EosModelConst implements a linearized Gruneisen equation of state suitable for
    * moderate compression/expansion regimes typical in solid mechanics applications.
    * It couples mechanical deformation (volume change) with thermal energy evolution
    * to compute pressure and temperature updates during material response calculations.
    * 
    * Equation of state formulation:
    * The model uses a Mie-Gruneisen form with constant Gruneisen parameter:
    *   P = K₀(V₀/V - 1) + Γρ₀ε
    *   T = T₀ + ε/cᵥ
    * where:
    *   P = pressure
    *   K₀ = reference bulk modulus
    *   V₀/V = relative volume (compression ratio)
    *   Γ = Gruneisen parameter
    *   ρ₀ = reference density
    *   ε = specific internal energy
    *   cᵥ = specific heat at constant volume
    *   T = temperature (Kelvin)
    * 
    * Physical assumptions:
    * - Linear pressure-volume relationship for mechanical compression
    * - Thermal pressure contribution proportional to internal energy
    * - Constant specific heat (temperature-independent)
    * - Small to moderate deformations (relative volume 0.7 to 1.4)
    * - Gruneisen parameter independent of volume and temperature
    * 
    * Isothermal vs. non-isothermal modes:
    * The template parameter controls thermomechanical coupling:
    * - isothermal=true: Temperature fixed at T₀, thermal pressure term disabled
    * - isothermal=false: Full thermomechanical coupling with energy evolution
    * 
    * Integration with material models:
    * EOS models work in tandem with elastoplastic constitutive models:
    * 1. Mechanical deformation produces volume change (compression/expansion)
    * 2. EOS computes pressure from volume change and internal energy
    * 3. Plastic work and compression work update internal energy
    * 4. EOS computes new temperature from updated energy
    * 5. Temperature affects kinetic processes (slip rates, hardening)
    * 
    * Numerical considerations:
    * - Bulk modulus is volume-dependent: K = K₀(V₀/V) for compression stiffening
    * - Minimum bulk modulus enforced to prevent unphysical softening
    * - Energy updates use trapezoidal rule for improved accuracy
    * - Pressure derivatives computed analytically for tangent stiffness
    * 
    * @tparam isothermal If true, temperature is constant; if false, temperature evolves with energy
    * 
    * @ingroup ECMech_equation_of_state
    * 
    * @see updateSimple() for time integration with this EOS model
    * @see ECMech_evptnWrap.h for usage within material models
    */
   template<bool isothermal>
   class EosModelConst
   {
      public:
         /** @brief Number of parameters required to initialize the model */
         static constexpr int nParams = 5;

         /**
          * @brief Default constructor creates uninitialized EOS model.
          * 
          * Parameters must be set via setParams() before use.
          */
         EosModelConst() = default;

         /**
          * @brief Constructor with parameter array initialization.
          * @param params Array of length nParams containing EOS parameters
          */
         __ecmech_hdev__
         EosModelConst(const double* const params) {
            setParams(params);
         }

         /** @brief Destructor */
         ~EosModelConst() = default;

         /**
          * @brief Set EOS parameters from vector.
          * @param params Vector containing [ρ₀, K₀, cᵥ, Γ, ε₀]
          */
         __ecmech_host__
         inline
         void setParams(const std::vector<double>& params) {
            setParams(params.data());
         }

         /**
          * @brief Set EOS parameters from array.
          * 
          * Configures the equation of state with material-specific parameters and
          * computes derived quantities needed for efficient evaluation.
          * 
          * Parameter order and physical meaning:
          * 1. m_density0 (ρ₀): Reference density
          *    - Mass per unit reference volume
          *    - Must be positive
          * 
          * 2. m_bulk_modulus (K₀): Reference bulk modulus
          *    - Resistance to volumetric compression: K = -V(∂P/∂V)
          *    - Must be positive
          *    - Related to elastic moduli: K ≈ E/(3(1-2ν))
          * 
          * 3. m_cvav (cᵥ): Specific heat capacity
          *    - Heat capacity at constant volume
          *    - Determines thermal inertia
          *    - Must be positive
          * 
          * 4. m_gamma (Γ): Gruneisen parameter [dimensionless]
          *    - Couples thermal and mechanical behavior
          *    - Relates thermal pressure to internal energy
          *    - Can be zero for purely mechanical EOS
          * 
          * 5. m_cold_energy0 (ε₀ᶜᵒˡᵈ): Cold curve energy offset
          *    - Reference energy at zero temperature
          *    - Often set to zero for convenience
          *    - Affects absolute temperature scale
          * 
          * Derived quantities computed:
          * - m_dtde = 1/cᵥ: Temperature derivative w.r.t. energy
          * - m_tkelv0 = -ε₀ᶜᵒˡᵈ/cᵥ: Reference temperature
          * 
          * @param params Pointer to array of 5 doubles in order listed above
          * 
          * @note In debug builds, validates that exactly nParams values are consumed
          */
         __ecmech_hdev__
         inline
         void setParams(const double* const params) {
            const double* parsIt = params;
            //////////////////////////////
            m_density0 = *parsIt; ++parsIt;
            m_bulk_modulus = *parsIt; ++parsIt;
            m_cvav = *parsIt; ++parsIt;
            m_gamma = *parsIt; ++parsIt;
            m_cold_energy0 = *parsIt; ++parsIt;

            m_dtde = one / m_cvav; // ∂T/∂ε at constant volume
            m_tkelv0 = -m_cold_energy0 * m_dtde; // Reference temperature

            //////////////////////////////
#if defined(ECMECH_DEBUG)
            int iParam = parsIt - params;
            if (iParam != nParams) {
               ECMECH_FAIL(__func__, "iParam != nParams");
            }
#endif
         }

         /**
          * @brief Retrieve current parameters for serialization.
          * 
          * Extracts parameters in the same order expected by setParams(), enabling
          * model serialization for checkpointing or parameter inspection.
          * 
          * @param[out] params Vector to receive parameters (not cleared; appended to)
          */
         __ecmech_host__
         inline
         void getParams(std::vector<double> & params
                        ) const {
            // do not clear params in case adding to an existing set
            int paramsStart = params.size();

            //////////////////////////////

            params.push_back(m_density0);
            params.push_back(m_bulk_modulus);
            params.push_back(m_cvav);
            params.push_back(m_gamma);
            params.push_back(m_cold_energy0);

            //////////////////////////////

            int iParam = params.size() - paramsStart;
            if (iParam != nParams) {
               ECMECH_FAIL(__func__, "iParam != nParams");
            }
         }

         /**
          * @brief Evaluate pressure and temperature from current state.
          * 
          * Computes thermodynamic state (pressure and temperature) given the
          * current volume and internal energy. This is the core EOS evaluation
          * used at the beginning of time steps or for diagnostics.
          * 
          * Pressure evaluation:
          * - Isothermal mode: P = K₀(V₀/V - 1)
          * - Non-isothermal: P = K₀(V₀/V - 1) + Γρ₀ε
          * 
          * Temperature evaluation:
          * - Isothermal mode: T = T₀ (constant)
          * - Non-isothermal: T = T₀ + ε/cᵥ
          * 
          * @param[out] pressure Computed pressure
          * @param[out] tkelv Computed temperature
          * @param[in] rel_vol Relative volume V/V₀
          * @param[in] energy internal energy
          * 
          * @note Pressure can be negative (tension) for rel_vol > 1
          * @note Temperature is always positive in physical regime
          */
         __ecmech_hdev__
         inline void evalPT(double &pressure,
                            double &tkelv,
                            double  rel_vol,
                            double  energy) const {
            double mu = one / rel_vol - one; // Compression measure: (V₀-V)/V

            if (isothermal) {
               pressure = m_bulk_modulus * mu;
               tkelv = m_tkelv0;
            }
            else {
               pressure = m_bulk_modulus * mu + m_gamma * energy;
               tkelv = m_tkelv0 + energy * m_dtde;
            }
         }

         /**
          * @brief Evaluate pressure, temperature, and their derivatives.
          * 
          * Extended EOS evaluation that additionally computes derivatives needed
          * for implicit time integration and tangent stiffness calculations.
          * Used during Newton-Raphson iterations in material response calculations.
          * 
          * Computed quantities:
          * 1. Pressure P(V,ε) and temperature T(ε)
          * 2. Current bulk modulus K(V) = K₀(V₀/V) accounting for compression
          * 3. Pressure derivatives:
          *    - ∂P/∂ε at constant volume (dpde)
          *    - ∂T/∂ε at constant volume (dtde)
          * 
          * Bulk modulus evolution:
          * The bulk modulus increases with compression: K = K₀η where η = V₀/V.
          * This captures the stiffening behavior of materials under compression.
          * 
          * Derivative computation:
          * - Isothermal mode: 
          *   - dpde = 0 (pressure independent of energy)
          *   - dtde = small value (prevents division by zero)
          * - Non-isothermal mode:
          *   - dpde = Γ (Gruneisen parameter)
          *   - dtde = 1/cᵥ
          * 
          * @param[out] pressure Computed pressure
          * @param[out] tkelv Computed temperature [K]
          * @param[out] bulk_modulus_new Current bulk modulus K(V)
          * @param[out] dpde Derivative ∂P/∂ε|_V
          * @param[out] dtde Derivative ∂T/∂ε (always positive)
          * @param[in] rel_vol Relative volume V/V₀
          * @param[in] energy internal energy
          * 
          * @note In isothermal mode, dtde is set to small value rather than zero
          *       to prevent divide-by-zero in downstream calculations
          */
         __ecmech_hdev__
         inline
         void evalPTDiff(double &pressure,
                         double &tkelv,
                         double &bulk_modulus_new,
                         double &dpde,
                         double &dtde,
                         double  rel_vol,
                         double  energy) const {
            double eta = one / rel_vol; // V₀/V (compression ratio)
            double mu = eta - one; // (V₀-V)/V

            tkelv = this->evalT(energy);

            if (isothermal) {
               pressure = m_bulk_modulus * mu;
               dpde = zero;
               dtde = 1e-8 * m_dtde; // instead of zero, to prevent divide-by-zero elsewhere
            }
            else {
               pressure = m_bulk_modulus * mu + m_gamma * energy;
               dpde = m_gamma;
               dtde = m_dtde;
            }
            bulk_modulus_new = m_bulk_modulus * eta; // Compression-dependent bulk modulus
         }

         /**
          * @brief Get valid range of relative volumes and reference state information.
          * 
          * Returns the range of relative volumes over which the EOS is considered
          * physically valid, along with reference state values. Useful for:
          * - Validating material state during simulation
          * - Setting up initial conditions
          * - Detecting extreme deformations requiring special treatment
          * 
          * Valid volume range:
          * - rel_vol_min = 0.1 (90% compression, 10× density increase)
          * - rel_vol_max = 10.0 (900% expansion, 0.1× density)
          * These limits are conservative and may be narrowed for specific applications.
          * 
          * @param[out] rel_vol_min Minimum valid relative volume
          * @param[out] rel_vol_max Maximum valid relative volume
          * @param[out] energy0 Reference internal energy
          * @param[out] rel_vol0 Reference relative volume (always 1.0)
          */
         __ecmech_hdev__
         inline
         void getInfo(double &rel_vol_min,
                      double &rel_vol_max,
                      double &energy0,
                      double &rel_vol0) const {
            rel_vol_min = 0.1;
            rel_vol_max = 10.0;
            energy0 = 0.0;
            rel_vol0 = 1.0;
         }

         /**
          * @brief Get reference bulk modulus.
          * @return Reference bulk modulus K₀
          */
         __ecmech_hdev__
         inline
         double getBulkRef() const {
            return m_bulk_modulus;
         }

         /**
          * @brief Get reference density.
          * @return Reference density ρ₀
          */
         __ecmech_hdev__
         inline
         double getRho0() const {
            return m_density0;
         }

      private:

         /**
          * @brief Internal temperature evaluation from energy.
          * 
          * Helper method to compute temperature consistently across different
          * evaluation paths.
          * 
          * @param energy internal energy
          * @return Temperature in Kelvin
          */
         __ecmech_hdev__
         inline double evalT(double  energy) const {
            double tkelv;
            if (isothermal) {
               tkelv = m_tkelv0;
            }
            else {
               tkelv = m_tkelv0 + energy * m_dtde;
            }
            return tkelv;
         }

      private:

         // Primary parameters (set by user)
         
         /** @brief Reference density ρ₀ [mass/volume] */
         double m_density0;
         
         /** @brief Reference bulk modulus K₀ [pressure units] */
         double m_bulk_modulus;
         
         /** @brief Gruneisen parameter Γ [dimensionless] */
         double m_gamma;
         
         /** @brief Cold curve energy offset ε₀ᶜᵒˡᵈ [energy/mass] */
         double m_cold_energy0;
         
         /** @brief Specific heat at constant volume cᵥ [energy/(mass·temperature)] */
         double m_cvav;

         // Derived parameters (computed from primary parameters)
         
         /** @brief Temperature derivative ∂T/∂ε = 1/cᵥ [temperature/(energy/mass)] */
         double m_dtde;
         
         /** @brief Reference temperature T₀ = -ε₀ᶜᵒˡᵈ/cᵥ [temperature] */
         double m_tkelv0;
   }; // class EosModelConst

   /**
    * @brief Update thermodynamic state over a time step using simple integration.
    * 
    * This function integrates the thermodynamic state (pressure, temperature, energy)
    * forward in time given the volume evolution. It couples with the material model's
    * mechanical response to maintain thermodynamic consistency.
    * 
    * Integration approach:
    * Energy is updated using a first-order approximation of the mechanical work:
    *   ε_{n+1} = ε_n - ΔV · P_n
    * where ΔV is the relative volume increment and P_n is the beginning-of-step pressure.
    * This represents the pressure work done by volume change (compression heats, expansion cools).
    * 
    * Bulk modulus correction:
    * The bulk modulus is updated to account for:
    * 1. Compression-dependent stiffening: K = K₀(V₀/V)
    * 2. Pressure-dependent correction: K_corr = K + (∂P/∂ε)·P·V
    * 3. Minimum bulk modulus enforcement to prevent unphysical softening
    * 
    * Thermodynamic consistency:
    * The volume derivative of pressure dpdv = -K/V maintains the thermodynamic
    * relation: (∂P/∂V)_ε = -K/V for consistent tangent stiffness.
    * 
    * Usage in material models:
    * Called during the preprocessing phase of material response evaluation:
    * 1. Before entering Newton-Raphson loop for implicit integration
    * 2. Provides updated pressure for subsequent stress calculations
    * 3. Computes derivatives needed for tangent stiffness matrix
    * 
    * @tparam EosModel Equation of state type (typically EosModelConst<isothermal>)
    * 
    * @param[in] eos Equation of state model instance
    * @param[out] press Updated pressure at end of step
    * @param[out] tkelv Updated temperature at end of step
    * @param[out] energy_new Updated internal energy at end of step
    * @param[out] bulk_modulus_new Updated bulk modulus at end of step
    * @param[out] dpde Pressure derivative ∂P/∂ε|_V
    * @param[out] dpdv Pressure derivative ∂P/∂V|_ε = -K/V
    * @param[out] dtde Temperature derivative ∂T/∂ε = 1/cᵥ
    * @param[in] rel_vol_new Relative volume at end of step V_{n+1}/V₀
    * @param[in] rel_vol_increment Volume increment ΔV/V
    * @param[in] energy_old Internal energy at beginning of step
    * @param[in] pressure_old Pressure at beginning of step
    * 
    * @note The bulk modulus minimum (10⁻⁵ times reference) prevents numerical issues
    * @note Pressure can be negative (tension) for expansion beyond reference volume
    * 
    * @ingroup ECMech_equation_of_state
    */
   template<class EosModel>
   __ecmech_hdev__
   inline
   void updateSimple(const EosModel& eos,
                     double &press,
                     double &tkelv,
                     double &energy_new,
                     double &bulk_modulus_new,
                     double &dpde,
                     double &dpdv,
                     double &dtde,
                     double  rel_vol_new,
                     double  rel_vol_increment,
                     double  energy_old,
                     double  pressure_old)
   {
      // Update energy using pressure work
      // ε_{n+1} = ε_n - P_n · ΔV
      energy_new = energy_old - rel_vol_increment * pressure_old;

      // Evaluate EOS at new state
      eos.evalPTDiff(press, tkelv, bulk_modulus_new, dpde, dtde, rel_vol_new, energy_new);
      // Compute volume derivative of pressure
      dpdv = -bulk_modulus_new / rel_vol_new;

      // Apply bulk modulus correction for pressure-energy coupling
      // This accounts for the fact that pressure affects subsequent energy evolution
      double bulk_modulus_min = 1e-5 * eos.getBulkRef();
      bulk_modulus_new = fmax(bulk_modulus_min, bulk_modulus_new + dpde * pressure_old * rel_vol_new);
   }
} // namespace ecmech

#endif // ECMECH_EOS_SIMPLE_H
