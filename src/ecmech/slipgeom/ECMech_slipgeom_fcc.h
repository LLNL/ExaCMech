/**
 * @file ECMech_slipgeom_fcc.h
 * @brief Slip geometry for face-centered cubic (FCC) crystals.
 * 
 * This file defines the slip system geometry for FCC crystal structures,
 * which exhibit slip on {111} planes in <110> directions. FCC metals include
 * aluminum, copper, nickel, gold, silver, and austenitic stainless steels.
 * 
 * **FCC slip systems**:
 * - **Total systems**: 12 (always active)
 * - **Slip planes**: {111} family (4 unique planes)
 * - **Slip directions**: <110> family (3 directions per plane)
 * - **Notation**: Often called "octahedral slip"
 * 
 * **Crystallographic details**:
 * Each {111} plane has three <110> slip directions:
 * - Plane (111): Directions [0 1̄ 1], [1 0̄ 1], [1̄ 1 0]
 * - Plane (1̄11): Directions [0 1 1], [1 0̄ 1̄], [1 1̄ 0]
 * - Plane (11̄1): Directions [0 1 1̄], [1 0 1], [1 1 0̄]
 * - Plane (111̄): Directions [0 1̄ 1̄], [1 0 1̄], [1̄ 1̄ 0]
 * 
 * Total: 4 planes × 3 directions = 12 slip systems
 * 
 * **Physical behavior**:
 * - FCC metals are typically ductile due to high symmetry
 * - All 12 systems have equal Schmid factors for uniaxial loading
 * - Low temperature-dependence of critical resolved shear stress
 * - Slip occurs predominantly on {111} planes even at high temperatures
 * 
 * **Implementation details**:
 * - SlipGeomFCC : public SlipGeom<12>
 * - dynamic = false (fixed slip systems, Schmid law applies)
 * - nParams = 0 (no adjustable parameters)
 * - Slip normals (m): In {111} directions (normalized)
 * - Slip directions (s): In <110> directions (normalized)
 * 
 * **Normalization**:
 * - m vectors: 1/√3 × [±1, ±1, ±1]
 * - s vectors: 1/√2 × [0, ±1, ±1] and permutations
 * - Ensures m·s = 0 (orthogonality verified in fillFromMS)
 * 
 * **Typical kinetics models used with FCC**:
 * - Voce hardening: Phenomenological isotropic hardening
 * - KMBalD: Dislocation-density-based hardening
 * - Rate-dependent power law plasticity
 * 
 * **Crystal orientation**:
 * - Slip systems defined in crystal lattice frame
 * - Rotation tensors used to map to sample frame
 * - Lattice rotations updated during deformation
 * 
 * @see SlipGeom for base class interface
 * @see matModel for integration into crystal plasticity framework
 * @see KineticsVocePL for common FCC kinetics model
 */

#pragma once

#include "ECMech_slipgeom_base.h"

namespace ecmech {


   /**
    * @brief Face-centered cubic (FCC) slip system geometry.
    * 
    * SlipGeomFCC implements the crystallographic slip systems for FCC crystal structures.
    * FCC metals slip primarily on {111} planes in <110> directions, giving 12 slip systems
    * from the combination of 4 {111} slip planes and 3 <110> directions per plane.
    * 
    * Slip system family:
    * - 12 octahedral slip systems {111}<110>
    * 
    * Crystallographic notation:
    * - Slip planes: {111} family (close-packed planes)
    *   * (111), (1̄11), (11̄1), (111̄)
    * - Slip directions: <110> family (close-packed directions)
    *   * [011̄], [101̄], [11̄0] and permutations with sign changes
    * 
    * Slip system enumeration:
    * Systems are organized by slip plane:
    * - Systems  0-2:  (111) plane
    * - Systems  3-5:  (11̄1) plane
    * - Systems  6-8:  (1̄11) plane
    * - Systems 9-11:  (1̄1̄1) plane
    * 
    * Each plane has 3 <110> slip directions, giving 4 × 3 = 12 total systems.
    * 
    * Coordinate convention:
    * Vectors defined in crystal frame where:
    * - x₁, x₂, x₃ aligned with cubic cell edges [100], [010], [001]
    * - All slip normals and directions expressed in this basis
    * 
    * Physical properties:
    * - FCC structure: Cu, Al, Ni, Au, Ag, Pb, etc.
    * - Close-packed {111} planes have lowest energy
    * - <110> directions have shortest Burgers vector magnitude
    * - All 12 systems are crystallographically equivalent by symmetry
    * 
    * Parameters:
    * No adjustable parameters - slip geometry is purely crystallographic.
    * 
    * @ingroup ECMech_slip_geometry
    * 
    * @see SlipGeom for base class interface
    * @see SlipGeomBCC for body-centered cubic geometry
    * @see SlipGeomHCP for hexagonal close-packed geometry
    */
   class SlipGeomFCC : public SlipGeom<12>
   {
      public:
         /** @brief Slip geometry does not depend on state */
         static const bool dynamic = false;
         /** @brief Number of parameters required */
         static constexpr int nParams = 0;

         /** @brief Default constructor */
         SlipGeomFCC() = default;
         /** @brief Destructor */
         __ecmech_hdev__
         ~SlipGeomFCC() {}

         /**
          * @brief Constructor with parameters.
          * @param params Parameter array (unused for FCC)
          */
         __ecmech_hdev__
         SlipGeomFCC(const double* const params) {
            setParams(params);
         }

         /**
          * @brief Initialize slip geometry from parameter vector.
          * @param params Parameter vector (empty for FCC)
          */
         __ecmech_host__
         void setParams(const std::vector<double> & params
                        )
         {
            setParams(params.data());
         }

         /**
          * @brief Initialize slip geometry from parameter array.
          * 
          * Constructs the 12 octahedral {111}<110> slip systems for FCC crystals.
          * 
          * Implementation:
          * 1. Defines slip plane normals m for 4 {111} planes
          * 2. Defines slip directions s for 3 <110> directions per plane
          * 3. Computes Schmid tensor components P and Q via fillFromMS()
          * 4. Stores normals and directions for reference
          * 
          * @param params Parameter array (unused, FCC geometry is fixed)
          */
         __ecmech_hdev__
         void setParams(const double* const)
         {
            // m = (/ sqr3i, sqr3i, sqr3i /)
            // s = (/ zero, sqr2i, -sqr2i /)
            //
            // do not yet bother with making slip systems from symmetry group -- just write them out
            const double P3 = sqr3i, M3 = -sqr3i;
            const double P2 = sqr2i, M2 = -sqr2i;
            const double Z = zero;
            //#Slip plane normal CUB111
            const double mVecs[ nslip * ecmech::ndim ] = {
               P3, P3, P3,
               P3, P3, P3,
               P3, P3, P3,
               P3, P3, M3,
               P3, P3, M3,
               P3, P3, M3,
               P3, M3, P3,
               P3, M3, P3,
               P3, M3, P3,
               P3, M3, M3,
               P3, M3, M3,
               P3, M3, M3};
            //#Slip direction CUB110
            const double sVecs[ nslip * ecmech::ndim ] = {
               Z,  P2, M2,
               P2, Z,  M2,
               P2, M2, Z,
               Z,  P2, P2,
               P2, Z,  P2,
               P2, M2, Z,
               Z,  P2, P2,
               P2, Z,  M2,
               P2, P2, Z,
               Z,  P2, M2,
               P2, Z,  P2,
               P2, P2, Z};

            fillFromMS(this->m_P_ref_vec, this->m_Q_ref_vec,
                        mVecs, sVecs, this->nslip);

            for (int i = 0; i < nslip * ecmech::ndim; i++) {
               m_s_ref_vec[i] = sVecs[i];
               m_m_ref_vec[i] = mVecs[i];
            }
         }

         /**
          * @brief Retrieve slip geometry parameters.
          * @param params Parameter vector (empty for FCC)
          */
         __ecmech_host__
         void getParams(std::vector<double> & /* params */
                        ) const {
            // do not clear params in case adding to an existing set
         }
         
   }; // SlipGeomFCC

}