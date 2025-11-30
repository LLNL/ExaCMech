/**
 * @file ECMech_slipgeom_hcp.h
 * @brief Slip geometry for hexagonal close-packed (HCP) crystals.
 * 
 * This file defines slip system geometries for HCP crystal structures, which
 * exhibit complex slip behavior due to their lower symmetry compared to cubic
 * crystals. HCP metals include titanium, magnesium, zinc, zirconium, beryllium,
 * and many rare earth elements.
 * 
 * **HCP crystal structure**:
 * - Lattice parameters: a (basal), c (height)
 * - c/a ratio: Material-dependent (ideal = √(8/3) ≈ 1.633)
 * - Lower symmetry than cubic crystals
 * - Anisotropic mechanical properties
 * 
 * **HCP slip system families**:
 * 
 * 1. **Basal <a> slip** (3 systems):
 *    - Slip plane: (0001) basal plane
 *    - Slip directions: <112̄0> family (3 directions)
 *    - Easiest slip mode in most HCP metals
 *    - Active at all temperatures
 *    - Cannot accommodate c-axis strain
 * 
 * 2. **Prismatic <a> slip** (3 systems):
 *    - Slip planes: {101̄0} prismatic planes (3 planes)
 *    - Slip directions: <112̄0> family
 *    - Important for general plasticity
 *    - CRSS typically higher than basal
 *    - Cannot accommodate c-axis strain
 * 
 * 3. **Pyramidal <a> slip** (6 systems):
 *    - Slip planes: {101̄1} first-order pyramidal planes
 *    - Slip directions: <112̄0> family (2 directions per plane)
 *    - Higher CRSS than basal and prismatic
 *    - Active at elevated temperatures or high stress
 *    - Cannot accommodate c-axis strain
 * 
 * 4. **Pyramidal <c+a> slip** (12 systems):
 *    - Slip planes: {101̄1} first-order pyramidal planes
 *    - Slip directions: <112̄3> family (4 directions per plane)
 *    - CRITICAL for c-axis strain accommodation
 *    - Highest CRSS of all slip modes
 *    - Essential for ductility and avoiding twinning
 *    - Often requires elevated temperature
 * 
 * Total: 3 + 3 + 6 + 12 = 24 slip systems
 * 
 * **Implementation: SlipGeomHCPaBRYcaY1**:
 * - SlipGeomHCPaBRYcaY1 : public SlipGeom<24>
 * - Name origin: "aBRYcaY1" = <a> Basal + pRismatic + pYramidal + <c+a> pyramidal type-Y1
 * - Corresponds to EVP_HCP_a_BRY_ca_Y1 (code 32) in Fortran version
 * - dynamic = false (fixed slip systems)
 * - nParams = 1 (c/a ratio)
 * 
 * **c/a ratio dependence**:
 * The c/a ratio affects:
 * - Pyramidal plane orientations
 * - Slip direction components
 * - Relative CRSS values
 * - Deformation mode selection
 * 
 * **Miller-Bravais indices**:
 * HCP uses 4-index notation [hkil] where i = -(h+k):
 * - Planes: (hkil) with h+k+i = 0
 * - Directions: [uvtw] with u+v+t = 0
 * - Conversion to Cartesian via miller_to_orthog_sngl()
 * 
 * **Coordinate system**:
 * - x₁-axis: Along a₁ basal vector
 * - x₂-axis: 30° from a₂ (between a₁ and a₂)
 * - x₃-axis: Along c-axis (perpendicular to basal plane)
 * 
 * **Physical behavior**:
 * - Strong plastic anisotropy
 * - Limited slip systems → twinning often active
 * - <c+a> slip critical for ductility
 * - Temperature-dependent slip mode selection
 * - Texture evolution strongly affects properties
 * 
 * **Typical kinetics models for HCP**:
 * - KMBalD with separate hardening for each slip family
 * - Temperature-dependent CRSS ratios
 * - Latent hardening between slip families
 * - Coupled slip-twinning models
 * 
 * **Hardening interactions**:
 * Different slip families exhibit different hardening interactions:
 * - Basal-basal: Strong latent hardening
 * - Basal-prismatic: Moderate interaction
 * - <a>-<c+a>: Weak interaction
 * - <c+a>-<c+a>: Important at high strains
 * 
 * **Limitations**:
 * - Does not include twinning systems
 * - Single pyramidal variant (other variants possible)
 * - Second-order pyramidal systems not included
 * - Non-Schmid effects not implemented
 * 
 * **Future extensions?**:
 * - SlipGeomHCPaBRYcaY2: Second-order <c+a> pyramidal {112̄2}<112̄3>
 * - Combined slip-twinning geometries
 * - Non-Schmid effects for HCP
 * 
 * @see SlipGeom for base class interface
 * @see miller_to_orthog_sngl for Miller index conversion
 * @see KineticsKMBalD for typical HCP kinetics
 */

#pragma once

#include "ECMech_slipgeom_base.h"

namespace ecmech {

   /**
    * @brief Hexagonal close-packed (HCP) slip system geometry with multiple slip families.
    * 
    * SlipGeomHCPaBRYcaY1 implements crystallographic slip systems for HCP crystal structures.
    * HCP metals exhibit highly anisotropic plastic behavior due to their non-cubic symmetry
    * and strong dependence on c/a ratio. Multiple slip families are required to accommodate
    * arbitrary plastic deformation.
    * 
    * Slip system families (24 total systems):
    * 
    * 1. **Basal <a> slip** (3 systems, indices 0-2):
    *    - Slip planes: {0001} basal plane (perpendicular to c-axis)
    *    - Slip directions: <11̄20> in-plane <a> directions
    *    - Easiest slip mode, typically has lowest CRSS
    *    - Cannot accommodate strain along c-axis
    * 
    * 2. **Prismatic <a> slip** (3 systems, indices 3-5):
    *    - Slip planes: {11̄00} prismatic planes (parallel to c-axis)
    *    - Slip directions: <11̄20> <a> directions
    *    - Moderately difficult, required for in-plane deformation
    *    - Cannot accommodate strain along c-axis
    * 
    * 3. **Pyramidal <a> slip** (6 systems, indices 6-11):
    *    - Slip planes: {10̄11} first-order pyramidal planes
    *    - Slip directions: <11̄20> <a> directions  
    *    - Higher CRSS than basal or prismatic
    *    - Still cannot accommodate c-axis strain
    * 
    * 4. **Pyramidal <c+a> slip** (12 systems, indices 12-23):
    *    - Slip planes: {10̄11} first-order pyramidal planes
    *    - Slip directions: <11̄23> <c+a> directions with c-component
    *    - Essential for c-axis strain accommodation
    *    - Typically highest CRSS, activated at high stress or temperature
    * 
    * Crystal structure:
    * HCP unit cell characterized by c/a ratio:
    * - Ideal c/a = √(8/3) ≈ 1.633
    * - c/a ratio affects slip plane orientations and slip resistance
    * 
    * Coordinate convention:
    * - x₁, x₂: Basal plane (perpendicular to c-axis)
    * - x₃: c-axis direction
    * - Slip geometry depends on c/a ratio parameter
    * 
    * Miller-Bravais indices:
    * HCP uses 4-index notation [uvtw] where t = -(u+v):
    * - Directions: <11̄20> for <a>, <11̄23> for <c+a>
    * - Planes: {0001} basal, {11̄00} prismatic, {10̄11} pyramidal
    * 
    * Parameter:
    * c/a ratio must be provided to correctly compute slip geometry for
    * pyramidal systems where plane orientations depend on lattice parameters.
    * 
    * @ingroup ECMech_slip_geometry
    * 
    * @see SlipGeom for base class interface
    * @see SlipGeomFCC for face-centered cubic geometry
    * @see SlipGeomBCC for body-centered cubic geometry
    */
   class SlipGeomHCPaBRYcaY1 : public SlipGeom<3 + 3 + 6 + 12>
   {
      public:
         /** @brief Slip geometry does not depend on state */
         static const bool dynamic = false;

         /** @brief Total number of slip systems */
         // 3  basal <a> + 3  prismatic <a> + 6  pyramidal <a> + 12 pyramidal <c+a>
         //static constexpr int nslip = 3 + 3 + 6 + 12;
         
         /** @brief Number of parameters required (c/a ratio) */
         static constexpr int nParams = 1;

         /** @brief Default constructor */
         SlipGeomHCPaBRYcaY1() = default;
         /** @brief Destructor */
         __ecmech_hdev__
         ~SlipGeomHCPaBRYcaY1() {}

         /**
          * @brief Constructor with parameters.
          * @param params Parameter array containing c/a ratio
          */
         __ecmech_hdev__
         SlipGeomHCPaBRYcaY1(const double* const params) {
            setParams(params);
         }


         /**
          * @brief Initialize from parameter vector.
          * @param params Parameter vector containing c/a ratio
          */
         __ecmech_host__
         void setParams(const std::vector<double> & params)
         {
            setParams(params.data());
         }

         /**
          * @brief Initialize HCP slip geometry from c/a ratio.
          * 
          * Constructs 24 slip systems across four families, with pyramidal plane orientations
          * computed from the lattice c/a ratio parameter.
          * 
          * Implementation sequence:
          * 1. Extract c/a ratio from parameters
          * 2. Compute pyramidal {10̄11} plane orientations (c/a dependent)
          * 3. Construct slip plane normals and directions for all 24 systems
          * 4. Compute Schmid tensors via fillFromMS()
          * 
          * Miller-Bravais to Cartesian conversion:
          * Uses miller_to_orthog_sngl() helper to convert 4-index [uvtw] notation
          * to 3D Cartesian coordinates in the crystal frame.
          * 
          * Slip system organization:
          * - Systems 0-2: Basal {0001}<11̄20>
          * - Systems 3-5: Prismatic {11̄00}<11̄20>
          * - Systems 6-11: Pyramidal {10̄11}<11̄20>
          * - Systems 12-23: Pyramidal {10̄11}<11̄23>
          * 
          * @param params Parameter array: [0] = c/a ratio
          */
         __ecmech_hdev__
         void setParams(const double* const params)
         {
            const double* parsIt = params;

            m_cOverA = *parsIt; ++parsIt;

            // pyramidal 10-11 1-210 depends on c/a
            //
            double m_ya[ecmech::ndim], s_ya[ecmech::ndim];
            {
               double an[ecmech::nMiller] = { one, zero, -one, one }; // plane
               double ab[ecmech::nMiller] = { one, -two, one, zero }; // direction
               //
               miller_to_orthog_sngl(an, ab,
                                     m_ya, s_ya,
                                     m_cOverA);
            }
            double m_ya_pp = sqrt(1.0 - m_ya[2] * m_ya[2]);

            // pyramidal 10-11 -1-123 depends on c/a
            //
            double m_y1ca[ecmech::ndim], s_y1ca[ecmech::ndim];
            {
               double an[ecmech::nMiller] = { one, zero, -one, one }; // plane
               double ab[ecmech::nMiller] = { -one, -one, two, three }; // direction
               //
               miller_to_orthog_sngl(an, ab,
                                     m_y1ca, s_y1ca,
                                     m_cOverA);
            }
            double m_y1ca_pp = sqrt(1.0 - m_y1ca[2] * m_y1ca[2]);
            double s_y1ca_pp = sqrt(1.0 - s_y1ca[2] * s_y1ca[2]);

            const double mVecs[ nslip * ecmech::ndim ] = {
               zero, zero, one,
               zero, zero, one,
               zero, zero, one,

               -halfsqr3, onehalf, zero,
               -halfsqr3, -onehalf, zero,
               zero, -one, zero,

               m_ya[0], m_ya[1], m_ya[2],
               m_ya[0], -m_ya[1], -m_ya[2],
               m_ya[0], m_ya[1], -m_ya[2],
               -m_ya[0], m_ya[1], -m_ya[2],
               zero, m_ya_pp, -m_ya[2],
               zero, -m_ya_pp, -m_ya[2],

               m_y1ca[0], m_y1ca[1], m_y1ca[2],
               m_y1ca[0], -m_y1ca[1], -m_y1ca[2],
               m_y1ca[0], m_y1ca[1], -m_y1ca[2],
               zero, m_y1ca_pp, -m_y1ca[2],
               -m_y1ca[0], m_y1ca[1], -m_y1ca[2],
               -m_y1ca[0], -m_y1ca[1], -m_y1ca[2],
               zero, -m_y1ca_pp, -m_y1ca[2],
               zero, m_y1ca_pp, m_y1ca[2],
               -m_y1ca[0], m_y1ca[1], m_y1ca[2],
               -m_y1ca[0], -m_y1ca[1], m_y1ca[2],
               m_y1ca[0], -m_y1ca[1], m_y1ca[2],
               zero, -m_y1ca_pp, m_y1ca[2]
            };
            const double sVecs[ nslip * ecmech::ndim ] = {
               onehalf, halfsqr3, zero,
               onehalf, -halfsqr3, zero,
               one, zero, zero,

               onehalf, halfsqr3, zero,
               onehalf, -halfsqr3, zero,
               one, zero, zero,

               s_ya[0], s_ya[1], zero,
               s_ya[0], -s_ya[1], zero,
               -s_ya[0], -s_ya[1], zero,
               -s_ya[0], s_ya[1], zero,
               -one, zero, zero,
               one, zero, zero,

               s_y1ca[0], s_y1ca[1], s_y1ca[2],
               s_y1ca[0], -s_y1ca[1], -s_y1ca[2],
               -s_y1ca_pp, zero, -s_y1ca[2],
               s_y1ca[0], s_y1ca[1], -s_y1ca[2],
               -s_y1ca[0], s_y1ca[1], -s_y1ca[2],
               s_y1ca_pp, zero, -s_y1ca[2],
               -s_y1ca[0], -s_y1ca[1], -s_y1ca[2],
               -s_y1ca[0], s_y1ca[1], s_y1ca[2],
               s_y1ca_pp, zero, s_y1ca[2],
               -s_y1ca[0], -s_y1ca[1], s_y1ca[2],
               -s_y1ca_pp, zero, s_y1ca[2],
               s_y1ca[0], -s_y1ca[1], s_y1ca[2]
            };

            fillFromMS(this->m_P_ref_vec, this->m_Q_ref_vec,
                       mVecs, sVecs, this->nslip);

            for (int i = 0; i < nslip * ecmech::ndim; i++) {
               m_s_ref_vec[i] = sVecs[i];
               m_m_ref_vec[i] = mVecs[i];
            }

#if defined(ECMECH_DEBUG)
            int iParam = parsIt - params;
            if (iParam != nParams) {
               ECMECH_FAIL(__func__, "iParam != nParams");
            }
#endif
         }

         /**
          * @brief Retrieve slip geometry parameters.
          * @param params Parameter vector to receive c/a ratio
          */
         __ecmech_host__
         void getParams(std::vector<double> & params
                        ) const {
#ifdef ECMECH_DEBUG
            // do not clear params in case adding to an existing set
            int paramsStart = params.size();
#endif
            params.push_back(m_cOverA);
#ifdef ECMECH_DEBUG
            assert((params.size() - paramsStart) == nParams);
#endif
         }

      private:
         /** @brief Lattice c/a ratio, controls pyramidal plane orientations */
         double m_cOverA;
    
   }; // SlipGeomHCPaBRYcaY1

}