/**
 * @file ECMech_slipgeom_base.h
 * @brief Base class and utilities for crystal slip geometry definitions.
 * 
 * This file provides the foundational infrastructure for defining slip systems
 * in crystal plasticity models. Slip systems consist of slip planes (normals)
 * and slip directions that define the crystallographic planes and directions
 * along which plastic deformation occurs via dislocation motion.
 * 
 * **Core concepts**:
 * - **Slip plane normal (m)**: Unit vector perpendicular to the slip plane
 * - **Slip direction (s)**: Unit vector along the slip direction (Burgers vector)
 * - **Schmid tensor**: Symmetric part P = ½(s⊗m + m⊗s) for resolved shear stress
 * - **Anti-symmetric part Q**: ½(s⊗m - m⊗s) for lattice spin contributions
 * 
 * **Key components**:
 * - **SlipGeom<nslip, nSlipExtra>**: Template base class for all slip geometries
 *   - nslip: Number of slip systems (compile-time constant)
 *   - nSlipExtra: Extra variables for non-Schmid effects (default 0)
 * 
 * - **fillFromMS()**: Utility function to compute P and Q tensors from m and s vectors
 *   - Converts slip plane normals and directions to Schmid factors
 *   - Handles deviatoric (5-component) and skew (3-component) representations
 * 
 * **Design patterns**:
 * - Template on number of slip systems for compile-time optimization
 * - Virtual interface for dynamic slip systems (non-Schmid effects)
 * - Static vs dynamic modes:
 *   - dynamic = false: Fixed slip systems (FCC, standard BCC, HCP)
 *   - dynamic = true: State-dependent slip systems (BCC non-Schmid)
 * 
 * **Memory layout**:
 * - m_P_ref_vec: [ntvec × nslip] Schmid factors for deviatoric stress
 * - m_Q_ref_vec: [nwvec × nslip] Anti-Schmid factors for lattice spin
 * - m_m_ref_vec: [ndim × nslip] Slip plane normals
 * - m_s_ref_vec: [ndim × nslip] Slip directions
 * 
 * **Usage in crystal plasticity**:
 * 1. Resolve Kirchhoff stress onto slip systems: τᵅ = P^α : τ
 * 2. Compute plastic deformation rate: D^p = Σ γ̇ᵅ P^α
 * 3. Compute plastic spin: W^p = Σ γ̇ᵅ Q^α
 * 4. Update lattice orientation and elastic strain
 * 
 * **Derived classes**:
 * - SlipGeomFCC: Face-centered cubic (12 systems)
 * - SlipGeomBCC: Body-centered cubic (12, 24, or 48 systems)
 * - SlipGeomBCCNonSchmid: BCC with non-Schmid effects (dynamic)
 * - SlipGeomBCCPencil: BCC with pencil glide non-Schmid (dynamic)
 * - SlipGeomHCPaBRYcaY1: Hexagonal close-packed (24 systems)
 * 
 * **Performance considerations**:
 * - Compile-time nslip enables loop unrolling and vectorization
 * - Deviatoric representation reduces from 6 to 5 components
 * - Precomputed P and Q tensors avoid repeated calculations
 * 
 * @see SlipGeomFCC for FCC implementation
 * @see SlipGeomBCC for BCC implementation
 * @see SlipGeomHCPaBRYcaY1 for HCP implementation
 * @see fillFromMS for tensor computation from Miller indices
 */

#pragma once

#include "ECMech_core.h"
#include "ECMech_util.h"

namespace ecmech {
   /**
    * @brief Compute Schmid tensor components from slip plane normals and directions.
    * 
    * Transforms crystallographic slip system descriptions (plane normal m, slip direction s)
    * into deviatoric (P) and spin (Q) tensor components used in plastic rate calculations.
    * 
    * For each slip system α, forms the Schmid tensor:
    *   S^α = s^α ⊗ m^α
    * 
    * Then decomposes into:
    * - P^α: Symmetric deviatoric part (5 components)
    * - Q^α: Skew-symmetric part (3 components as axial vector)
    * 
    * @param[out] P Symmetric deviatoric Schmid components, shape [ntvec, nslip]
    * @param[out] Q Skew-symmetric Schmid components, shape [nwvec, nslip]
    * @param[in] mVecs Slip plane normals, shape [nslip, ndim]
    * @param[in] sVecs Slip directions, shape [nslip, ndim]
    * @param[in] nslip Number of slip systems
    */
   __ecmech_hdev__
   static void
   fillFromMS(double* const P, // ntvec * nslip
              double* const Q, // nwvec * nslip
              const double* const mVecs, // nslip * ndim
              const double* const sVecs, // nslip * ndim
              const int nslip)
   {
      for (int iSlip = 0; iSlip<nslip; ++iSlip) {
         const double* mVec = &(mVecs[iSlip * ecmech::ndim]);
         const double* sVec = &(sVecs[iSlip * ecmech::ndim]);
#ifndef NO_CHECKS
         if (fabs(vecsyadotb<ecmech::ndim>(mVec, sVec)) > idp_eps_sqrt) {
            ECMECH_FAIL(__func__, "internal error");
         }
#endif

         // CALL vec_x_vect_mn(crys%vecs(:,is),crys%vecm(:,is),crys%t_ref(:,:,is),DIMS,DIMS)
         double T_ref[ ecmech::ndim * ecmech::ndim ];
         vecsMaTb<ndim>(T_ref, sVec, mVec);

         double P_vecd[ ecmech::ntvec ];
         double Q_veccp[ ecmech::nwvec ];
         matToPQ(P_vecd, Q_veccp, T_ref);

         for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
            P[ECMECH_NM_INDX(iTvec, iSlip, ecmech::ntvec, nslip)] = P_vecd[iTvec];
         }

         for (int iWvec = 0; iWvec < ecmech::nwvec; ++iWvec) {
            Q[ECMECH_NM_INDX(iWvec, iSlip, ecmech::nwvec, nslip)] = Q_veccp[iWvec];
         }

         //
         // in some approaches, it is useful to form the outer product of P_vecd with itself, for tangent stiffness contributions
      }
   }

   /**
    * @brief Base template class for slip system geometry.
    * 
    * SlipGeom provides the crystallographic description of slip systems for a material.
    * It stores and computes the geometric relationships between slip plane normals,
    * slip directions, and the corresponding Schmid tensors needed for crystal plasticity.
    * 
    * Core functionality:
    * - Storage of slip plane normals (m) and slip directions (s) for all systems
    * - Computation of symmetric (P) and skew (Q) Schmid tensor components
    * - Resolution of stress onto slip systems to compute resolved shear stresses
    * - Provision of geometric data for plastic strain rate calculations
    * 
    * Schmid tensor decomposition:
    * For slip system α with plane normal m^α and slip direction s^α:
    *   S^α = s^α ⊗ m^α (dyadic product)
    * 
    * Decomposed as:
    *   S^α = P^α + Q^α
    * where:
    * - P^α: Symmetric deviatoric part (drives plastic stretching)
    * - Q^α: Skew-symmetric part (drives plastic spin/rotation)
    * 
    * Crystal kinematics:
    * Plastic deformation rate: D^p = Σ_α γ̇^α P^α
    * Plastic spin:             W^p = Σ_α γ̇^α Q^α
    * where γ̇^α are slip system shearing rates from kinetics.
    * 
    * Coordinate convention:
    * All vectors defined in crystal frame where elastic constitutive law applies.
    * Rotations to sample frame handled externally via quaternions.
    * 
    * Template parameters:
    * @tparam num_slip Number of slip systems for this crystal structure
    * @tparam num_extra Number of additional slip-related quantities (default 0)
    * 
    * @ingroup ECMech_slip_geometry
    * 
    * @see ECMech_slipgeom_fcc.h for FCC implementations
    * @see ECMech_slipgeom_bcc.h for BCC implementations
    * @see ECMech_slipgeom_hcp.h for HCP implementations
    */
   template<size_t num_slip, size_t num_extra = 0>
   class SlipGeom {
      public:
         /** @brief Number of slip systems */
         static constexpr int nslip = num_slip;
         /** @brief Number of extra slip-related quantities */
         static constexpr int nSlipExtra = num_extra;

         /** @brief Virtual destructor */
         __ecmech_hdev__
         virtual ~SlipGeom(){}

         /**
          * @brief Access symmetric Schmid tensor components.
          * @return Pointer to P array, shape [ntvec, nslip]
          */
         __ecmech_hdev__ inline virtual const double* getP() const { return m_P_ref_vec; }
         /**
          * @brief Access skew Schmid tensor components.
          * @return Pointer to Q array, shape [nwvec, nslip]
          */
         __ecmech_hdev__ inline virtual const double* getQ() const { return m_Q_ref_vec; }
         /**
          * @brief Access slip plane normals.
          * @return Pointer to m array, shape [nslip, ndim]
          */
         __ecmech_hdev__ inline const double* getM() const { return m_m_ref_vec; }
         /**
          * @brief Access slip directions.
          * @return Pointer to s array, shape [nslip, ndim]
          */
         __ecmech_hdev__ inline const double* getS() const { return m_s_ref_vec; }

         /**
          * @brief Compute stress-dependent Schmid tensors (for non-Schmid models).
          * 
          * Default implementation copies reference values. Overridden by derived classes
          * that implement stress-dependent slip geometry.
          * 
          * @param[out] chia Additional geometric parameters
          * @param[out] P_vec Symmetric Schmid components
          * @param[out] Q_vec Skew Schmid components  
          * @param[in] SvecP Stress state (unused in base class)
          */
         __ecmech_hdev__ inline virtual void getPQ(double* /* chia */, 
                                                   double* P_vec, 
                                                   double* Q_vec, 
                                                   const double* const /* SvecP = nullptr */) const 
         {
             for (int iTvec = 0; iTvec < ecmech::ntvec * nslip; ++iTvec) {
                 P_vec[iTvec] = m_P_ref_vec[iTvec];
             }
             for (int iWvec = 0; iWvec < ecmech::nwvec * nslip; ++iWvec) {
                 Q_vec[iWvec] = m_Q_ref_vec[iWvec];
             }
         }

         /**
          * @brief Compute resolved shear stresses on slip systems.
          * 
          * Projects stress tensor onto slip systems using Schmid tensors:
          *   τ^α = σ : P^α = Σ_ij σ_ij P^α_ij
          * 
          * @param[out] taua Resolved shear stress per system, length nslip
          * @param[in] kirchoff Kirchhoff stress (deviatoric 5-vector)
          * @param[in] P_vec Symmetric Schmid components
          */
         __ecmech_hdev__ inline virtual void evalRSS(double* taua, 
                                                     const double* const kirchoff, 
                                                     const double* P_vec) const
         {
             // resolve stress onto slip systems
             vecsVaTM<ecmech::ntvec, nslip>(taua, kirchoff, P_vec);
         }
       
      protected:
         /** @brief Slip plane normals for all systems, shape [nslip, ndim] */
         double m_m_ref_vec[ ecmech::ndim * nslip ];
         /** @brief Slip directions for all systems, shape [nslip, ndim] */
         double m_s_ref_vec[ ecmech::ndim * nslip ];
         /** @brief Symmetric Schmid components for all systems, shape [nslip, ntvec] */
         double m_P_ref_vec[ ecmech::ntvec * nslip ];
         /** @brief Skew Schmid components for all systems, shape [nslip, nwvec] */
         double m_Q_ref_vec[ ecmech::nwvec * nslip ];
   };
}