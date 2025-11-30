/**
 * @file ECMech_slipgeom_bcc.h
 * @brief Slip geometry for body-centered cubic (BCC) crystals.
 * 
 * This file defines slip system geometries for BCC crystal structures, which
 * can activate multiple slip families depending on temperature and stress state.
 * BCC metals include iron, chromium, tungsten, molybdenum, tantalum, and ferritic steels.
 * 
 * **BCC slip system families**:
 * 
 * 1. **{110}<111> systems** (12 systems, base, always active):
 *    - Primary slip family in BCC metals
 *    - Slip planes: {110} family (6 unique planes)
 *    - Slip directions: <111> family (2 directions per plane)
 *    - Most commonly active at all temperatures
 *    - Often the only active family at low to moderate stress
 * 
 * 2. **{112}<111> systems** (12 systems, optional):
 *    - Additional slip family active at higher stresses/temperatures
 *    - Slip planes: {112} family
 *    - Slip directions: <111> family
 *    - Systems 12-23 (when nslip >= 24)
 *    - Important for polycrystal deformation and texture evolution
 * 
 * 3. **{123}<111> systems** (24 systems, optional):
 *    - High temperature or very high stress slip family
 *    - Slip planes: {123} family
 *    - Slip directions: <111> family
 *    - Systems 24-47 (when nslip = 48)
 *    - Rarely dominant but contributes to overall ductility
 * 
 * **Implementation variants**:
 * 
 * - **SlipGeomBCC<nslip>**: Template class for Schmid-based BCC
 *   - nslip = 12: Base {110}<111> systems only
 *   - nslip = 24: Base + {112}<111> pencil glide
 *   - nslip = 48: All three families
 *   - dynamic = false (fixed Schmid factors)
 *   - nParams = 0 (no parameters)
 * 
 * - **SlipGeomBCCNonSchmid**: Non-Schmid effects for BCC plasticity
 *   - nslip = 12 (base systems only)
 *   - dynamic = true (stress-state-dependent slip)
 *   - nParams = 3 (omega parameters for non-Schmid terms)
 *   - Captures twinning/anti-twinning asymmetry
 *   - Important for low-temperature BCC behavior
 * 
 * - **SlipGeomBCCPencil**: Combined pencil glide with non-Schmid effects
 *   - Similar to BCC24 but with state-dependent geometry
 *   - Uses chi angle to modulate slip resistance
 *   - Accounts for loading direction effects on slip
 * 
 * **Physical behavior**:
 * - Strong temperature dependence of yield stress (Peierls barrier)
 * - Twinning/anti-twinning asymmetry at low temperatures
 * - Non-Schmid effects due to core structure of screw dislocations
 * - Multiple slip families enable high-temperature ductility
 * 
 * **Normalization**:
 * - {110} planes: 1/√2 × [0, ±1, ±1] and permutations
 * - {112} planes: 1/√6 × [±1, ±2, ±2] and permutations (factors of √6)
 * - {123} planes: 1/√14 × [±1, ±2, ±3] and permutations
 * - <111> directions: 1/√3 × [±1, ±1, ±1]
 * 
 * **Variable naming**:
 * - nslipAddBase = 12: Base {110}<111> systems
 * - nslipAddPGa = 12: {112}<111> systems
 * - nslipAddPGb = 24: {123}<111> systems
 * 
 * **Typical kinetics models for BCC**:
 * - BCCMD: Dislocation-density model with thermal activation
 * - KMBalD: Temperature-dependent hardening
 * - Rate-dependent models capturing Peierls stress
 * 
 * **Non-Schmid effects**:
 * For BCC metals, especially at low temperatures:
 * - Resolved shear stress modified by stress-state
 * - Twinning vs anti-twinning directions have different CRSS
 * - Captured through omega parameters in SlipGeomBCCNonSchmid
 * 
 * @see SlipGeom for base class interface
 * @see SlipGeomBCCNonSchmid for non-Schmid BCC plasticity
 * @see KineticsBCCMD for BCC kinetics model
 */

#pragma once

#include "ECMech_slipgeom_base.h"

namespace ecmech {

   /**
    * @brief Body-centered cubic (BCC) slip system geometry with configurable slip families.
    * 
    * SlipGeomBCC implements crystallographic slip systems for BCC crystal structures with
    * support for multiple slip families of varying activity levels. BCC metals exhibit more
    * complex slip behavior than FCC due to less clearly defined slip planes and temperature-
    * dependent slip system activation.
    * 
    * Slip system families (activated progressively):
    * 
    * 1. **Base {110}<111> systems** (12 systems, always included):
    *    - Most commonly active at low to moderate temperatures
    *    - Slip planes: {110} family
    *    - Slip directions: <111> family
    *    - Systems 0-11
    * 
    * 2. **{112}<111> systems** (12 systems, optional):
    *    - Additional slip family active at higher stresses/temperatures
    *    - Slip planes: {112} family
    *    - Slip directions: <111> family  
    *    - Systems 12-23 (when nslip >= 24)
    * 
    * 3. **{123}<111> systems** (24 systems, optional):
    *    - High temperature or very high stress slip family
    *    - Slip planes: {123} family
    *    - Slip directions: <111> family
    *    - Systems 24-47 (when nslip = 48)
    * 
    * Template parameter nSlipTmplt:
    * Determines which slip families are included:
    * - nSlipTmplt = 12: Base {110}<111> only
    * - nSlipTmplt = 24: Base + {112}<111>
    * - nSlipTmplt = 48: Base + {112}<111> + {123}<111>
    * 
    * Coordinate convention:
    * Vectors defined in crystal frame with:
    * - x₁, x₂, x₃ aligned with cubic cell edges [100], [010], [001]
    * - All normals and directions in this basis
    * 
    * Parameters:
    * No adjustable parameters - slip geometry is purely crystallographic.
    * 
    * @tparam nSlipTmplt Number of slip systems (must be 12, 24, or 48)
    * 
    * @ingroup ECMech_slip_geometry
    * 
    * @see SlipGeomFCC for face-centered cubic geometry
    * @see SlipGeomBCCNonSchmid for non-Schmid stress-dependent variant
    * @see SlipGeomBCCPencil for pencil glide variant
    */
   template<int nSlipTmplt>
   class SlipGeomBCC : public SlipGeom<nSlipTmplt>
   {
      private:
         /** @brief Number of {110}<111> base systems */
         static constexpr int nslipAddBase = 12;
         /** @brief Number of {112}<111> systems */
         static constexpr int nslipAddPGa = 12;
         /** @brief Number of {123}<111> systems */
         static constexpr int nslipAddPGb = 24;

      public:
         /** @brief Slip geometry does not depend on state */
         static const bool dynamic = false;
         /** @brief Total number of slip systems */
         static constexpr int nslip = nSlipTmplt;
         /** @brief Number of parameters required */
         static constexpr int nParams = 0;

         /** @brief Total systems with base only */
         static constexpr int nslipBase = nslipAddBase;
         /** @brief Total systems with base + {112}<111> systems */
         static constexpr int nslipPGa = nslipAddBase + nslipAddPGa;
         /** @brief Total systems with all families */
         static constexpr int nslipPGb = nslipAddBase + nslipAddPGa + nslipAddPGb;

         /**
          * @brief Default constructor with compile-time validation.
          * 
          * Asserts that template parameter is one of the supported values.
          */
         __ecmech_hdev__
         SlipGeomBCC() {
            assert(nslip == nslipBase || nslip == nslipPGa || nslip == nslipPGb);
         }
         /** @brief Destructor */
         __ecmech_hdev__
         ~SlipGeomBCC() {}

         /**
          * @brief Constructor with parameters.
          * @param params Parameter array (unused for BCC)
          */
         __ecmech_hdev__
         SlipGeomBCC(const double* const params) {
            setParams(params);
         }

         /**
          * @brief Initialize from parameter vector.
          * @param params Parameter vector (empty for BCC)
          */
         __ecmech_host__
         void setParams(const std::vector<double> & params)
         {
            setParams(params.data());
         }


         /**
          * @brief Initialize BCC slip geometry from parameters.
          * 
          * Constructs slip systems for BCC structure based on template parameter nslip.
          * Each slip family is added conditionally based on the total system count.
          * 
          * Slip system construction:
          * 1. **Base {110}<111>** (always included, 12 systems):
          *    - Slip direction: <111> in 4 equivalent variants
          *    - Slip planes: {110} in 3 variants per direction
          * 
          * 2. **{112}<111>** (included if nslip >= 24, adds 12 systems):
          *    - Same <111> directions
          *    - Higher-index {112} slip planes
          * 
          * 3. **{123}<111>** (included if nslip = 48, adds 24 systems):
          *    - Same <111> directions
          *    - Highest-index {123} slip planes
          * 
          * Implementation uses helper lambda add_vec_data() to sequentially
          * populate the master arrays mVecs and sVecs, then computes Schmid
          * tensors via fillFromMS().
          * 
          * @param params Parameter array (unused, BCC geometry is fixed)
          */
         __ecmech_hdev__
         void setParams(const double* const)
         {
            double mVecs[nslipPGb * ecmech::ndim] = {};
            double sVecs[nslipPGb * ecmech::ndim] = {};

            auto add_vec_data = [=] (const double* const array_src, double* const array_dst, const size_t length) {
               for(size_t ivec = 0; ivec < length; ivec++)
               {
                  array_dst[ivec] = array_src[ivec];
               }
            };

            {
               // m = (/ zero, sqr2i, -sqr2i /)
               // s = (/ sqr3i, sqr3i, sqr3i /)
               const int nslipThese = nslipAddBase;
               //
               // do not yet bother with making slip systems from symmetry group -- just write them out
               const double P3 = sqr3i, M3 = -sqr3i;
               const double P2 = sqr2i, M2 = -sqr2i;
               const double Z = zero;
               //#Slip direction CUB111
               const double sVecsThese[ nslip * ecmech::ndim ] = {
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
               //#Slip plane normal CUB110
               const double mVecsThese[ nslip * ecmech::ndim ] = {
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

               add_vec_data(mVecsThese, mVecs, nslipThese * ecmech::ndim);
               add_vec_data(sVecsThese, sVecs, nslipThese * ecmech::ndim);
            }

            if (nslip >= nslipPGa) {
               const double twSqr6i = 2.0 * sqr6i;

               // 12 {112}<111> slip systems
               const int nslipThese = nslipAddPGa;

               const double mVecsThese[ nslipThese * ecmech::ndim ] = {
                  -twSqr6i, sqr6i, sqr6i,
                  sqr6i, -twSqr6i, sqr6i,
                  sqr6i, sqr6i, -twSqr6i,
                  -sqr6i, -twSqr6i, sqr6i,
                  twSqr6i, sqr6i, sqr6i,
                  -sqr6i, sqr6i, -twSqr6i,
                  twSqr6i, -sqr6i, sqr6i,
                  -sqr6i, twSqr6i, sqr6i,
                  -sqr6i, -sqr6i, -twSqr6i,
                  sqr6i, twSqr6i, sqr6i,
                  -twSqr6i, -sqr6i, sqr6i,
                  sqr6i, -sqr6i, -twSqr6i,
               };
               const double sVecsThese[ nslipThese * ecmech::ndim ] = {
                  sqr3i, sqr3i, sqr3i,
                  sqr3i, sqr3i, sqr3i,
                  sqr3i, sqr3i, sqr3i,
                  -sqr3i, sqr3i, sqr3i,
                  -sqr3i, sqr3i, sqr3i,
                  -sqr3i, sqr3i, sqr3i,
                  -sqr3i, -sqr3i, sqr3i,
                  -sqr3i, -sqr3i, sqr3i,
                  -sqr3i, -sqr3i, sqr3i,
                  sqr3i, -sqr3i, sqr3i,
                  sqr3i, -sqr3i, sqr3i,
                  sqr3i, -sqr3i, sqr3i,
               };
               add_vec_data(mVecsThese, &mVecs[nslipAddBase], nslipThese * ecmech::ndim);
               add_vec_data(sVecsThese, &sVecs[nslipAddBase], nslipThese * ecmech::ndim);
            }

            if (nslip >= nslipPGb) {
               const double mPg2a = 1.0 / sqrt(14.0);
               const double mPg2b = 2.0 / sqrt(14.0);
               const double mPg2c = 3.0 / sqrt(14.0);

               // 24 {123}<111> slip systems
               const int nslipThese = nslipAddPGb;

               const double mVecsThese[ nslipThese * ecmech::ndim ] = {
                  mPg2c, -mPg2a, -mPg2b,
                  -mPg2b, mPg2c, -mPg2a,
                  -mPg2a, -mPg2b, mPg2c,
                  mPg2a, mPg2c, -mPg2b,
                  -mPg2c, -mPg2b, -mPg2a,
                  mPg2b, -mPg2a, mPg2c,
                  -mPg2c, mPg2a, -mPg2b,
                  mPg2b, -mPg2c, -mPg2a,
                  mPg2a, mPg2b, mPg2c,
                  -mPg2a, -mPg2c, -mPg2b,
                  mPg2c, mPg2b, -mPg2a,
                  -mPg2b, mPg2a, mPg2c,
                  -mPg2a, mPg2c, mPg2b,
                  mPg2c, -mPg2b, mPg2a,
                  -mPg2b, -mPg2a, -mPg2c,
                  -mPg2c, -mPg2a, mPg2b,
                  mPg2b, mPg2c, mPg2a,
                  mPg2a, -mPg2b, -mPg2c,
                  mPg2a, -mPg2c, mPg2b,
                  -mPg2c, mPg2b, mPg2a,
                  mPg2b, mPg2a, -mPg2c,
                  mPg2c, mPg2a, mPg2b,
                  -mPg2b, -mPg2c, mPg2a,
                  -mPg2a, mPg2b, -mPg2c,
               };
               const double sVecsThese[ nslipThese * ecmech::ndim ] = {
                  sqr3i, sqr3i, sqr3i,
                  sqr3i, sqr3i, sqr3i,
                  sqr3i, sqr3i, sqr3i,
                  -sqr3i, sqr3i, sqr3i,
                  -sqr3i, sqr3i, sqr3i,
                  -sqr3i, sqr3i, sqr3i,
                  -sqr3i, -sqr3i, sqr3i,
                  -sqr3i, -sqr3i, sqr3i,
                  -sqr3i, -sqr3i, sqr3i,
                  sqr3i, -sqr3i, sqr3i,
                  sqr3i, -sqr3i, sqr3i,
                  sqr3i, -sqr3i, sqr3i,
                  sqr3i, sqr3i, -sqr3i,
                  sqr3i, sqr3i, -sqr3i,
                  sqr3i, sqr3i, -sqr3i,
                  -sqr3i, sqr3i, -sqr3i,
                  -sqr3i, sqr3i, -sqr3i,
                  -sqr3i, sqr3i, -sqr3i,
                  -sqr3i, -sqr3i, -sqr3i,
                  -sqr3i, -sqr3i, -sqr3i,
                  -sqr3i, -sqr3i, -sqr3i,
                  sqr3i, -sqr3i, -sqr3i,
                  sqr3i, -sqr3i, -sqr3i,
                  sqr3i, -sqr3i, -sqr3i,
               };
               add_vec_data(mVecsThese, &mVecs[nslipAddPGa], nslipThese * ecmech::ndim);
               add_vec_data(sVecsThese, &sVecs[nslipAddPGa], nslipThese * ecmech::ndim);
            }

            fillFromMS(this->m_P_ref_vec, this->m_Q_ref_vec, mVecs, sVecs, this->nslip);

            for (int i = 0; i < nslip * ecmech::ndim; i++) {
               this->m_s_ref_vec[i] = sVecs[i];
               this->m_m_ref_vec[i] = mVecs[i];
            }
         }


         /**
          * @brief Retrieve slip geometry parameters.
          * @param params Parameter vector (empty for BCC)
          */
         __ecmech_host__
         void getParams(std::vector<double> & /* params */
                        ) const {
            // do not clear params in case adding to an existing set
         }

   }; // SlipGeomBCC

   /**
    * @brief BCC pencil glide slip geometry with stress-dependent slip plane rotation.
    * 
    * SlipGeomBCCPencil implements a simplified BCC slip geometry with 4 primary {111}<111> slip
    * systems plus 4 additional "extra" quantities for stress-dependent pencil glide behavior.
    * This model captures the phenomenon where BCC slip planes can rotate within the pencil glide
    * zone formed by the intersection of {110} and {112} planes.
    * 
    * Pencil glide mechanism:
    * In BCC crystals at certain conditions, the slip plane can rotate continuously between
    * crystallographic planes while maintaining the <111> slip direction. The "pencil" is the
    * zone of possible slip planes containing the fixed slip direction.
    * 
    * Slip systems (4 primary):
    * - Slip directions: Four {111} directions (one per octant)
    * - Slip planes: Initially {110} planes, but rotate based on stress state
    * - Additional 4 "extra" quantities track pencil glide parameters (nSlipExtra = 4)
    * 
    * Stress-dependent behavior:
    * Unlike standard BCC, slip plane normals are recomputed each time based on:
    * 1. Current stress state (Cauchy stress tensor)
    * 2. Peach-Koehler force on dislocation
    * 3. Maximum resolved shear stress principle (MRSSP)
    * 4. Chi angle (orientation within pencil glide zone)
    * 
    * Dynamic flag:
    * Set to true because slip geometry changes with stress state, requiring
    * recomputation of P and Q tensors during analysis.
    * 
    * Coordinate convention:
    * - Crystal frame aligned with cubic cell edges
    * - Slip directions fixed in <111> family
    * - Slip plane normals computed dynamically
    * 
    * Parameters:
    * No adjustable parameters - geometry computed from stress state.
    * 
    * @ingroup ECMech_slip_geometry
    * 
    * @see SlipGeomBCC for standard BCC slip geometry
    * @see SlipGeomBCCNonSchmid for non-Schmid stress effects
    */
   class SlipGeomBCCPencil : public SlipGeom<4, 4>
   {
      public:
         /** @brief Slip geometry depends on stress state */
         static const bool dynamic = true;
         /** @brief Number of parameters required */
         static constexpr int nParams = 0;

         /** @brief Default constructor */
         SlipGeomBCCPencil() = default;
         /** @brief Destructor */
         __ecmech_hdev__
         ~SlipGeomBCCPencil() {}

         /**
          * @brief Constructor with parameters.
          * @param params Parameter array (unused)
          */
         __ecmech_hdev__
         SlipGeomBCCPencil(const double* const params) {
            setParams(params);
         }

         /**
          * @brief Initialize from parameter vector.
          * @param params Parameter vector (empty)
          */
         __ecmech_host__
         void setParams(const std::vector<double> & params)
         {
            setParams(params.data());
         }

         /**
          * @brief Initialize BCC pencil glide slip geometry.
          * 
          * Sets up 4 base {111} slip directions with initial {110} slip plane normals.
          * Actual slip planes will be recomputed dynamically based on stress state via getPQ().
          * 
          * Slip system enumeration:
          * - System 0: [111] direction
          * - System 1: [-111] direction  
          * - System 2: [1-11] direction
          * - System 3: [11-1] direction
          * 
          * Initial planes are {110} family, but these serve only as reference - actual
          * active planes determined by stress-dependent pencil glide mechanism.
          * 
          * @param params Parameter array (unused)
          */
         __ecmech_hdev__
         void setParams(const double* const)
         {
            // s = (/ sqr3i, sqr3i, sqr3i /)
            //
            const double P3 = sqr3i, M3 = -sqr3i;
            const double P2 = sqr2i, M2 = -sqr2i;
            const double Z = zero;
         
            const double sVecs[ nslip * ecmech::ndim ] = {
               P3, P3, P3,
               M3, P3, P3,
               P3, M3, P3,
               P3, P3, M3};
               
            const double mVecs[ nslip * ecmech::ndim ] = {
               Z, M2, P2,
               Z, M2, P2,
               Z, P2, P2,
               Z, P2, P2};

            fillFromMS(this->m_P_ref_vec, this->m_Q_ref_vec,
                       mVecs, sVecs, this->nslip);

            for (int i = 0; i < nslip * ecmech::ndim; i++) {
               m_s_ref_vec[i] = sVecs[i];
               m_m_ref_vec[i] = mVecs[i];
            }
         }

         /**
          * @brief Retrieve slip geometry parameters.
          * @param params Parameter vector (empty)
          */
         __ecmech_host__
         void getParams(std::vector<double> & /* params */
                        ) const {
            // do not clear params in case adding to an existing set
         }

         /**
          * @brief Compute stress-dependent Schmid tensors via pencil glide mechanism.
          * 
          * Determines optimal slip plane orientation for each <111> slip direction based on
          * current stress state. This implements the pencil glide model where slip planes
          * rotate to maximize resolved shear stress.
          * 
          * Algorithm:
          * 1. Convert deviatoric stress vector to full Cauchy stress tensor
          * 2. For each slip system:
          *    a. Compute Peach-Koehler force: f_PK = (σ·b) × b
          *    b. If |f_PK| > ε: slip plane normal m = b × f_PK (perpendicular to force)
          *    c. If |f_PK| ≤ ε: use reference plane normal (no rotation)
          *    d. Compute MRSSP angle χ (orientation in pencil glide zone)
          *    e. Fold χ into twinning/anti-twinning (T/AT) region [-30°, 30°]
          * 3. Compute Schmid tensors P, Q from updated (m, s) pairs
          * 4. Store χ angles in chia array for extra slip quantities
          * 
          * Twinning/anti-twinning asymmetry:
          * Chi angle is folded to [-π/6, π/6] to capture T/AT directional preference
          * in BCC slip, where slip resistance depends on sense of shear.
          * 
          * @param[out] chia Chi angles (MRSSP orientation) for nslip systems
          * @param[out] P_vec Updated symmetric Schmid components
          * @param[out] Q_vec Updated skew Schmid components
          * @param[in] SvecP Stress state (deviatoric 6-vector + pressure)
          */
         __ecmech_hdev__ inline void getPQ(double* chia, 
                                           double* P_vec, 
                                           double* Q_vec, 
                                           const double* const SvecP) const override final
         {
             double eps = 1e-10;
             double mVecs[nslip * ecmech::ndim];
             
             double S[ecmech::ndim * ecmech::ndim];
             // Svec: 11' 22' 33' 23 31 12 p
             S[ECMECH_NN_INDX(0, 0, 3)] = SvecP[0] + ecmech::onethird * SvecP[6];
             S[ECMECH_NN_INDX(1, 1, 3)] = SvecP[1] + ecmech::onethird * SvecP[6];
             S[ECMECH_NN_INDX(2, 2, 3)] = SvecP[2] + ecmech::onethird * SvecP[6];
             S[ECMECH_NN_INDX(1, 2, 3)] = S[ECMECH_NN_INDX(2, 1, 3)] = SvecP[3];
             S[ECMECH_NN_INDX(2, 0, 3)] = S[ECMECH_NN_INDX(0, 2, 3)] = SvecP[4];
             S[ECMECH_NN_INDX(0, 1, 3)] = S[ECMECH_NN_INDX(1, 0, 3)] = SvecP[5];
             
             for (int iSlip = 0; iSlip < nslip; ++iSlip) {
                 const double* sVec = &m_s_ref_vec[iSlip * ecmech::ndim];
                 
                 // PK force direction
                 double fpk[ecmech::ndim] = {0.0};
                 double Sb[ecmech::ndim];
                 vecsVMa<ecmech::ndim>(Sb, S, sVec);
                 if (vecNorm<ecmech::ndim>(Sb) > eps) {
                     vecCrossProd(fpk, Sb, sVec);
                 }
                 
                 // Normal direction
                 double* mVec = &mVecs[iSlip * ecmech::ndim];
                 if (vecNorm<ecmech::ndim>(fpk) > eps) {
                     vecCrossProd(mVec, sVec, fpk);
                     vecsVNormalize<ecmech::ndim>(mVec);
                 } else {
                     for (int i = 0; i < ecmech::ndim; i++)
                        mVec[i] = m_m_ref_vec[iSlip * ecmech::ndim + i];
                 }
                 
                 // MRSSP angle
                 double n0Vec[ecmech::ndim] = { //n0 = 1/sqrt(2)*(2*b[0],-b[1],-b[2])
                      2.0*sqr2i*sVec[0],
                     -1.0*sqr2i*sVec[1],
                     -1.0*sqr2i*sVec[2]
                 };
                 double t0Vec[ecmech::ndim];
                 vecCrossProd(t0Vec, n0Vec, sVec);
                 double fx = vecsyadotb<ecmech::ndim>(fpk, t0Vec);
                 double fy = vecsyadotb<ecmech::ndim>(fpk, n0Vec);
                 double chi = atan2(fy, fx)-M_PI/6.0;
                 // Fold into T/AT primary region (-30:30)
                 if (chi > 1.0*M_PI/6.0 && chi <= 3.0*M_PI/6.0) {
                     chi = M_PI/3.0-chi;
                 } else if (chi > 3.0*M_PI/6.0 && chi <= 5.0*M_PI/6.0) {
                     chi -= 2.0*M_PI/3.0;
                 } else if (chi >= -7.0*M_PI/6.0 && chi < -5.0*M_PI/6.0) {
                     chi = -M_PI-chi;
                 } else if (chi >= -5.0*M_PI/6.0 && chi < -3.0*M_PI/6.0) {
                     chi += 2.0*M_PI/3.0;
                 } else if (chi >= -3.0*M_PI/6.0 && chi < -1.0*M_PI/6.0) {
                     chi = -M_PI/3.0-chi;
                 }
                 chia[iSlip] = chi;
                 m_chia[iSlip] = chia[iSlip];
             }
             
             fillFromMS(P_vec, Q_vec, mVecs, m_s_ref_vec, nslip);
             fillFromMS(m_P_vec, m_Q_vec, mVecs, m_s_ref_vec, nslip);
         }
         /**
          * @brief Access updated symmetric Schmid components.
          * @return Pointer to P array (stress-dependent)
          */
         __ecmech_hdev__ inline virtual const double* getP() const override final { return m_P_vec; }
         /**
          * @brief Access updated skew Schmid components.
          * @return Pointer to Q array (stress-dependent)
          */
         __ecmech_hdev__ inline virtual const double* getQ() const override final { return m_Q_vec; }
         /**
          * @brief Retrieve extra slip quantities (chi angles).
          * @param[out] chia Array to receive chi angles, length nslip
          */
         __ecmech_hdev__
         inline
         void getExtras(double* const chia) const
         {
            for (size_t islip = 0; islip < nslip; islip++) {
               chia[islip] = m_chia[islip];
            }
         }
      private:
         /** @brief Stress-dependent symmetric Schmid components */
         mutable double m_P_vec[ ecmech::ntvec * nslip ];
         /** @brief Stress-dependent skew Schmid components */
         mutable double m_Q_vec[ ecmech::nwvec * nslip ];
         /** @brief Chi angles (MRSSP orientation in pencil glide zone) */
         mutable double m_chia[ nslip ];

   }; // SlipGeomBCCPencil
   
   /**
    * @brief BCC slip geometry with non-Schmid stress effects.
    * 
    * SlipGeomBCCNonSchmid implements 12 {110}<111> BCC slip systems with non-Schmid contributions
    * to resolved shear stress. Non-Schmid effects capture deviations from classical Schmid law where
    * stress components not aligned with the slip system still influence slip resistance.
    * 
    * Non-Schmid behavior in BCC:
    * BCC metals (especially at low temperatures) exhibit asymmetric slip behavior where the critical
    * resolved shear stress depends not just on τ_RSS = σ:P but also on stress components perpendicular
    * to the slip plane and parallel to the twinning direction.
    * 
    * Modified resolved shear stress:
    *   τ_eff = τ_Schmid + ω₁·τ₁ + ω₂·τ₂ + ω₃·τ₃
    * where:
    * - τ_Schmid: Classical Schmid stress (σ:P)
    * - τ₁: Non-Schmid component 1 (rotated plane, same slip direction)
    * - τ₂: Non-Schmid component 2 (same plane, perpendicular direction)
    * - τ₃: Non-Schmid component 3 (combined rotation)
    * - ω₁, ω₂, ω₃: Weighting parameters (material-dependent)
    * 
    * Slip direction sense selection:
    * For each slip system, both +b and -b directions are evaluated. The direction producing
    * higher effective resolved shear stress is selected, accounting for twinning/anti-twinning
    * asymmetry inherent in non-Schmid behavior.
    * 
    * Slip systems:
    * 12 {110}<111> systems identical to SlipGeomBCC<12>, but with stress-dependent
    * Schmid tensors recomputed based on non-Schmid contributions.
    * 
    * Parameters (3 required):
    * - ω₁: Weight for rotated plane contribution
    * - ω₂: Weight for cross-slip direction contribution  
    * - ω₃: Weight for combined rotation contribution
    * 
    * When |ω₁| + |ω₂| + |ω₃| < ε, model reverts to classical Schmid behavior for efficiency.
    * 
    * Dynamic flag:
    * Set to true because effective slip geometry depends on stress state through
    * non-Schmid projections and slip direction selection.
    * 
    * @ingroup ECMech_slip_geometry
    * 
    * @see SlipGeomBCC for standard Schmid-law BCC geometry
    * @see SlipGeomBCCPencil for pencil glide model
    */
   class SlipGeomBCCNonSchmid : public SlipGeom<12>
   {
      public:
         /** @brief Slip geometry depends on stress state */
         static const bool dynamic = true;
         /** @brief Number of parameters (3 omega weights) */
         static constexpr int nParams = 3;
         /** @brief No extra slip quantities */
         static constexpr size_t nSlipExtra = 0;

         /** @brief Default constructor */
         SlipGeomBCCNonSchmid() = default;
         /** @brief Destructor */
         __ecmech_hdev__
         ~SlipGeomBCCNonSchmid() {}


         /**
          * @brief Constructor with parameters.
          * @param params Parameter array [ω₁, ω₂, ω₃]
          */
         __ecmech_hdev__
         SlipGeomBCCNonSchmid(const double* const params) {
            setParams(params);
         }

         /**
          * @brief Initialize from parameter vector.
          * @param params Parameter vector containing [ω₁, ω₂, ω₃]
          */
         __ecmech_host__
         void setParams(const std::vector<double> & params)
         {
            setParams(params.data());
         }

         /**
          * @brief Initialize non-Schmid BCC slip geometry.
          * 
          * Sets up 12 base {110}<111> slip systems and extracts non-Schmid weighting
          * parameters. Checks if parameters are effectively zero to enable isotropic
          * (classical Schmid) optimization.
          * 
          * Slip plane construction:
          * Planes are chosen to be consistent with twinning/anti-twinning (T/AT) directions,
          * ensuring proper asymmetry in non-Schmid response.
          * 
          * Isotropic check:
          * If sum(|ω_i|) < √ε, sets m_isotropic = true to bypass non-Schmid calculations
          * and use classical Schmid projections for computational efficiency.
          * 
          * @param params Parameter array: [0]=ω₁, [1]=ω₂, [2]=ω₃
          */
         __ecmech_hdev__
         void setParams(const double* const params)
         {
            const double* parsIt = params;

            m_omegas[0] = *parsIt; ++parsIt;
            m_omegas[1] = *parsIt; ++parsIt;
            m_omegas[2] = *parsIt; ++parsIt;
            const double sum_omegas = fabs(m_omegas[0]) + fabs(m_omegas[1]) + fabs(m_omegas[2]);
            m_isotropic = sum_omegas < ecmech::idp_eps_sqrt;

            const double P3 = sqr3i, M3 = -sqr3i;
            const double P2 = sqr2i, M2 = -sqr2i;
            const double Z = zero;
            
            const double sVecs[ nslip * ecmech::ndim ] = {
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
           
            // This list of planes has been generated to be 
            // consistent with T/AT directions
            const double mVecs[ nslip * ecmech::ndim ] = {
                 P2, M2,  Z,
                 Z,  P2, M2,
                 M2,  Z, P2,
                 P2,  Z, P2,
                 M2, P2,  Z,
                 Z,  M2, M2,
                 P2,  Z, M2,
                 M2, M2,  Z,
                 Z,  P2, P2,
                 P2, P2,  Z,
                 Z,  M2, P2,
                 M2,  Z, M2};
            
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
          * @brief Retrieve non-Schmid parameters.
          * @param params Parameter vector to receive [ω₁, ω₂, ω₃]
          */
         __ecmech_host__
         void getParams(std::vector<double> & params
                        ) const {
#ifdef ECMECH_DEBUG
            // do not clear params in case adding to an existing set
            int paramsStart = params.size();
#endif
            params.push_back(m_omegas[0]);
            params.push_back(m_omegas[1]);
            params.push_back(m_omegas[2]);

#ifdef ECMECH_DEBUG
            assert((params.size() - paramsStart) == nParams);
#endif
         }

         /**
          * @brief Compute non-Schmid resolved shear stress projection.
          * 
          * For each slip system, evaluates both +b and -b slip directions with non-Schmid
          * contributions, selecting the direction with higher effective RSS. This captures
          * T/AT asymmetry in BCC slip.
          * 
          * Algorithm for each system:
          * 1. Loop over slip direction senses: s and -s
          * 2. For each sense:
          *    a. Schmid contribution: τ_S = σ : P(m,s)
          *    b. Non-Schmid 1: τ_NS1 = ω₁ · σ : P(m',s) where m' rotated by 30°
          *    c. Non-Schmid 2: τ_NS2 = ω₂ · σ : P(m,s×m) (cross-slip direction)
          *    d. Non-Schmid 3: τ_NS3 = ω₃ · σ : P(m',(s×m')) (combined)
          *    e. Total: τ_eff = τ_S + τ_NS1 - τ_NS2 + τ_NS3
          * 3. Select sense with max(τ_eff)
          * 4. Compute and store P, Q for selected sense
          * 
          * Rotation for non-Schmid components:
          * m' = 0.5·m + √3/2·(s×m) represents 30° rotation about slip direction,
          * corresponding to geometric relationship between {110} and {112} planes.
          * 
          * @param[out] taua Effective resolved shear stress per system
          * @param[out] P_vec Symmetric Schmid components (if fill_PQ=true)
          * @param[out] Q_vec Skew Schmid components (if fill_PQ=true)
          * @param[in] kirchoff Kirchhoff stress (deviatoric 5-vector)
          * @param[in] fill_PQ Whether to fill P_vec and Q_vec arrays
          */
         __ecmech_hdev__ inline void NSprojection(double* taua,
                                                  double* P_vec, 
                                                  double* Q_vec, 
                                                  const double* const kirchoff,
                                                  bool fill_PQ) const
         {
             // Resolve stress onto slip systems
             // Compute RSS considering both senses of the slip direction
             // and select the most favorable direction
             for (int iSlip = 0; iSlip < nslip; ++iSlip) {
                 double P_s[2 * ecmech::ntvec];
                 double Q_s[2 * ecmech::nwvec];
                 double tau_s[2] = { 0.0 };
                 
                 for (int iS = 0; iS < 2; ++iS) {
                     
                     double P_tmp[ecmech::ntvec];
                     double Q_tmp[ecmech::nwvec];
                     
                     const double* mVec = &m_m_ref_vec[iSlip * ecmech::ndim];
                     double sVec[ecmech::ndim];
                     for (int i = 0; i < ecmech::ndim; i++)
                         sVec[i] = (1 - 2*iS) * m_s_ref_vec[iSlip * ecmech::ndim + i];
                         
                     // Schmid
                     fillFromMS(P_tmp, Q_tmp, mVec, sVec, 1);
                     for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                         P_s[iS * ecmech::ntvec + iTvec] = P_tmp[iTvec];
                         tau_s[iS] += kirchoff[iTvec] * P_tmp[iTvec];
                     }
                     for (int iWvec = 0; iWvec < ecmech::nwvec; ++iWvec) {
                         Q_s[iS * ecmech::nwvec + iWvec] = Q_tmp[iWvec];
                     }
                     
                     // Non-Schmid
                     double smVec[ecmech::ndim];
                     double mpVec[ecmech::ndim];
                     vecCrossProd(smVec, sVec, mVec);
                     for (int i = 0; i < ecmech::ndim; i++)
                         mpVec[i] = 0.5*mVec[i] + 0.8660254037844386*smVec[i];
                         
                     // Omega 1
                     fillFromMS(P_tmp, Q_tmp, mpVec, sVec, 1);
                     for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                         tau_s[iS] += m_omegas[0] * kirchoff[iTvec] * P_tmp[iTvec];
                     }
                     
                     // Omega 2
                     fillFromMS(P_tmp, Q_tmp, mVec, smVec, 1);
                     for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                         tau_s[iS] -= m_omegas[1] * kirchoff[iTvec] * P_tmp[iTvec];
                     }
                     
                     // Omega 3
                     double mpsVec[ecmech::ndim];
                     vecCrossProd(mpsVec, mpVec, sVec);
                     fillFromMS(P_tmp, Q_tmp, mpVec, mpsVec, 1);
                     for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                         tau_s[iS] += m_omegas[2] * kirchoff[iTvec] * P_tmp[iTvec];
                     }
                 }
                 
                 // Keep highest value
                 int iS = (int)(tau_s[1] > tau_s[0]);
                 
                 taua[iSlip] = tau_s[iS];
                 if (taua[iSlip] < 0.0) taua[iSlip] = 0.0;
                 
                 if (fill_PQ) {
                     for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                         P_vec[ECMECH_NM_INDX(iTvec, iSlip, ecmech::ntvec, nslip)] = P_s[iS * ecmech::ntvec + iTvec];
                         m_P_vec[ECMECH_NM_INDX(iTvec, iSlip, ecmech::ntvec, nslip)] = P_s[iS * ecmech::ntvec + iTvec];
                     }
                     for (int iWvec = 0; iWvec < ecmech::nwvec; ++iWvec) {
                         Q_vec[ECMECH_NM_INDX(iWvec, iSlip, ecmech::nwvec, nslip)] = Q_s[iS * ecmech::nwvec + iWvec];
                         m_Q_vec[ECMECH_NM_INDX(iWvec, iSlip, ecmech::nwvec, nslip)] = Q_s[iS * ecmech::nwvec + iWvec];
                     }
                 }
             }
         }

         /**
          * @brief Compute stress-dependent Schmid tensors with non-Schmid effects.
          * 
          * If isotropic (all ω ≈ 0), returns reference Schmid tensors.
          * Otherwise, computes full non-Schmid projection from stress state.
          * 
          * @param[out] chia Unused for this class
          * @param[out] P_vec Symmetric Schmid components
          * @param[out] Q_vec Skew Schmid components
          * @param[in] SvecP Stress state for non-Schmid evaluation
          */
         __ecmech_hdev__ inline void getPQ(double* /*chia*/,
                                           double* P_vec, 
                                           double* Q_vec, 
                                           const double* const SvecP) const override final
         {
            if (m_isotropic) {
                  for (int iTvec = 0; iTvec < ecmech::ntvec * nslip; ++iTvec) {
                        P_vec[iTvec] = m_P_ref_vec[iTvec];
                        m_P_vec[iTvec] = m_P_ref_vec[iTvec];
                  }
                  for (int iWvec = 0; iWvec < ecmech::nwvec * nslip; ++iWvec) {
                        Q_vec[iWvec] = m_Q_ref_vec[iWvec];
                        m_Q_vec[iWvec] = m_Q_ref_vec[iWvec];
                  }
            } else {
               // we need to reverse the stress first...
               double kirchoff[ecmech::nsvec];
               kirchoff[iSvecS] = -sqr3 * SvecP[iSvecP];
               kirchoff[0] = sqr2i * SvecP[0] - sqr2i * SvecP[1];
               kirchoff[1] = - sqr3b2 * SvecP[0] - sqr3b2 * SvecP[1];
               kirchoff[4] = sqr2 * SvecP[3]; // 23
               kirchoff[3] = sqr2 * SvecP[4]; // 31
               kirchoff[2] = sqr2 * SvecP[5]; // 12
               
               double taua[nslip];
               NSprojection(taua, P_vec, Q_vec, kirchoff, true);
            }
         }

         /**
          * @brief Compute resolved shear stress with non-Schmid effects.
          * 
          * If isotropic, uses classical Schmid projection.
          * Otherwise, applies full non-Schmid model.
          * 
          * @param[out] taua Resolved shear stress per system
          * @param[in] kirchoff Kirchhoff stress
          * @param[in] P_vec Unused (non-Schmid computes internally)
          */
         __ecmech_hdev__ inline void evalRSS(double* taua, 
                                             const double* const kirchoff, 
                                             const double* /*P_vec*/) const override final
         {
            if (m_isotropic) {
               SlipGeom::evalRSS(taua, kirchoff, m_P_ref_vec);
            } else {
               NSprojection(taua, NULL, NULL, kirchoff, false);
            }
         }

         /**
          * @brief Access stress-dependent symmetric Schmid components.
          * @return Pointer to P array
          */
         __ecmech_hdev__ inline virtual const double* getP() const override final { return m_P_vec; }
         /**
          * @brief Access stress-dependent skew Schmid components.
          * @return Pointer to Q array
          */
         __ecmech_hdev__ inline virtual const double* getQ() const override final { return m_Q_vec; }
         /**
          * @brief No extra slip quantities for this class.
          * @param[out] Unused
          */
         __ecmech_hdev__
         inline
         void getExtras(double* const) const {}
      private:
         /** @brief Stress-dependent symmetric Schmid components */
         mutable double m_P_vec[ ecmech::ntvec * nslip ];
         /** @brief Stress-dependent skew Schmid components */
         mutable double m_Q_vec[ ecmech::nwvec * nslip ];
         /** @brief Non-Schmid weighting parameters [ω₁, ω₂, ω₃] */
         double m_omegas[3];
         /** @brief Flag indicating classical Schmid behavior (all ω ≈ 0) */
         bool m_isotropic = false;

   }; // SlipGeomBCCNonSchmid

}