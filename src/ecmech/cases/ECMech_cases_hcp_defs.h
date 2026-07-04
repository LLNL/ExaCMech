/**
 * @file ECMech_cases_hcp_defs.h
 * @brief Glossary of the concrete HCP (hexagonal close-packed) material models available
 * through the string names accepted by `makeMatModelHCP`/`makeMatModel`.
 *
 * Only one HCP model is currently defined; the corresponding `"evptn_HCP_A"` string name
 * (used in `cases/ECMech_cases_hcp.cxx`) is what a calling code passes to
 * `ecmech::makeMatModel`/`ecmech::makeMatModelHCP` to obtain it.
 */

#pragma once

#include "ECMech_cases.h"

// We currently don't have a ton of HCP use cases but that could change in the future
// Provide a set of defines that users can bring in if they'd like / have a need for them
namespace ecmech {

    /** @brief HCP slip geometry with basal, prismatic, pyramidal-<a>, and pyramidal-<c+a> families (24 slip systems total), parameterized by the lattice c/a ratio. @see SlipGeomHCPaBRYcaY1 in slipgeom/ECMech_slipgeom_hcp.h */
    using SlipGeom_HCP_A = SlipGeomHCPaBRYcaY1;
    /** @brief Kocks-Mecking balanced thermally-activated kinetics with the athermal/thermal-activation split enabled, the MTS `p`/`q` exponents pegged to 1 (simplified form), and per-slip-system hardening parameters (one dislocation-density variable per slip system, matching #SlipGeom_HCP_A's slip-system count) — needed since HCP slip families (basal/prismatic/pyramidal) have very different critical resolved shear stresses. @see KineticsKMBalD in kinetics/ECMech_kinetics_KMBalD.h */
    using Kin_HCP_A = KineticsKMBalD<true, true, true, true, SlipGeom_HCP_A::nslip>;
    /** @brief `"evptn_HCP_A"`: HCP slip (basal/prismatic/pyramidal-<a>/pyramidal-<c+a>) with per-family Kocks-Mecking kinetics (Kin_HCP_A) and #EVPTN_hex thermoelasticity. */
    using matModelEvptn_HCP_A = evptn::matModel<SlipGeom_HCP_A, Kin_HCP_A, EVPTN_hex, EOS_const_model>;

}
