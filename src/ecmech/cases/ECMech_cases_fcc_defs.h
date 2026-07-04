/**
 * @file ECMech_cases_fcc_defs.h
 * @brief Glossary of the concrete FCC (face-centered cubic) material models available
 * through the string names accepted by `makeMatModelFCC`/`makeMatModel`.
 *
 * Each `matModelEvptn_FCC_*` alias below is a distinct kinetics-model choice on top of
 * the fixed 12-system FCC slip geometry (SlipGeomFCC); the corresponding `evptn_FCC_*`
 * string name (used in `cases/ECMech_cases_fcc.cxx`) is what a calling code passes to
 * `ecmech::makeMatModel`/`ecmech::makeMatModelFCC` to obtain one.
 *
 * All FCC models here use #EVPTN_cubic thermoelasticity and #EOS_const_model.
 */

#pragma once

#include "ECMech_cases.h"

// Provide a set of defines that users can bring in if they'd like / have a need for them
namespace ecmech {
    // Make all of our Orowan models use the logrithmic formulation for better stability.
    /** @brief Orowan dislocation-density kinetics on the FCC slip geometry (SlipGeomFCC), using the logarithmic hardening-update form for numerical stability. @see KineticsOrowanD in kinetics/ECMech_kinetics_OrowanD.h */
    using Kin_OroD_Iso_FCC = KineticsOrowanD<false, false, false, true, false, 1, SlipGeomFCC, true>;

    /** @brief `"evptn_FCC_A"`: FCC slip with linear Voce (power-law + saturation) hardening. */
    using matModelEvptn_FCC_A = evptn::matModel<SlipGeomFCC, Kin_Voce, EVPTN_cubic, EOS_const_model >;
    /** @brief `"evptn_FCC_AH"`: same as matModelEvptn_FCC_A but with the nonlinear Voce hardening variant. */
    using matModelEvptn_FCC_AH = evptn::matModel<SlipGeomFCC, Kin_VoceNL, EVPTN_cubic, EOS_const_model >;
    /** @brief `"evptn_FCC_B"`: FCC slip with Kocks-Mecking balanced thermally-activated (MTS-like) kinetics, athermal/thermal-activation split disabled, and a single dislocation-density hardening variable. */
    using matModelEvptn_FCC_B = evptn::matModel<SlipGeomFCC, Kin_KMBalD_FFF, EVPTN_cubic, EOS_const_model >;
    /** @brief `"evptn_FCC_C"`: FCC slip with isotropic Orowan dislocation-density kinetics (Kin_OroD_Iso_FCC). */
    using matModelEvptn_FCC_C = evptn::matModel<SlipGeomFCC, Kin_OroD_Iso_FCC, EVPTN_cubic, EOS_const_model >;

}
