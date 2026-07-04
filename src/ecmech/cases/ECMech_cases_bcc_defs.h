/**
 * @file ECMech_cases_bcc_defs.h
 * @brief Glossary of the concrete BCC (body-centered cubic) material models available
 * through the string names accepted by `makeMatModelBCC`/`makeMatModel`.
 *
 * Each `matModelEvptn_BCC_*` alias below is a distinct combination of slip geometry and
 * kinetics model; the corresponding `evptn_BCC_*` string name (used in
 * `cases/ECMech_cases_bcc.cxx`, `cases/ECMech_cases_bcc_oro.cxx`, and
 * `cases/ECMech_cases_bcc_oro_big.cxx`) is what a calling code passes to
 * `ecmech::makeMatModel`/`ecmech::makeMatModelBCC` to obtain one.
 *
 * All BCC models here use #EVPTN_cubic thermoelasticity and #EOS_const_model.
 */

#pragma once

#include "ECMech_cases.h"

// Provide a set of defines that users can bring in if they'd like / have a need for them
namespace ecmech {

    /** @brief Standard 12 {110}<111> BCC slip systems. @see SlipGeomBCC in slipgeom/ECMech_slipgeom_bcc.h */
    using SlipGeom_BCC_A = SlipGeomBCC<12>;
    /** @brief Extended 24-system BCC slip geometry: the 12 {110}<111> systems plus 12 {112}<111> systems. @see SlipGeomBCC in slipgeom/ECMech_slipgeom_bcc.h */
    using SlipGeom_BCC_B = SlipGeomBCC<24>;

    // Make all of our Orowan models use the logrithmic formulation for better stability.
    /** @brief Orowan dislocation-density kinetics on #SlipGeom_BCC_A (12 slip systems) with an isotropic forest-interaction matrix, using the logarithmic hardening-update form for numerical stability. @see KineticsOrowanD in kinetics/ECMech_kinetics_OrowanD.h */
    using Kin_OroD_Iso_BCC = KineticsOrowanD<true, false, false, true, false, 1, SlipGeom_BCC_A, true>;
    /** @brief Same as #Kin_OroD_Iso_BCC but with a full anisotropic forest-interaction matrix instead of an isotropic one. */
    using Kin_OroD_Aniso_BCC = KineticsOrowanD<true, false, false, false, false, 1, SlipGeom_BCC_A, true>;
    /** @brief Same as #Kin_OroD_Iso_BCC but on the 24-system #SlipGeom_BCC_B geometry. */
    using Kin_OroD_Iso_BCC_24 = KineticsOrowanD<true, false, false, true, false, 1, SlipGeom_BCC_B, true>;
    /** @brief Same as #Kin_OroD_Aniso_BCC but on the 24-system #SlipGeom_BCC_B geometry. */
    using Kin_OroD_Aniso_BCC_24 = KineticsOrowanD<true, false, false, false, false, 1, SlipGeom_BCC_B, true>;
    /** @brief Anisotropic Orowan kinetics paired with the non-Schmid BCC slip geometry (SlipGeomBCCNonSchmid), for capturing non-Schmid (twinning/anti-twinning asymmetry) effects on slip activation. */
    using Kin_OroD_Aniso_BCC_NS = KineticsOrowanD<true, false, false, false, false, 1, SlipGeomBCCNonSchmid, true>;
    /** @brief Mobile-dislocation-density (MD) pencil-glide kinetics, paired with the stress-dependent pencil-glide slip geometry (SlipGeomBCCPencil). @see KineticsBCCMD in kinetics/ECMech_kinetics_BCCMD.h */
    using Kin_BCC_MD = KineticsBCCMD<SlipGeomBCCPencil>;

    /** @brief `"evptn_BCC_A"`: 12-system BCC slip with linear Voce (power-law + saturation) hardening. */
    using matModelEvptn_BCC_A = evptn::matModel<SlipGeom_BCC_A, Kin_Voce, EVPTN_cubic, EOS_const_model>;
    /** @brief `"evptn_BCC_AH"`: same as matModelEvptn_BCC_A but with the nonlinear Voce hardening variant. */
    using matModelEvptn_BCC_AH = evptn::matModel<SlipGeom_BCC_A, Kin_VoceNL, EVPTN_cubic, EOS_const_model>;
    /** @brief `"evptn_BCC_B"`: 12-system BCC slip with Kocks-Mecking balanced thermally-activated (MTS-like) kinetics and a single dislocation-density hardening variable. */
    using matModelEvptn_BCC_B = evptn::matModel<SlipGeom_BCC_A, Kin_KMBalD_TFF, EVPTN_cubic, EOS_const_model>;
    /** @brief `"evptn_BCC_C"`: 12-system BCC slip with isotropic Orowan dislocation-density kinetics (Kin_OroD_Iso_BCC). */
    using matModelEvptn_BCC_C = evptn::matModel<SlipGeom_BCC_A, Kin_OroD_Iso_BCC, EVPTN_cubic, EOS_const_model>;
    /** @brief `"evptn_BCC_C_24"`: 24-system BCC slip with isotropic Orowan dislocation-density kinetics (Kin_OroD_Iso_BCC_24). */
    using matModelEvptn_BCC_C_24 = evptn::matModel<SlipGeom_BCC_B, Kin_OroD_Iso_BCC_24, EVPTN_cubic, EOS_const_model>;
    /** @brief `"evptn_BCC_D"`: 12-system BCC slip with anisotropic Orowan dislocation-density kinetics (Kin_OroD_Aniso_BCC). */
    using matModelEvptn_BCC_D = evptn::matModel<SlipGeom_BCC_A, Kin_OroD_Aniso_BCC, EVPTN_cubic, EOS_const_model>;
    /** @brief `"evptn_BCC_D_24"`: 24-system BCC slip with anisotropic Orowan dislocation-density kinetics (Kin_OroD_Aniso_BCC_24). */
    using matModelEvptn_BCC_D_24 = evptn::matModel<SlipGeom_BCC_B, Kin_OroD_Aniso_BCC_24, EVPTN_cubic, EOS_const_model>;
    /** @brief `"evptn_BCC_E"`: non-Schmid BCC slip geometry with anisotropic Orowan dislocation-density kinetics (Kin_OroD_Aniso_BCC_NS). */
    using matModelEvptn_BCC_E = evptn::matModel<SlipGeomBCCNonSchmid, Kin_OroD_Aniso_BCC_NS, EVPTN_cubic, EOS_const_model>;
    /** @brief `"evptn_BCC_MD"`: stress-dependent pencil-glide BCC slip geometry with mobile-dislocation-density (MD) kinetics (Kin_BCC_MD). */
    using matModelEvptn_BCC_MD = evptn::matModel<SlipGeomBCCPencil, Kin_BCC_MD, EVPTN_cubic, EOS_const_model>;

    /** @brief Build one of the non-Orowan BCC models (A, AH, B, MD) from its string name; see `cases/ECMech_cases_bcc.cxx`. */
    __ecmech_host__
    matModelBase* makeMatModelBCCNorm(const std::string &modelName);
    /** @brief Build one of the 12-slip-system Orowan BCC models (C, D, E) from its string name; see `cases/ECMech_cases_bcc_oro.cxx`. */
    __ecmech_host__
    matModelBase* makeMatModelBCCOro(const std::string &modelName);
    /** @brief Build one of the 24-slip-system Orowan BCC models (C_24, D_24) from its string name; see `cases/ECMech_cases_bcc_oro_big.cxx`. */
    __ecmech_host__
    matModelBase* makeMatModelBCCOroBig(const std::string &modelName);

}

