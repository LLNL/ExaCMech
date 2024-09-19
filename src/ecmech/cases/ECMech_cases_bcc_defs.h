#pragma once

#include "ECMech_cases.h"

// Provide a set of defines that users can bring in if they'd like / have a need for them
namespace ecmech {

    using SlipGeom_BCC_A = SlipGeomBCC<12>;
    using SlipGeom_BCC_B = SlipGeomBCC<24>;

    // Make all of our Orowan models use the logrithmic formulation for better stability.
    using Kin_OroD_Iso_BCC = KineticsOrowanD<true, false, false, true, false, 1, SlipGeom_BCC_A, true>;
    using Kin_OroD_Aniso_BCC = KineticsOrowanD<true, false, false, false, false, 1, SlipGeom_BCC_A, true>;
    using Kin_OroD_Iso_BCC_24 = KineticsOrowanD<true, false, false, true, false, 1, SlipGeom_BCC_B, true>;
    using Kin_OroD_Aniso_BCC_24 = KineticsOrowanD<true, false, false, false, false, 1, SlipGeom_BCC_B, true>;
    using Kin_OroD_Aniso_BCC_NS = KineticsOrowanD<true, false, false, false, false, 1, SlipGeomBCCNonSchmid, true>;
    using Kin_BCC_MD = KineticsBCCMD<SlipGeomBCCPencil>;

    using matModelEvptn_BCC_A = evptn::matModel<SlipGeom_BCC_A, Kin_Voce, EVPTN_cubic, EOS_const_model>;
    using matModelEvptn_BCC_AH = evptn::matModel<SlipGeom_BCC_A, Kin_VoceNL, EVPTN_cubic, EOS_const_model>;
    using matModelEvptn_BCC_B = evptn::matModel<SlipGeom_BCC_A, Kin_KMBalD_TFF, EVPTN_cubic, EOS_const_model>;
    using matModelEvptn_BCC_C = evptn::matModel<SlipGeom_BCC_A, Kin_OroD_Iso_BCC, EVPTN_cubic, EOS_const_model>;
    using matModelEvptn_BCC_C_24 = evptn::matModel<SlipGeom_BCC_B, Kin_OroD_Iso_BCC_24, EVPTN_cubic, EOS_const_model>;
    using matModelEvptn_BCC_D = evptn::matModel<SlipGeom_BCC_A, Kin_OroD_Aniso_BCC, EVPTN_cubic, EOS_const_model>;
    using matModelEvptn_BCC_D_24 = evptn::matModel<SlipGeom_BCC_B, Kin_OroD_Aniso_BCC_24, EVPTN_cubic, EOS_const_model>;
    using matModelEvptn_BCC_E = evptn::matModel<SlipGeomBCCNonSchmid, Kin_OroD_Aniso_BCC_NS, EVPTN_cubic, EOS_const_model>;
    using matModelEvptn_BCC_MD = evptn::matModel<SlipGeomBCCPencil, Kin_BCC_MD, EVPTN_cubic, EOS_const_model>;

    __ecmech_host__
    matModelBase* makeMatModelBCCNorm(const std::string &modelName);
    __ecmech_host__
    matModelBase* makeMatModelBCCOro(const std::string &modelName);
    __ecmech_host__
    matModelBase* makeMatModelBCCOroBig(const std::string &modelName);

}

