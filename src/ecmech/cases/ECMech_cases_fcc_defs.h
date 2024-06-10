#pragma once

#include "ECMech_cases.h"

// Provide a set of defines that users can bring in if they'd like / have a need for them
namespace ecmech {

    using Kin_OroD_Iso_FCC = KineticsOrowanD<false, false, false, true, false, 1, SlipGeomFCC>;

    using matModelEvptn_FCC_A = evptn::matModel<SlipGeomFCC, Kin_Voce, EVPTN_cubic, EOS_const_model >;
    using matModelEvptn_FCC_AH = evptn::matModel<SlipGeomFCC, Kin_VoceNL, EVPTN_cubic, EOS_const_model >;
    using matModelEvptn_FCC_B = evptn::matModel<SlipGeomFCC, Kin_KMBalD_FFF, EVPTN_cubic, EOS_const_model >;
    using matModelEvptn_FCC_C = evptn::matModel<SlipGeomFCC, Kin_OroD_Iso_FCC, EVPTN_cubic, EOS_const_model >;

}