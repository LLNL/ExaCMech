#pragma once

#include "ECMech_cases.h"

// We currently don't have a ton of HCP use cases but that could change in the future
// Provide a set of defines that users can bring in if they'd like / have a need for them
namespace ecmech {

    using SlipGeom_HCP_A = SlipGeomHCPaBRYcaY1;
    using Kin_HCP_A = KineticsKMBalD<true, true, true, true, SlipGeom_HCP_A::nslip>;
    using matModelEvptn_HCP_A = evptn::matModel<SlipGeom_HCP_A, Kin_HCP_A, EVPTN_hex, EOS_const_model>;

}