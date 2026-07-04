/**
 * @file ECMech_cases_bcc_oro.cxx
 * @brief Factory for the 12-slip-system Orowan-kinetics BCC models (`evptn_BCC_C`,
 * `evptn_BCC_D`, `evptn_BCC_E`); dispatched to from `makeMatModelBCC` in
 * ECMech_cases_bcc.cxx.
 */

#include "ECMech_cases_bcc_defs.h"

namespace ecmech {

/**
 * @brief Build one of the 12-slip-system Orowan BCC models (C, D, E) from its exact
 * name.
 * @param modelName One of `"evptn_BCC_C"` (isotropic forest interaction),
 * `"evptn_BCC_D"` (anisotropic forest interaction), or `"evptn_BCC_E"` (non-Schmid slip
 * geometry with anisotropic forest interaction); see ECMech_cases_bcc_defs.h.
 * @return Newly heap-allocated `matModelBase*`, or `nullptr` if `modelName` is not one of
 * the above.
 */
__ecmech_host__
matModelBase* makeMatModelBCCOro(const std::string &modelName) {
    matModelBase* matModel = nullptr;

    if (modelName == "evptn_BCC_C") {
        auto mmECMEvptn = new ecmech::matModelEvptn_BCC_C();
        matModel = dynamic_cast<ecmech::matModelBase*>(mmECMEvptn);
    }
    else if (modelName == "evptn_BCC_D") {
        auto mmECMEvptn = new ecmech::matModelEvptn_BCC_D();
        matModel = dynamic_cast<ecmech::matModelBase*>(mmECMEvptn);
    }
    else if (modelName == "evptn_BCC_E") {
        auto mmECMEvptn = new ecmech::matModelEvptn_BCC_E();
        matModel = dynamic_cast<ecmech::matModelBase*>(mmECMEvptn);
    }

    return matModel;
}
}