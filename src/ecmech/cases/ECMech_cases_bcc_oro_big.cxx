/**
 * @file ECMech_cases_bcc_oro_big.cxx
 * @brief Factory for the 24-slip-system ("big") Orowan-kinetics BCC models
 * (`evptn_BCC_C_24`, `evptn_BCC_D_24`); dispatched to from `makeMatModelBCC` in
 * ECMech_cases_bcc.cxx.
 */

#include "ECMech_cases_bcc_defs.h"

namespace ecmech {

/**
 * @brief Build one of the 24-slip-system Orowan BCC models (C_24, D_24) from its exact
 * name.
 * @param modelName One of `"evptn_BCC_C_24"` (isotropic forest interaction) or
 * `"evptn_BCC_D_24"` (anisotropic forest interaction); see ECMech_cases_bcc_defs.h.
 * @return Newly heap-allocated `matModelBase*`, or `nullptr` if `modelName` is not one of
 * the above.
 */
__ecmech_host__
matModelBase* makeMatModelBCCOroBig(const std::string &modelName) {
    matModelBase* matModel = nullptr;

    if (modelName == "evptn_BCC_C_24") {
        auto mmECMEvptn = new ecmech::matModelEvptn_BCC_C_24();
        matModel = dynamic_cast<ecmech::matModelBase*>(mmECMEvptn);
    }
    else if (modelName == "evptn_BCC_D_24") {
        auto mmECMEvptn = new ecmech::matModelEvptn_BCC_D_24();
        matModel = dynamic_cast<ecmech::matModelBase*>(mmECMEvptn);
    }

    return matModel;
}
}