#include "ECMech_cases_bcc_defs.h"

namespace ecmech {

/**
* @brief These are not the only possible cases -- they are here as a convenience
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