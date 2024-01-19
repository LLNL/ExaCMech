#include "ECMech_cases.h"
#include "ECMech_cases_fcc_defs.h"

namespace ecmech {

/**
* @brief These are not the only possible cases -- they are here as a convenience
*/
__ecmech_host__
matModelBase* makeMatModelFCC(const std::string &modelName) {
    matModelBase* matModel = nullptr;

    if (modelName == "evptn_FCC_A") {
        auto mmECMEvptn = new ecmech::matModelEvptn_FCC_A();
        matModel = dynamic_cast<ecmech::matModelBase*>(mmECMEvptn);
    }
    else if (modelName == "evptn_FCC_AH") {
        auto mmECMEvptn = new ecmech::matModelEvptn_FCC_AH();
        matModel = dynamic_cast<ecmech::matModelBase*>(mmECMEvptn);
    }
    else if (modelName == "evptn_FCC_B") {
        auto mmECMEvptn = new ecmech::matModelEvptn_FCC_B();
        matModel = dynamic_cast<ecmech::matModelBase*>(mmECMEvptn);
    }

    return matModel;
}
}