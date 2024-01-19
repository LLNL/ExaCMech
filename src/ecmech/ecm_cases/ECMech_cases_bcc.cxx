#include "ECMech_cases.h"
#include "ECMech_cases_bcc_defs.h"

namespace ecmech {

/**
* @brief These are not the only possible cases -- they are here as a convenience
*/
__ecmech_host__
matModelBase* makeMatModelBCC(const std::string &modelName) {
    matModelBase* matModel = nullptr;

    if (modelName == "evptn_BCC_A") {
        auto mmECMEvptn = new ecmech::matModelEvptn_BCC_A();
        matModel = dynamic_cast<ecmech::matModelBase*>(mmECMEvptn);
    }
    else if (modelName == "evptn_BCC_AH") {
        auto mmECMEvptn = new ecmech::matModelEvptn_BCC_AH();
        matModel = dynamic_cast<ecmech::matModelBase*>(mmECMEvptn);
    }
    else if (modelName == "evptn_BCC_B") {
        auto mmECMEvptn = new ecmech::matModelEvptn_BCC_B();
        matModel = dynamic_cast<ecmech::matModelBase*>(mmECMEvptn);
    }
    else if (modelName == "evptn_BCC_C") {
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
    else if (modelName == "evptn_BCC_C_24") {
        auto mmECMEvptn = new ecmech::matModelEvptn_BCC_C_24();
        matModel = dynamic_cast<ecmech::matModelBase*>(mmECMEvptn);
    }
    else if (modelName == "evptn_BCC_D_24") {
        auto mmECMEvptn = new ecmech::matModelEvptn_BCC_D_24();
        matModel = dynamic_cast<ecmech::matModelBase*>(mmECMEvptn);
    }
    else if (modelName == "evptn_BCC_MD") {
        auto mmECMEvptn = new ecmech::matModelEvptn_BCC_MD();
        matModel = dynamic_cast<ecmech::matModelBase*>(mmECMEvptn);
    }

    return matModel;
}
}