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

/**
* @brief These are not the only possible cases -- they are here as a convenience
*/
__ecmech_host__
std::tuple<size_t, size_t>
modelNumParamsHistFCC(const std::string_view &modelName) {

    if (modelName == "evptn_FCC_A") {
        return std::tuple(ecmech::matModelEvptn_FCC_A::nParams, ecmech::matModelEvptn_FCC_A::numHist);
    }
    else if (modelName == "evptn_FCC_AH") {
        return std::tuple(ecmech::matModelEvptn_FCC_AH::nParams, ecmech::matModelEvptn_FCC_AH::numHist);
    }
    else if (modelName == "evptn_FCC_B") {
        return std::tuple(ecmech::matModelEvptn_FCC_B::nParams, ecmech::matModelEvptn_FCC_B::numHist);
    } else if (modelName == "evptn_FCC_C") {
        return std::tuple(ecmech::matModelEvptn_FCC_C::nParams, ecmech::matModelEvptn_FCC_C::numHist);
    }

    return std::tuple(0, 0);
}

}