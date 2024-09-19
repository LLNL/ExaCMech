#include "ECMech_cases_fcc_defs.h"
#include "ECMech_cases_util.h"


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
std::map<std::string, size_t>
modelParamIndexMapFCC(const std::string_view &modelName) {
    if (modelName == "evptn_FCC_A") {
        return NumParamIndexInfo<ecmech::matModelEvptn_FCC_A>().m_maps;
    }
    else if (modelName == "evptn_FCC_AH") {
        return NumParamIndexInfo<ecmech::matModelEvptn_FCC_AH>().m_maps;
    }
    else if (modelName == "evptn_FCC_B") {
        return NumParamIndexInfo<ecmech::matModelEvptn_FCC_B>().m_maps;
    } else if (modelName == "evptn_FCC_C") {
        return NumParamIndexInfo<ecmech::matModelEvptn_FCC_C>().m_maps;
    }

    return {};
}

}