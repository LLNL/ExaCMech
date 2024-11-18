#include "ECMech_cases_bcc_defs.h"
#include "ECMech_cases_util.h"

namespace ecmech {

/**
* @brief These are not the only possible cases -- they are here as a convenience
*/
__ecmech_host__
matModelBase* makeMatModelBCC(const std::string &modelName) {
    matModelBase* matModel = nullptr;

    std::array<std::string, 4> norm = {"evptn_BCC_A", "evptn_BCC_AH", "evptn_BCC_B", "evptn_BCC_MD"};
    std::array<std::string, 3> oro = {"evptn_BCC_C", "evptn_BCC_D", "evptn_BCC_E"};
    std::array<std::string, 3> oro_big = {"evptn_BCC_C_24", "evptn_BCC_D_24"};

    auto find_case = [=] (auto& string_comp, auto& string_array) -> bool {
        auto it = std::find_if(string_array.begin(), string_array.end(),
                        [&](const auto st)
                        { return st.find(string_comp) != std::string::npos; });
        return (it != string_array.end());
    };

    if (find_case(modelName, norm)) {
        matModel = makeMatModelBCCNorm(modelName);
    }
    else if (find_case(modelName, oro)) {
        matModel = makeMatModelBCCOro(modelName);
    }
    else if (find_case(modelName, oro_big)) {
        matModel = makeMatModelBCCOroBig(modelName);
    }

    return matModel;
}

__ecmech_host__
matModelBase* makeMatModelBCCNorm(const std::string &modelName) {
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
    else if (modelName == "evptn_BCC_MD") {
        auto mmECMEvptn = new ecmech::matModelEvptn_BCC_MD();
        matModel = dynamic_cast<ecmech::matModelBase*>(mmECMEvptn);
    }

    return matModel;
}

/**
* @brief These are not the only possible cases -- they are here as a convenience
*/
__ecmech_host__
std::map<std::string, size_t>
modelParamIndexMapBCC(const std::string_view &modelName) {

    if (modelName == "evptn_BCC_A") {
        return NumParamIndexInfo<ecmech::matModelEvptn_BCC_A>().m_maps;
    }
    else if (modelName == "evptn_BCC_AH") {
        return NumParamIndexInfo<ecmech::matModelEvptn_BCC_AH>().m_maps;
    }
    else if (modelName == "evptn_BCC_B") {
        return NumParamIndexInfo<ecmech::matModelEvptn_BCC_B>().m_maps;
    } else if (modelName == "evptn_BCC_C") {
        return NumParamIndexInfo<ecmech::matModelEvptn_BCC_C>().m_maps;
    }
    else if (modelName == "evptn_BCC_D") {
        return NumParamIndexInfo<ecmech::matModelEvptn_BCC_D>().m_maps;
    } else if (modelName == "evptn_BCC_E") {
        return NumParamIndexInfo<ecmech::matModelEvptn_BCC_E>().m_maps;
    }
    else if (modelName == "evptn_BCC_C_24") {
       return NumParamIndexInfo<ecmech::matModelEvptn_BCC_C_24>().m_maps;
    }
    else if (modelName == "evptn_BCC_D_24") {
        return NumParamIndexInfo<ecmech::matModelEvptn_BCC_D_24>().m_maps;
    }
    else if (modelName == "evptn_BCC_MD") {
        return NumParamIndexInfo<ecmech::matModelEvptn_BCC_MD>().m_maps;
    }
    return {};
}

}