#include "ECMech_cases.h"
#include "ECMech_cases_bcc_defs.h"

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

/**
* @brief These are not the only possible cases -- they are here as a convenience
*/
__ecmech_host__
std::tuple<size_t, size_t>
modelNumParamsHistBCC(const std::string_view &modelName) {

    if (modelName == "evptn_BCC_A") {
        return std::tuple(ecmech::matModelEvptn_BCC_A::nParams, ecmech::matModelEvptn_BCC_A::numHist);
    }
    else if (modelName == "evptn_BCC_AH") {
        return std::tuple(ecmech::matModelEvptn_BCC_AH::nParams, ecmech::matModelEvptn_BCC_AH::numHist);
    }
    else if (modelName == "evptn_BCC_B") {
        return std::tuple(ecmech::matModelEvptn_BCC_B::nParams, ecmech::matModelEvptn_BCC_B::numHist);
    } else if (modelName == "evptn_BCC_C") {
        return std::tuple(ecmech::matModelEvptn_BCC_C::nParams, ecmech::matModelEvptn_BCC_C::numHist);
    }
    else if (modelName == "evptn_BCC_D") {
        return std::tuple(ecmech::matModelEvptn_BCC_D::nParams, ecmech::matModelEvptn_BCC_D::numHist);
    } else if (modelName == "evptn_BCC_E") {
        return std::tuple(ecmech::matModelEvptn_BCC_E::nParams, ecmech::matModelEvptn_BCC_E::numHist);
    }
    else if (modelName == "evptn_BCC_C_24") {
        return std::tuple(ecmech::matModelEvptn_BCC_C_24::nParams, ecmech::matModelEvptn_BCC_C_24::numHist);
    }
    else if (modelName == "evptn_BCC_D_24") {
        return std::tuple(ecmech::matModelEvptn_BCC_D_24::nParams, ecmech::matModelEvptn_BCC_D_24::numHist);
    }
    else if (modelName == "evptn_BCC_MD") {
        return std::tuple(ecmech::matModelEvptn_BCC_MD::nParams, ecmech::matModelEvptn_BCC_MD::numHist);
    }
    return std::tuple(0, 0);
}

}