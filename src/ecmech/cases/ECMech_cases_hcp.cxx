#include "ECMech_cases_hcp_defs.h"
#include "ECMech_cases_util.h"


namespace ecmech {

__ecmech_host__
matModelBase* makeMatModelHCP(const std::string &modelName) {
    matModelBase* matModel = nullptr;
    if (modelName == "evptn_HCP_A") {
        auto mmECMEvptn = new ecmech::matModelEvptn_HCP_A();
        matModel = dynamic_cast<ecmech::matModelBase*>(mmECMEvptn);
    }

    return matModel;
}

__ecmech_host__
std::map<std::string, size_t>
modelParamIndexMapHCP(const std::string_view &modelName) {
    if (modelName == "evptn_HCP_A") {
        return NumParamIndexInfo<ecmech::matModelEvptn_HCP_A>().m_maps;
    }
    return {};
}

}