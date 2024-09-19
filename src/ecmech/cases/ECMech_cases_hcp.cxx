#include "ECMech_cases.h"
#include "ECMech_cases_hcp_defs.h"

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
std::tuple<size_t, size_t>
modelNumParamsHistHCP(const std::string_view &modelName) {
    if (modelName == "evptn_HCP_A") {
        return std::tuple(ecmech::matModelEvptn_HCP_A::nParams, ecmech::matModelEvptn_HCP_A::numHist);
    }
    return std::tuple(0, 0);
}

}