/**
 * @file ECMech_cases_hcp.cxx
 * @brief HCP model-name factory and parameter/history-index lookup, dispatched to from
 * `makeMatModel`/`modelParamIndexMap` in ECMech_cases.cxx.
 */

#include "ECMech_cases_hcp_defs.h"
#include "ECMech_cases_util.h"


namespace ecmech {

/**
 * @brief Build an HCP material model from its exact name.
 * @param modelName Currently only `"evptn_HCP_A"` is recognized; see
 * ECMech_cases_hcp_defs.h.
 * @return Newly heap-allocated `matModelBase*`, or `nullptr` if `modelName` is not
 * recognized.
 */
__ecmech_host__
matModelBase* makeMatModelHCP(const std::string &modelName) {
    matModelBase* matModel = nullptr;
    if (modelName == "evptn_HCP_A") {
        auto mmECMEvptn = new ecmech::matModelEvptn_HCP_A();
        matModel = dynamic_cast<ecmech::matModelBase*>(mmECMEvptn);
    }

    return matModel;
}

/**
 * @brief Look up parameter-count and history-array-index information for an HCP model by
 * its exact name.
 * @param modelName Currently only `"evptn_HCP_A"` is recognized.
 * @return Lookup map as described in ECMech_cases.h's modelParamIndexMap, or an empty
 * map if `modelName` is not recognized.
 */
__ecmech_host__
std::map<std::string, size_t>
modelParamIndexMapHCP(const std::string_view &modelName) {
    if (modelName == "evptn_HCP_A") {
        return NumParamIndexInfo<ecmech::matModelEvptn_HCP_A>().m_maps;
    }
    return {};
}

}