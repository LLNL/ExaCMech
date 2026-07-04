/**
 * @file ECMech_cases_fcc.cxx
 * @brief FCC model-name factory and parameter/history-index lookup, dispatched to from
 * `makeMatModel`/`modelParamIndexMap` in ECMech_cases.cxx.
 */

#include "ECMech_cases_fcc_defs.h"
#include "ECMech_cases_util.h"


namespace ecmech {

/**
 * @brief Build an FCC material model from its exact name.
 * @param modelName One of `"evptn_FCC_A"` (linear Voce), `"evptn_FCC_AH"` (nonlinear
 * Voce), `"evptn_FCC_B"` (Kocks-Mecking), or `"evptn_FCC_C"` (Orowan dislocation
 * density); see ECMech_cases_fcc_defs.h for what each represents physically.
 * @return Newly heap-allocated `matModelBase*`, or `nullptr` if `modelName` is not one of
 * the above.
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
    else if (modelName == "evptn_FCC_C") {
        auto mmECMEvptn = new ecmech::matModelEvptn_FCC_C();
        matModel = dynamic_cast<ecmech::matModelBase*>(mmECMEvptn);
    }

    return matModel;
}

/**
 * @brief Look up parameter-count and history-array-index information for an FCC model by
 * its exact name.
 * @param modelName One of `"evptn_FCC_A"`, `"evptn_FCC_AH"`, `"evptn_FCC_B"`, or
 * `"evptn_FCC_C"`; see ECMech_cases_fcc_defs.h.
 * @return Lookup map as described in ECMech_cases.h's modelParamIndexMap, or an empty
 * map if `modelName` is not recognized.
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