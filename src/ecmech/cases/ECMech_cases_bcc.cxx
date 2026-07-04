/**
 * @file ECMech_cases_bcc.cxx
 * @brief BCC model-name dispatch: groups the `evptn_BCC_*` names into the "normal"
 * (non-Orowan), 12-slip-system Orowan, and 24-slip-system Orowan sub-families and
 * routes to the matching `makeMatModelBCC*` implementation (in this file,
 * ECMech_cases_bcc_oro.cxx, and ECMech_cases_bcc_oro_big.cxx respectively).
 */

#include "ECMech_cases_bcc_defs.h"
#include "ECMech_cases_util.h"

namespace ecmech {

/**
 * @brief Build a BCC material model from its name.
 *
 * Classifies `modelName` into one of three groups -- `norm` (A, AH, B, MD), `oro` (C, D,
 * E), or `oro_big` (C_24, D_24) -- by checking whether the canonical name for any model
 * in that group contains `modelName` as a substring, then delegates to
 * makeMatModelBCCNorm/makeMatModelBCCOro/makeMatModelBCCOroBig accordingly. See
 * ECMech_cases_bcc_defs.h for what each concrete model name means physically.
 * @param modelName Model name to build, e.g. `"evptn_BCC_A"`.
 * @return Newly heap-allocated `matModelBase*`, or `nullptr` if `modelName` does not
 * match any known BCC model.
 */
__ecmech_host__
matModelBase* makeMatModelBCC(const std::string &modelName) {
    matModelBase* matModel = nullptr;

    std::array<std::string, 4> norm = {"evptn_BCC_A", "evptn_BCC_AH", "evptn_BCC_B", "evptn_BCC_MD"};
    std::array<std::string, 3> oro = {"evptn_BCC_C", "evptn_BCC_D", "evptn_BCC_E"};
    std::array<std::string, 2> oro_big = {"evptn_BCC_C_24", "evptn_BCC_D_24"};

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

/**
 * @brief Build one of the non-Orowan BCC models (A, AH, B, MD) from its exact name.
 * @param modelName One of `"evptn_BCC_A"`, `"evptn_BCC_AH"`, `"evptn_BCC_B"`,
 * `"evptn_BCC_MD"`; see ECMech_cases_bcc_defs.h for what each represents physically.
 * @return Newly heap-allocated `matModelBase*`, or `nullptr` if `modelName` is not one of
 * the above.
 */
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
 * @brief Look up parameter-count and history-array-index information for a BCC model by
 * its exact name.
 * @param modelName One of the `evptn_BCC_*` names declared in ECMech_cases_bcc_defs.h.
 * @return Lookup map as described in ECMech_cases.h's modelParamIndexMap, or an empty
 * map if `modelName` is not recognized.
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