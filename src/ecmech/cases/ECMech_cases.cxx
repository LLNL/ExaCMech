/**
 * @file ECMech_cases.cxx
 * @brief Top-level model-name dispatch: routes a model name to the FCC/BCC/HCP family
 * based on which of those substrings it contains, then delegates to the matching
 * family's `makeMatModel*`/`modelParamIndexMap*` implementation.
 */

#include "ECMech_cases.h"

namespace ecmech {
/**
 * @brief Build a material model from its name.
 *
 * Dispatches on whether `modelName` contains "FCC", "BCC", or "HCP" and forwards to
 * makeMatModelFCC/makeMatModelBCC/makeMatModelHCP accordingly. This list of families is
 * not exhaustive by design -- new crystal-structure families can be added as additional
 * `else if` branches without needing to change this function's signature.
 * @param modelName Model name to build, e.g. `"evptn_FCC_A"`.
 * @return Newly heap-allocated `matModelBase*`.
 * @throws (via ECMECH_FAIL) if `modelName` contains none of "FCC"/"BCC"/"HCP", or if the
 * matching family does not recognize the specific name (returns `nullptr` internally,
 * which is caught here and turned into a failure).
 */
__ecmech_host__
matModelBase* makeMatModel(const std::string &modelName) {
   matModelBase* matModel = nullptr;

   if (modelName.find("FCC") != std::string::npos) {
      matModel = makeMatModelFCC(modelName);
   }
   else if (modelName.find("BCC") != std::string::npos) {
      matModel = makeMatModelBCC(modelName);
   }
   else if (modelName.find("HCP") != std::string::npos) {
      matModel = makeMatModelHCP(modelName);
   }

   if (matModel == nullptr) {
      std::string msg = std::string("model name not recognized : ") + modelName;
      ECMECH_FAIL(__func__, msg.c_str());
   }

   return matModel;
}

/**
 * @brief Look up parameter-count and history-array-index information for a model by
 * name.
 *
 * Uses the same "FCC"/"BCC"/"HCP" substring dispatch as makeMatModel, forwarding to
 * modelParamIndexMapFCC/modelParamIndexMapBCC/modelParamIndexMapHCP.
 * @param modelName Model name to look up, e.g. `"evptn_FCC_A"`.
 * @return Lookup map as described in ECMech_cases.h's modelParamIndexMap.
 * @throws (via ECMECH_FAIL) if `modelName` contains none of "FCC"/"BCC"/"HCP". Unlike
 * makeMatModel, an unrecognized name *within* a matched family is not itself treated as
 * a failure here -- the family-level function simply returns an empty map.
 */
__ecmech_host__
std::map<std::string, size_t>
modelParamIndexMap(const std::string_view &modelName) {
   if (modelName.find("FCC") != std::string_view::npos) {
      return modelParamIndexMapFCC(modelName);
   }
   else if (modelName.find("BCC") != std::string_view::npos) {
      return modelParamIndexMapBCC(modelName);
   }
   else if (modelName.find("HCP") != std::string_view::npos) {
      return modelParamIndexMapHCP(modelName);
   }
   else {
      std::string msg = std::string("model name not recognized : ") + std::string(modelName);
      ECMECH_FAIL(__func__, msg.c_str());
   }
   return {};
}
} // namespace ecmech
