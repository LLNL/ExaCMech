#include "ECMech_cases.h"

namespace ecmech {
/**
 * @brief These are not the only possible cases -- they are here as a convenience
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
