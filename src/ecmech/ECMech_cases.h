// -*-c++-*-

#include "evptn/ECMech_evptn.h"
#include "kinetics/ECMech_kinetics.h"
#include "ECMech_slipgeom.h"
#include "ECMech_evptnWrap.h"

namespace ecmech {

   // constant non-isothermal EOS model
   using EOS_const_model = EosModelConst<false>;

   // some common kinetic forms
   using Kin_Voce   = KineticsVocePL<false>;
   using Kin_VoceNL = KineticsVocePL<true>;
   using Kin_KMBalD_TFF = KineticsKMBalD<true, false, false, false, 1>;
   using Kin_KMBalD_FFF = KineticsKMBalD<false, false, false, false, 1>;

   // shortening this as the type is a bit long
   using EVPTN_cubic = evptn::ThermoElastNCubic;
   using EVPTN_hex   = evptn::ThermoElastNHexag;

   // keeping only the most basic stuff here as all the typedefs have been moved elsewhere...

   __ecmech_host__
   matModelBase* makeMatModel(const std::string &modelName);
   __ecmech_host__
   matModelBase* makeMatModelFCC(const std::string &modelName);
   __ecmech_host__
   matModelBase* makeMatModelBCC(const std::string &modelName);
   __ecmech_host__
   matModelBase* makeMatModelHCP(const std::string &modelName);   

   __ecmech_host__
   std::tuple<size_t, size_t>
   modelNumParamsHist(const std::string_view &modelName);
   __ecmech_host__
   std::tuple<size_t, size_t>
   modelNumParamsHistFCC(const std::string_view &modelName);
   __ecmech_host__
   std::tuple<size_t, size_t>
   modelNumParamsHistBCC(const std::string_view &modelName);
   __ecmech_host__
   std::tuple<size_t, size_t>
   modelNumParamsHistHCP(const std::string_view &modelName);

}
