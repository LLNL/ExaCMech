#include "ECMech_evptnWrap.h"
#include "ECMech_eosSimple.h"

namespace ecmech {

   typedef matModelEvptn< SlipGeomFCC, KineticsVocePL, ThermoElastNCubic, EosModelConst<false> > matModelEvptn_FCC_A ;

} // namespace ecmech
