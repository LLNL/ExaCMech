#include "ECMech_evptnWrap.h"
#include "ECMech_eosSimple.h"
#include "ECMech_evptn.h"
#include "ECMech_kinetics.h"
#include "ECMech_slipgeom.h"

namespace ecmech {

   // TO_DO : not sure that this is really useful
   typedef evptn::matModel< SlipGeomFCC, KineticsVocePL, evptn::ThermoElastNCubic, EosModelConst<false> > matModelEvptn_FCC_A ;

} // namespace ecmech
