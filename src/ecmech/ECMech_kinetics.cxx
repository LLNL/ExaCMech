#include "ECMech_kinetics.h"

namespace ecmech {
   typedef KineticsKMBalD<true, false, false> KineticsKMBalD_PeierlsFF;
   typedef KineticsKMBalD<false, true, true> KineticsKMBalD_TT;
} // namespace ecmech
