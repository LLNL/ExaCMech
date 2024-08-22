#include "ECMech_kinetics.h"

namespace ecmech {
   typedef KineticsKMBalD<true, false, false, false, 1> KineticsKMBalD_PeierlsFF;
   typedef KineticsKMBalD<false, true, true, false, 1> KineticsKMBalD_TT;
} // namespace ecmech
