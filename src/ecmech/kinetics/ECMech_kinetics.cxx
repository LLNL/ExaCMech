/**
 * @file ECMech_kinetics.cxx
 * @brief A couple of additional named `KineticsKMBalD` convenience typedefs, not
 * currently wired up to any `cases/` model but available for direct use.
 */

#include "ECMech_kinetics.h"

namespace ecmech {
   /**
    * @brief Kocks-Mecking balanced thermally-activated kinetics with the
    * athermal/thermal-activation split enabled (the athermal part associated with the
    * Peierls barrier) and the MTS `p`/`q` exponents left general (not pegged to 1).
    * @see KineticsKMBalD in ECMech_kinetics_KMBalD.h
    */
   typedef KineticsKMBalD<true, false, false, false, 1> KineticsKMBalD_PeierlsFF;
   /**
    * @brief Kocks-Mecking balanced thermally-activated kinetics with the
    * athermal/thermal-activation split disabled and the MTS `p`/`q` exponents pegged to
    * 1 (simplified form).
    * @see KineticsKMBalD in ECMech_kinetics_KMBalD.h
    */
   typedef KineticsKMBalD<false, true, true, false, 1> KineticsKMBalD_TT;
} // namespace ecmech
