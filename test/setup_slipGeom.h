/**
 * @file setup_slipGeom.h
 *
 * @brief Fragment (see `setup_base.h` for the general `setup_*.h` inclusion pattern)
 * that sets slip-geometry parameters for `SlipGeomFCC`/`SlipGeom_BCC_A` -- both have
 * `nParams == 0`, so this is a no-op that exists purely so every test case can
 * `#include` a slip-geometry setup fragment uniformly regardless of crystal structure.
 */
{
   std::vector<double> paramsThese; // no parameters for FCC or BCC
#ifdef STACK_PARAMS
   params.insert(params.end(), paramsThese.begin(), paramsThese.end());
#else
   slipGeom.setParams(paramsThese);
#endif
}
