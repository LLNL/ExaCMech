/**
 * @file setup_elastn.h
 *
 * @brief Fragment (see `setup_base.h` for the general `setup_*.h` inclusion pattern)
 * that sets cubic elastic constants (`ThermoElastNCubic`, see `ECMech_elastic.h`), for
 * either an already-declared `elastN` object or the flat `params` vector, depending on
 * `STACK_PARAMS`.
 */
{
   // Cubic elastic constants [C11, C12, C44] -- see ThermoElastNCubic::setParams in
   // ECMech_elastic.h for the parameter order and stability constraints these satisfy.
   double c11 = 300e-2, c12 = 100e-2, c44 = 100e-2;
   std::vector<double> paramsThese { c11, c12, c44 };
#ifdef STACK_PARAMS
   params.insert(params.end(), paramsThese.begin(), paramsThese.end());
#else
   elastN.setParams(paramsThese);
#endif
}
