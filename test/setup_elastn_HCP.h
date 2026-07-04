/**
 * @file setup_elastn_HCP.h
 *
 * @brief Fragment (see `setup_base.h` for the general `setup_*.h` inclusion pattern)
 * that sets hexagonal elastic constants (`ThermoElastNHexag`, see `ECMech_elastic.h`),
 * for either an already-declared `elastN` object or the flat `params` vector, depending
 * on `STACK_PARAMS`.
 */
{
   // Hexagonal elastic constants [C11, C12, C13, C33, C44, g_vecd2] -- see
   // ThermoElastNHexag::setParams in ECMech_elastic.h for the parameter order; g_vecd2
   // is the anisotropic Gruneisen (thermal) parameter, zeroed out here.
   double c11 = 162.4e-2, c12 = 92e-2, c13 = 69e-2, c33 = 180.7e-2, c44 = 46.7e-2;
   double g_vecd2 = 0.0;
   std::vector<double> paramsThese { c11, c12, c13, c33, c44, g_vecd2 };
#ifdef STACK_PARAMS
   params.insert(params.end(), paramsThese.begin(), paramsThese.end());
#else
   elastN.setParams(paramsThese);
#endif
}
