/**
 * @file setup_slipGeom_HCP.h
 *
 * @brief Fragment (see `setup_base.h` for the general `setup_*.h` inclusion pattern)
 * that sets the single HCP slip-geometry parameter (`SlipGeom_HCP_A::nParams == 1`, the
 * c/a lattice ratio; see `ECMech_slipgeom_hcp.h`), for either an already-declared
 * `slipGeom` object or the flat `params` vector, depending on `STACK_PARAMS`.
 */
{
   double covera = 1.59773226818; // c/a ratio (dimensionless)
   std::vector<double> paramsThese { covera };
#ifdef STACK_PARAMS
   params.insert(params.end(), paramsThese.begin(), paramsThese.end());
#else
   slipGeom.setParams(paramsThese);
#endif
}
