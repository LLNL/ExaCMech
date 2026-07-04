/**
 * @file setup_eos.h
 *
 * @brief Fragment (see `setup_base.h` for the general `setup_*.h` inclusion pattern)
 * that sets the simple EOS model (`EosModelConst`, see `ECMech_eosSimple.h`), for
 * either an already-declared `eos` object or the flat `params` vector, depending on
 * `STACK_PARAMS`.
 *
 * The two modes genuinely differ here, not just in *where* the values go: in the direct
 * (non-`STACK_PARAMS`) case, `EosModelConst::setParams` takes all 5 parameters `[ρ₀,
 * K₀, cᵥ, Γ, ε₀]` (with the bulk modulus read from the already-configured `elastN`
 * object via `getBulkMod()`, since a real caller wouldn't otherwise know it yet). But
 * `matModel::initFromParams` (`ECMech_evptnWrap.h`) computes `ρ₀`/`K₀`/`cᵥ` itself
 * (`ρ₀`/`cᵥ` from `density0`/`cvav` above, `K₀` from the elastic model it just built) --
 * so in `STACK_PARAMS` mode this fragment only pushes the *remaining* `nParamsEOS =
 * EosModel::nParams - mmodel->nParamsEOSHave` parameters (`Γ`, `ε₀`) onto `params`, and
 * uses a `-1.0` placeholder for the bulk modulus since it's never actually read in that
 * mode.
 */
{
#ifdef STACK_PARAMS
   double bulk_modulus = -1.0; // dummy
#else
   double bulk_modulus = elastN.getBulkMod();
#endif
   double gamma = 1.7;
   double cold_energy0 = -cvav * 300.;
   const std::vector<double> paramsThese { density0, bulk_modulus, cvav, gamma, cold_energy0 };
#ifdef STACK_PARAMS
   int nParamsEOS = paramsThese.size() - mmodel->nParamsEOSHave; // nasty complexity to match what happens in matModel
   for (int iP = 0; iP<nParamsEOS; ++iP) {
      params.push_back(paramsThese[mmodel->nParamsEOSHave + iP]);
   }

#else
   eos.setParams(paramsThese);
#endif
}
