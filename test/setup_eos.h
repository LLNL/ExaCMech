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
