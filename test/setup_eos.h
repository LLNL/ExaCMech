{
#ifdef STACK_PARAMS
   double bulkMod = -1.0;     // dummy
#else
   double bulkMod = elastN.getBulkMod();
#endif
   double gamma = 1.7;
   double ec0 = -cvav * 300.;
   const std::vector<double> paramsThese{ rho0, bulkMod, cvav, gamma, ec0 };
#ifdef STACK_PARAMS
   int nParamsEOS = paramsThese.size() - mmodel->nParamsEOSHave;   // nasty complexity to match what happens in matModel
   for (int iP = 0; iP<nParamsEOS; ++iP) {
      params.push_back(paramsThese[mmodel->nParamsEOSHave + iP]);
   }

#else
   eos.setParams(paramsThese);
#endif
}
