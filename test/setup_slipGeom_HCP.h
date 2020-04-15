{
   double covera = 1.59773226818;
   std::vector<double> paramsThese { covera };
#ifdef STACK_PARAMS
   params.insert(params.end(), paramsThese.begin(), paramsThese.end());
#else
   slipGeom.setParams(paramsThese);
#endif
}
