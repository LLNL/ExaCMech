{
   std::vector<double> paramsThese; // no parameters for FCC
#ifdef STACK_PARAMS
   params.insert(params.end(), paramsThese.begin(), paramsThese.end());
#else
   slipGeom.setParams(paramsThese);
#endif
}
