{
   std::vector<double> paramsThese; // no parameters for FCC or BCC
#ifdef STACK_PARAMS
   params.insert(params.end(), paramsThese.begin(), paramsThese.end());
#else
   slipGeom.setParams(paramsThese);
#endif
}
