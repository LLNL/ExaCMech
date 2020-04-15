{
   double c11 = 300e-2, c12 = 100e-2, c44 = 100e-2;
   std::vector<double> paramsThese { c11, c12, c44 };
#ifdef STACK_PARAMS
   params.insert(params.end(), paramsThese.begin(), paramsThese.end());
#else
   elastN.setParams(paramsThese);
#endif
}
