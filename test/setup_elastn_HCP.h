{
   double c11 = 162.4e-2, c12 = 92e-2, c13 = 69e-2, c33 = 180.7e-2, c44 = 46.7e-2;
   double g_vecd2 = 0.0;
   std::vector<double> paramsThese { c11, c12, c13, c33, c44, g_vecd2 };
#ifdef STACK_PARAMS
   params.insert(params.end(), paramsThese.begin(), paramsThese.end());
#else
   elastN.setParams(paramsThese);
#endif
}
