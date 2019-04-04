   {
#ifdef STACK_PARAMS
      real8 bulkMod = -1.0 ; // dummy 
#else
      real8 bulkMod = elastN.getBulkMod() ;
#endif
      real8 gamma = 1.7 ;
      real8 ec0 = -cvav*300. ; 
      const std::vector<real8> paramsThese{ rho0, bulkMod, cvav, gamma, ec0 };
#ifdef STACK_PARAMS
      int nParamsEOS = paramsThese.size()-mmodel->nParamsEOSHave ; // nasty complexity to match what happens in matModel
      for ( int iP=0; iP<nParamsEOS; ++iP ) {
         params.push_back(paramsThese[mmodel->nParamsEOSHave+iP]) ;
      }
#else
      eos.setParams( paramsThese ) ;
#endif      
   }
