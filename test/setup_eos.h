   {
      real8 bulkMod = elastN.getBulkMod() ;
      real8 gamma = 1.7 ;
      real8 ec0 = -cvav*300. ; 
      const std::vector<real8> params{ rho0, bulkMod, cvav, gamma, ec0 };
      eos.setParams( params ) ;
   }
