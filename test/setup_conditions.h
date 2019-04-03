   real8 d_svec_kk_sm[ecmech::nsvp] = {-0.5, -0.5, 1.0,
                                       0.0, 0.0, 0.0,
                                       0.0} ;
   vecsVsa<ecmech::nsvp>(d_svec_kk_sm, sqr2b3) ; // scale so that d_vecd_sm comes out to unit magnitude
   //
   real8 d_vecd_sm[ecmech::ntvec] ;
   svecToVecd(d_vecd_sm, d_svec_kk_sm) ;
   
   real8 dt = 1e-1 ;

   real8 w_veccp_sm[ecmech::nwvec] = {0.0, 0.0, 0.5} ;

   real8 volRatio[ecmech::nvr] = {1.0, 1.0, 0.0, 0.0} ;
