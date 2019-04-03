#include "SNLS_TrDLDenseG.h"

#include "ECMech_cases.h"
#include "ECMech_evptnWrap.h"

int main(int argc, char *argv[])
{
   int outputLevel = 1 ;
   if ( argc > 1 ) {
      outputLevel = atoi(argv[1]) ;
   }

   using namespace ecmech ;

   typedef Kin_KMBalD_FFF Kinetics  ;
   typedef SlipGeomFCC SlipGeom ;
   typedef evptn::ThermoElastNCubic ThermoElastN ;
   typedef EosModelConst<false> EosModel ;

   ecmech::evptn::matModel<SlipGeom, Kinetics, ThermoElastN, EosModel> mmodel ;

#include "setup_base.h"
   std::vector<int>           opts ; // none
   std::vector<std::string>   strs ; // none
   std::vector<real8>         params{ rho0, cvav, tolerance } ;
   {
   
      // NOTE : not using these instances other than as a hack to provide parameters
      Kinetics kinetics(SlipGeom::nslip) ;
#include "setup_kin_KMBalD_FFF.h"

      ThermoElastN elastN ;
#include "setup_elastn.h"
      
      EosModel eos ;
#include "setup_eos.h"
      
      kinetics.getParams( params ) ;
      elastN.getParams( params ) ;

      std::vector<real8> eosParams ;
      eos.getParams( eosParams ) ;
      int nParamsEOS = eosParams.size()-mmodel.nParamsEOSHave ; // nasty complexity to match what happens in matModel
      for ( int iP=0; iP<nParamsEOS; ++iP ) {
         params.push_back(eosParams[mmodel.nParamsEOSHave+iP]) ;
      }
   }
   //
   mmodel.initFromParams( opts, params, strs ) ;
   mmodel.complete() ;

   std::vector<real8>       hist_vec ;
   {
      std::vector<std::string> names ;
      std::vector<bool>        plot ;
      std::vector<bool>        state ;
      mmodel.getHistInfo(names, hist_vec, plot, state ) ;
   }
   real8* hist = &(hist_vec[0]);
   static const int iHistLbGdot = mmodel.iHistLbGdot ;
   real8* gdot = &(hist[iHistLbGdot]) ; 
   
   const int nPassed = 1 ; // do not change this without some serious changes elsewhere
   mmodel.setOutputLevel( outputLevel ) ;
   //
   {
#include "setup_conditions.h"      
      real8 eInt[ecmech::ne] = {0.0} ;
      real8 stressSvecP[ecmech::nsvp] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                         0.0};
      real8 tkelv[nPassed] ;
      real8 sdd[ecmech::nsdd*nPassed] ;
      mmodel.getResponse( dt, d_svec_kk_sm, w_veccp_sm, volRatio,
                          eInt, stressSvecP, hist, tkelv, sdd, nullptr,
                          nPassed ) ;

      std::cout << "Function evaluations: " << hist[evptn::iHistA_nFEval] << std::endl ;
   
      std::cout << "Slip system shearing rates : " ;
      printVec<SlipGeom::nslip>(gdot, std::cout) ;
   }
      
}

