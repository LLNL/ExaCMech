// -*-c++-*-

#ifndef ecmech_evptnWrap_include
#define ecmech_evptnWrap_include

#include "ECMech_core.h"

#ifdef __cuda_host_only__
#include <sstream>
#endif

#include "ECMech_matModelBase.h"
#include "ECMech_evptn.h"

namespace ecmech {

namespace evptn {

//
// template on the specifics of the crystal model ;
// but have a base class so that the templating can stop here
//   
template< class SlipGeom, class Kinetics, class ThermoElastN, class EosModel >
class matModel : public matModelBase
{
public :

   static const int iHistLbGdot = NumHist<SlipGeom, Kinetics, ThermoElastN, EosModel>::iHistLbGdot ;
   static const int numHist     = NumHist<SlipGeom, Kinetics, ThermoElastN, EosModel>::numHist ;

   // this are assumed to go in first
   static const int nParamsEOSHave = 3 ; // number that get from 'elsewhere'
         
   // constructor 
   matModel()
      : matModelBase(),
        _kinetics(SlipGeom::nslip),
        _outputLevel(0) {
   };

   __ecmech_hdev__
   void setOutputLevel( int outputLevel ) { _outputLevel = outputLevel ; } ;
   
   using matModelBase::initFromParams ;
   __ecmech_hdev__
   void initFromParams(const std::vector<int>         & opts,
                       const std::vector<real8>       & pars,
                       const std::vector<std::string> & strs ) {

      const int nParamsEOS = EosModel::nParams-nParamsEOSHave ; 
      int nParams =
         2 + 1 + // rho0, cvav, tolerance
         Kinetics::nParams +
         ThermoElastN::nParams +
         nParamsEOS ;
      
      if ( pars.size() != (unsigned int)nParams ) {
         ECMECH_FAIL(__func__,"wrong number of pars") ;
      }
      if ( opts.size() != 0 ) {
         ECMECH_FAIL(__func__,"wrong number of opts") ;
      }
      if ( strs.size() != 0 ) {
         ECMECH_FAIL(__func__,"wrong number of strs") ;
      }
      
      std::vector<real8>::const_iterator parsIt = pars.begin() ;
      
      _rho0 = *parsIt; ++parsIt;
      _cvav = *parsIt; ++parsIt;

      _tolerance = *parsIt; ++parsIt;      

      {
         const std::vector<real8> paramsThese(parsIt, parsIt+ThermoElastN::nParams) ;
         _elastN.setParams( paramsThese ) ; parsIt += ThermoElastN::nParams ;
      }
      {
         const std::vector<real8> paramsThese(parsIt, parsIt+Kinetics::nParams) ;
         _kinetics.setParams( paramsThese ) ; parsIt += Kinetics::nParams ;
      }
      {
         real8 bulkMod = _elastN.getBulkMod() ;
         std::vector<real8> paramsThese(EosModel::nParams) ;
         paramsThese[0] = _rho0 ;
         paramsThese[1] = bulkMod ;
         paramsThese[2] = _cvav ;
         std::copy(parsIt, parsIt+nParamsEOS, paramsThese.begin()+nParamsEOSHave) ;
         
         _eosModel.setParams( paramsThese ) ; parsIt += nParamsEOS ;

         _eosModel.getEV0( _e0, _v0 ) ;
      }

      int iParam = parsIt - pars.begin();
      if ( iParam != nParams ) {
         ECMECH_FAIL(__func__,"wrong number of params") ;
      }

#ifdef __cuda_host_only__
      //////////////////////////////
      
      _rhvNames.clear();
      _rhvVals.clear();
      _rhvPlot.clear();
      _rhvState.clear();

      _rhvNames.push_back("shrate_eff") ; _rhvVals.push_back(0.) ; _rhvPlot.push_back(true) ; _rhvState.push_back(true) ; // iHistA_shrateEff
      _rhvNames.push_back("shr_eff") ; _rhvVals.push_back(0.) ; _rhvPlot.push_back(true) ; _rhvState.push_back(true) ; // iHistA_shrEff
      _rhvNames.push_back("n_feval") ; _rhvVals.push_back(0.) ; _rhvPlot.push_back(true) ; _rhvState.push_back(true) ; // iHistA_nFEval
      // numHistAux
      //
      for ( int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec ) {
         std::ostringstream os;
         os << "t" << iTvec+1 ;
         _rhvNames.push_back(os.str()) ; _rhvVals.push_back(0.) ; _rhvPlot.push_back(true) ; _rhvState.push_back(true) ;
      }
      //
      {
         real8 qVal = 1.0 ;
         for ( int iQ = 0; iQ < ecmech::qdim; ++iQ ) {
            std::ostringstream os;
            os << "quat_" << iQ+1 ;
            _rhvNames.push_back(os.str()) ; _rhvVals.push_back(qVal) ; _rhvPlot.push_back(true) ; _rhvState.push_back(true) ;
            qVal = 0.0 ;
         }
      }
      //
      _kinetics.getHistInfo( _rhvNames, _rhvVals, _rhvPlot, _rhvState ) ;
      //
      for ( int iSlip = 0; iSlip < SlipGeom::nslip; ++iSlip ) {
         std::ostringstream os;
         os << "shrate_" << iSlip+1 ;
         _rhvNames.push_back(os.str()) ; _rhvVals.push_back(0.) ; _rhvPlot.push_back(true) ; _rhvState.push_back(true) ;
      }
      //
      if ( _rhvNames.size() != numHist ) {
         ECMECH_FAIL(__func__,"mismatch in numHist") ;
      }
#endif
      
   };
   
   __ecmech_host__
   void getParams(std::vector<int>         & opts,
                  std::vector<real8>       & pars,
                  std::vector<std::string> & strs ) const {
      // ...*** ;
      opts.clear() ;
      strs.clear() ;
      pars.clear() ;
      ECMECH_FAIL(__func__,"getParams not yet implemented") ;
   };

   __ecmech_host__
   void logParameters( std::ostringstream & oss ) const {
      // ...*** ;
      oss << "evptn constitutive model" << std::endl ;
      ECMECH_FAIL(__func__,"logParameters not yet implemented") ;
   };

   using matModelBase::getResponse ;
   __ecmech_hdev__
   void getResponse(const real8  & dt          , 
                    const real8  * defRateV    ,      
                    const real8  * spinV       ,         
                    const real8  * volRatioV   ,     
                          real8  * eIntV       ,         
                          real8  * stressSvecPV,  
                          real8  * histV       ,         
                          real8  * tkelvV      ,        
                          real8  * sddV        ,          
                          real8  * mtanSDV     ,      
                    const int    & nPassed      ) const {
      
      real8 *mtanSDThis = nullptr ;
      //
      for ( int i=0; i<nPassed; ++i ) {
         if ( mtanSDV != nullptr ) {
            mtanSDThis = &mtanSDV[ecmech::nsvec2*i];
         }
         getResponseSngl<SlipGeom, Kinetics, ThermoElastN, EosModel>( _slipGeom, _kinetics, _elastN, _eosModel,
                                                                      dt,
                                                                      _tolerance,
                                                                      &defRateV[ecmech::nsvp*i],
                                                                      &spinV[ecmech::ndim*i],
                                                                      &volRatioV[ecmech::nvr*i],
                                                                      &eIntV[ecmech::ne*i],
                                                                      &stressSvecPV[ecmech::nsvp*i],
                                                                      &histV[numHist*i],
                                                                      tkelvV[i],
                                                                      &sddV[ecmech::nsdd*i],
                                                                      mtanSDThis,
                                                                      _outputLevel );
      }
   };

#ifdef __cuda_host_only__
   __ecmech_host__
   void getHistInfo(std::vector<std::string> & names,
                    std::vector<real8>       & vals,
                    std::vector<bool>        & plot,
                    std::vector<bool>        & state) const {

      if ( _rhvNames.size() != numHist ) {
         ECMECH_FAIL(__func__, "have not yet set up history information") ;
      }
      names.resize(numHist); std::copy( _rhvNames.begin(), _rhvNames.end(), names.begin() ) ;
      vals.resize (numHist); std::copy(  _rhvVals.begin(),  _rhvVals.end(),  vals.begin() ) ;
      plot.resize (numHist); std::copy(  _rhvPlot.begin(),  _rhvPlot.end(),  plot.begin() ) ;
      state.resize(numHist); std::copy( _rhvState.begin(), _rhvState.end(), state.begin() ) ;

   };
#endif
   
   int getNumHist( ) const {
      return numHist ;
   } ;

private:
   
   SlipGeom _slipGeom ;
   Kinetics _kinetics ;
   ThermoElastN _elastN ;
   EosModel _eosModel ;

   real8 _tolerance ;
   int _outputLevel ;

#ifdef __cuda_host_only__
   // TO_DO : it is a good idea to have these host-only?
   std::vector<std::string> _rhvNames;
   std::vector<real8>       _rhvVals;
   std::vector<bool>        _rhvPlot;
   std::vector<bool>        _rhvState;
#endif
   
}; // class matModel 

} // namespace ecmech

} //  namespace evptn 

#endif // ecmech_evptnWrap_include
