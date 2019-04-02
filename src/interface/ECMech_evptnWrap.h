// -*-c++-*-

#ifndef ecmech_evptnWrap_include
#define ecmech_evptnWrap_include

namespace ecmech {

//
// template on the specifics of the crystal model ;
// but have a base class so that the templating can stop here
//   
template< class SlipGeom, class Kinetics, class ThermoElastN, class EosModel >
class matModelEvptn : public matModelBase
{
public :
   
   // constructor 
   matModelEvptn()
      : matModelBase() {
      _slipGeom     = new SlipGeom() ;
      _kinetics     = new Kinetics(SlipGeom::nslip) ;
      _elastN       = new ThermoElastN() ;
      _eosModel     = new EosModel() ;

      // see n_rsv_matmod in F90 code
      _numHist = ecmech::evptn_iHistLbH + Kinetics::nH + SlipGeom::nslip ;
         
   };

   // destructor 
   ~matModelEvptn() {
      delete _slipGeom ;
      delete _kinetics ;
      delete _elastN ;
      delete _eosModel ;
   };

   void initFromParams(const std::vector<int>    & opts,
                       const std::vector<real8>  & pars,
                       const std::vector<string> & strs ) {

      const int nParamsEOSHave = 3 ; // number that get from 'elsewhere'
      const int nParamsEOS = EosModel::nParams-nParamsEOSHave ; 
      int nParams =
         2 + 1 + // rho0, cvav, tolerance
         Kinetics::nParams +
         ThermoElastN::nParams +
         nParamsEOS
      
      if ( pars.size() != nParams ) {
         ECMECH_FAIL(__func__,"wrong number of pars") ;
      }
      if ( opts.size() != 0 ) {
         ECMECH_FAIL(__func__,"wrong number of opts") ;
      }
      if ( strs.size() != 0 ) {
         ECMECH_FAIL(__func__,"wrong number of opts") ;
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
         const std::vector<real8> paramsThese(EosModel::nParams) ;
         paramsThese[0] = _rho0 ;
         paramsThese[1] = bulkMod ;
         paramsThese[2] = _cvav ;
         std::copy(parsIt, parsIt+nParamsEOS, paramsThese.begin()+nParamsEOSHave) ;
         
         _eosModel.setParams( paramsThese ) ; parsIt += nParamsEOS ;

         _eosModel.getEV0( _e0, _v0 ) ;
      }

      int iParam = parsIt - params.begin();
      if ( iParam != nParams ) {
         ECMECH_FAIL(__func__,"wrong number of params") ;
      }

      //////////////////////////////
      
      _rhvNames.clear();
      _rhvVals.clear();
      _rhvPlot.clear();
      _rhvState.clear();

      _rhvNames.push_back("shrate_eff") ; _rhvVals.push_back(0.) ; _rhvPlot.push_back(true) ; _rhvState.push_back(true) ;
      _rhvNames.push_back("shr_eff") ; _rhvVals.push_back(0.) ; _rhvPlot.push_back(true) ; _rhvState.push_back(true) ;
      //
      for ( int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec ) {
         std::string name = "t" << iTvec+1 ;
         _rhvNames.push_back(name) ; _rhvVals.push_back(0.) ; _rhvPlot.push_back(true) ; _rhvState.push_back(true) ;
      }
      //
      {
         real8 qVal = 1.0 ;
         for ( int iQ = 0; iQ < ecmech::qdim; ++iQ ) {
            std::string name = "quat_" << iQ+1 ;
            _rhvNames.push_back(name) ; _rhvVals.push_back(qVal) ; _rhvPlot.push_back(true) ; _rhvState.push_back(true) ;
            qVal = 0.0 ;
         }
      }
      //
      _kinetics.getHistInfo( _rhvNames, _rhvVals, _rhvPlot, _rhvState ) ;
      //
      for ( int iSlip = 0; iSlip < SlipGeom::nslip; ++iSlip ) {
         std::string name = "shrate_" << iSlip+1 ;
         _rhvNames.push_back(name) ; _rhvVals.push_back(0.) ; _rhvPlot.push_back(true) ; _rhvState.push_back(true) ;
      }
      //
      if ( _rhvNames.size() != _numHist ) {
         ECMECH_FAIL(__func__,"mismatch in numHist") ;
      }

   };
   
   void getParams(std::vector<int>    & opts,
                  std::vector<real8>  & pars,
                  std::vector<string> & strs ) const {
      // ...*** ;
      ECMECH_FAIL(__func__,"getParams not yet implemented") ;
   };

   void logParameters( std::ostringstream & oss ) {
      // ...*** ;
      ECMECH_FAIL(__func__,"logParameters not yet implemented") ;
   };

   void getResponse(const real8  & dt          , 
                    const real8  & curTime     ,
                    const real8  * defRateV    ,      
                    const real8  * spinV       ,         
                    const real8  * volRatioV   ,     
                          real8  * eIntV       ,         
                          real8  * stressSvecPV,  
                          real8  * histV       ,         
                          real8  * tkelvV      ,        
                          real8  * sddV        ,          
                          real8  * mtanSDV     ,      
                    const int    & nHist       ,     
                    const int    & nPassed      ) const {
      
      real8 *mtanSDThis = nullptr ;
      if ( nHist != _numHist ) {
         ECMECH_FAIL(__func__,"numHist mismatch") ;         
      }
      //
      for ( int i=0; i<nPassed; ++i ) {
         if ( mtanSDV != nullptr ) {
            mtanSDThis = &mtanSDV[ecmech::nsvec2*i];
         }
         getResponseEvptnSngl<SlipGeom, Kinetics, ThermoElastN, EosModel>( _slipGeom, _kinetics, _elastN, _eosModel,
                                                                           dt,
                                                                           curTime,
                                                                           _tolerance,
                                                                           &defRateV[ecmech::nsvp*i],
                                                                           &spinV[ecmech::ndim*i],
                                                                           &volRatioV[ecmech::nvr*i],
                                                                           &eIntV[ecmech::ne*i],
                                                                           &stressSvecPV[ecmech::nsvp*i],
                                                                           &histV[_numHist*i],
                                                                           tkelvV[i],
                                                                           &sddV[ecmech::nsdd*i],
                                                                           mtanSDThis );
      }
   };

   void getHistInfo(std::vector<std::string> & names,
                    std::vector<real8>       & vals,
                    std::vector<bool>        & plot,
                    std::vector<bool>        & state) {

      if ( _rhvNames.size() != _numHist ) {
         ECMECH_FAIL(__func__, "have not yet set up history information") ;
      }
      names.resize(_numHist); std::copy( _rhvNames.begin(), _rhvNames.end(), names.begin() ) ;
      vals.resize (_numHist); std::copy(  _rhvVals.begin(),  _rhvVals.end(),  vals.begin() ) ;
      plot.resize (_numHist); std::copy(  _rhvPlot.begin(),  _rhvPlot.end(),  plot.begin() ) ;
      state.resize(_numHist); std::copy( _rhvState.begin(), _rhvState.end(), state.begin() ) ;

   };
   
   int getNumHist( ) const {
      if ( _numHist < 0 ) {
         ECMECH_FAIL(__func__,"numHist < 0") ;
      }
      return _numHist ;
   } ;

private:
   
   SlipGeom *_slipGeom ;
   Kinetics *_kinetics ;
   ThermoElastN *_elastN ;
   EosModel *_eosModel ;

   int _numHist ;

   real8 _tolerance ;

   std::vector<std::string> _rhvNames;
   std::vector<real8>       _rhvVals;
   std::vector<bool>        _rhvPlot;
   std::vector<bool>        _rhvState;
   
} // namespace ecmech

#endif // ecmech_evptnWrap_include
