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
      _numHist = 
         1 + // effective shearing rate
         1 + // accumulated shear
         ecmech::ntvec + ecmech::qdim + Kinetics::nH +
         SlipGeom::nslip ;
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

      int nParams = 4 + Kinetics::nParams + ThermoElastN::nParams + EosModel::nParams ;
      
      if ( pars.size() < nParams ) {
         ECMECH_FAIL(__func__,'wrong number of pars') ;
      }
      if ( opts.size() != 0 ) {
         ECMECH_FAIL(__func__,'wrong number of opts') ;
      }
      if ( strs.size() != 0 ) {
         ECMECH_FAIL(__func__,'wrong number of opts') ;
      }
      
      _rhvNames.clear();
      _rhvVals.clear();
      _rhvPlot.clear();
      _rhvState.clear();
      
      std::vector<real8>::const_iterator parsIt = pars.begin() ;
      
      _rho0 = *parsIt; ++parsIt;
      _v0   = *parsIt; ++parsIt;
      _e0   = *parsIt; ++parsIt;
      _cvav = *parsIt; ++parsIt;

      _elastN.setParams( parsIt )   ; parsIt += ThermoElastN::nParams ;
      _kinetics.setParams( parsIt ) ; parsIt += Kinetics::nParams ;
      _eosModel.setParams( parsIt ) ; parsIt += EosModel::nParams ;

      _kinetics.getHistInfo( _rhvNames, _rhvVals, _rhvPlot, _rhvState ) ;
...
      _rhvNames.clear();
      _rhvVals.clear();
      _rhvPlot.clear();
      _rhvState.clear();

         int nPars = parsIt - pars.begin();
   if ( ((int) pars.size() - nPars) == (int)(rhvVals.size()) ) {
      // assume remaining entries are initial history variable values
      for ( unsigned int iH=0; iH < rhvVals.size(); iH++ ) {
         rhvVals[iH] = *parsIt; ++parsIt;
      }
   }      
   else if ( nPars != (int) pars.size() ) 
   {
      std::ostringstream os;
      os << "params should have had length " << nPars << " but instead was " << pars.size() ;
      std::string msg(os.str());      

      MS_FAIL(msg.c_str());

   }
   };
   
   void getParams(std::vector<int>    & opts,
                  std::vector<real8>  & pars,
                  std::vector<string> & strs ) const {
      // ...*** ;
      ECMECH_FAIL(__func__,'getParams not yet implemented') ;
   };

   void logParameters( std::ostringstream & oss ) {
      // ...*** ;
      ECMECH_FAIL(__func__,'logParameters not yet implemented') ;
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
                    const real8  * elemLenV    ,      
                          real8  * sddV        ,          
                          real8  * mtanSDV     ,      
                    const int    & nHist       ,     
                    const int    & nPassed      ) const {
      ...*** ;
   };

   void getHistInfo(std::vector<std::string> & histNames,
                    std::vector<real8>       & initVals,
                    std::vector<bool>        & plot,
                    std::vector<bool>        & state) {
      ...*** ;
   };
   
   int getNumHist( ) const {
      if ( _numHist < 0 ) {
         ECMECH_FAIL(__func__,'numHist < 0') ;
      }
      return _numHist ;
   } ;

private:
   
   SlipGeom *_slipGeom ;
   Kinetics *_kinetics ;
   ThermoElastN *_elastN ;
   EosModel *_eosModel ;

   int _numHist ;

   std::vector<std::string> _rhvNames;
   std::vector<real8>       _rhvVals;
   std::vector<bool>        _rhvPlot;
   std::vector<bool>        _rhvState;
   
} // namespace ecmech

#endif // ecmech_evptnWrap_include
