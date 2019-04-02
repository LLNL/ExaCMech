// -*-c++-*-

#ifndef ECMECH_EOS_SIMPLE_H
#define ECMECH_EOS_SIMPLE_H

#include "ECMech_core.h"

#include <string>
#include <vector>

namespace ecmech {

template< bool isothermal >
class EosModelConst
{
public:
   static const int nParams = 4 ;
   
   // constructor
   __ecmech_hdev__
   EosModelConst() {};

   __ecmech_host__
   inline void setParams( const std::vector<real8> & params // const real8* const params
                          ) {

      std::vector<real8>::const_iterator parsIt = params.begin();

      //////////////////////////////

      _rho0    = *parsIt; ++parsIt ;
      _bulkMod = *parsIt; ++parsIt ;
      _cvav    = *parsIt; ++parsIt ;
      _gamma   = *parsIt; ++parsIt ;
      _ec0     = *parsIt; ++parsIt ;

      _dtde = one / _cvav ;
      _tK0  = -_ec0 * _dtde ;
      
      //////////////////////////////

      int iParam = parsIt - params.begin();
      if ( iParam != nParams ) {
         ECMECH_FAIL(__func__,"iParam != nParams");
      }
      
   };

   __ecmech_hdev__
   inline void evalPT( real8 &p,
                       real8 &tK,
                       real8  v,
                       real8  e ) const {
      
      real8 mu = one / v - one ;
      
      if ( isothermal ) {
         p  = _bulkMod * mu ;
         tK = _tK0 ;
      }
      else {
         p  = _bulkMod * mu + _gamma * e ;
         tK = _tK0 + e * _dtde ;
      }
      
   }
   
   __ecmech_hdev__
   inline void evalPTDiff( real8 &p,
                           real8 &tK,
                           real8 &bulkNew,
                           real8 &dpde,
                           real8  v,
                           real8  e ) const {

      real8 eta = one / v ;
      real8 mu  = eta - one ;

      if ( isothermal ) {
         p  = _bulkMod * mu ;
         tK = _tK0 ;
         dpde = zero ;
      }
      else {
         p  = _bulkMod * mu + _gamma * e ;
         tK = _tK0 + e * _dtde ;
         dpde = _gamma ;
      }
      bulkNew = _bulkMod * eta ;
      
   }
   
   __ecmech_hdev__
   inline void getEV0( real8 &e0,
                       real8 &v0 ) const {
      e0 = 0.0 ;
      v0 = 1.0 ;
   }
   
private:

   // parameters
   real8 _rho0, _bulkMod, _gamma, _ec0, _cvav ;

   // derived from parameters
   real8 _dtde, _tK0 ;

}; // class KineticsVocePL

template< class EosModel >
__ecmech_hdev__
inline void updateSimple( const EosModel& eos,
                          real8 &p,
                          real8 &tK,
                          real8 &eNew,
                          real8 &bulkNew, 
                          const real8* volRatio,
                          real8  eOld,
                          real8  pOld) {

   real8 vOld = volRatio[0] ;
   real8 vNew = volRatio[1] ;
   real8 delv = volRatio[3] ;

   eNew = eOld - delv * pOld ;

   real8 dpde ;
   eos.evalPTDiff( p, tK, bulkNew, dpde, vNew, eNew ) ;
   // dpdv = - bulkNew / vNew ; // not needed
   
   bulkNew = bulkNew + dpde * pOld * vNew ;
      
}
   
   
} // namespace ecmech

#endif  // ECMECH_EOS_SIMPLE_H
