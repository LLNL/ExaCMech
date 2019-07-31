// -*-c++-*-

#ifndef ECMECH_SLIPGEOM_H
#define ECMECH_SLIPGEOM_H

#include "ECMech_core.h"

namespace ecmech {

class SlipGeomFCC
{
public:
   
   static const int nslip = 12 ;
   // constructor and destructor
   __ecmech_hdev__  SlipGeomFCC();
   __ecmech_hdev__ ~SlipGeomFCC(){};

   __ecmech_hdev__ inline const double* getP() const { return _P_ref_vec; } ;
   __ecmech_hdev__ inline const double* getQ() const { return _Q_ref_vec; } ;
   
private:
   double _P_ref_vec[ ecmech::ntvec * nslip ] ;
   double _Q_ref_vec[ ecmech::nwvec * nslip ] ;

};

} // namespace ecmech

#endif  // ECMECH_SLIPGEOM_H
