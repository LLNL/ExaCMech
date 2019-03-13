// -*-c++-*-

#ifndef ECMECH_SLIPGEOM_H
#define ECMECH_SLIPGEOM_H

namespace ecmech {

class SlipGeomFCCA
{
public:
   
   static const int nslip = 12 ;
   // constructor and destructor
   __ecmech_hdev__  SlipGeomFCCA();
   __ecmech_hdev__ ~SlipGeomFCCA(){};

   __ecmech_hdev__ inline const real8* const getP() const { return _P_ref_vec; } ;
   __ecmech_hdev__ inline const real8* const getQ() const { return _Q_ref_vec; } ;
   
private:
   real8 _P_ref_vec[ ecmech::ntvec * nslip ] ;
   real8 _Q_ref_vec[ ecmech::nwvec * nslip ] ;

}

} // namespace ecmech

#endif  // ECMECH_SLIPGEOM_H
