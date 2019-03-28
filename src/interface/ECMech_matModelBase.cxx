#include "ECMech_matModelBase.h"

real8
ecmech::matModelBase::getRhoRef() const
{
  if ( _rho0 < 0.0 ) { // want to be able to call this before _complete
     mslib::UninitializedException ex;
     ex.setNote("rho0 does not appear to have been set.");
     ex.addFileInfo(__FILE__, __LINE__, "getRhoRef");
     throw ex;
  }

  real8 retval = _rho0;
  return retval;
}
