// -*-c++-*-

#ifndef ECMech_matModelBase_include
#define ECMech_matModelBase_include

#include <string>
#include <vector>

#include "ECMech_core.h"

namespace ecmech {
   // **************** //
   // Class Definition //
   // **************** //

   class matModelBase
   {
      protected:
         bool  _complete;
         double _rho0, _cvav, _v0, _e0, _bulkRef;

         // constructor
         __ecmech_host__
         matModelBase() :
            _complete(false),
            _rho0(-1.0),
            _cvav(-1.0),
            _v0(-1.0),
            _e0(-1.0),
            _bulkRef(-1.0)
         {};

      public:
         // deconstructor
         __ecmech_host__
         virtual ~matModelBase() {};

         virtual void initFromParams(const std::vector<int>& opts,
                                     const std::vector<double>& pars,
                                     const std::vector<std::string>& strs,
                                     void* call_back = nullptr) = 0;

         virtual void getParams(std::vector<int>& opts,
                                std::vector<double>& pars,
                                std::vector<std::string>& strs) const = 0;

         /**
          * @brief log parameters, including history information; more human-readable than getParams output
          */
         virtual void logParameters(std::ostringstream& oss) const = 0;

         /**
          * @brief Request response information for a group of host-code
          * points. The stress and history variables are updated over finite
          * time step of size dt.
          *
          * For arguments tha are of length x*nPassed, indexing is fastest along x
          *
          * The interface is always for 3D deformation -- in 2D some of the
          * deformation rate (defRateV) and spin (spinV) will be zero. For
          * anisotropic materials, the stress response can still be fully
          * populated with non-zeros. If the stress is used to encode state
          * for the given material model (which depends on the details of the
          * implementation), then all stress components should still be
          * stored even in 2D.
          *
          * Deviatoric stress and history variables are generally updated in
          * a first-order time integration (implicit; eg,
          * backward-Euler). Equation of state may be done either third-order
          * or first-order. When closed-form analytic ODE solutions are
          * available, some of the model implementations use that closed-form
          * solution for the state update. But those solutions assumes a
          * constant value of some arguments over the time step.
          *
          * This interface may not therefore be compatible with certain
          * integration schemes that want an instantaneous rate of evolution
          * for state variables and stress components. The library is meant
          * for operating over a finite time step size dt. Although it is
          * possible to do things like make two getResponse type calls to
          * evaluate at both the half-step and full-step dt.
          *
          * @param[in] nPassed : Number of host-code points passed
          *
          * @param[in] dt : Time step size
          *
          * @param[in] defRateV : Components of the deformation rate (symmetric part of the velocity gradient)
          * length nsvp*nPassed
          * Voigt ordering
          * first six components are the deviatoric part (zero trace)
          * along nsvp :
          *    [dxx, dyy, dzz, dyz, dxz, dxy, vdov]
          * for vdov, see volRatio; but note that sometimes other expressions are used for vdov;
          *    for example in implicit global time stepping
          *
          * @param[in] spinV : Components of the spin (skew part of the velocity gradient)
          * length ndim*nPassed
          * Voigt ordering
          *    wxx = (L32-L23)/2
          *    wyy = (L13-L31)/2
          *    wzz = (L21-L12)/2
          * along ndim :
          *    [wxx, wyy, wzz]
          *
          * @param volRatio[in] : information about volume evolution
          * length nvr*nPassed
          * along nvr :
          *    [vOld, vNew, vdov, delv]
          * vOld -- relative volume at beginning of time step
          * vNew -- relative volume at end of time step
          * vdov = delv / (dt * 0.5*(vNew+vOld)) -- volumetric strain rate
          * delv = vNew - vOld -- increment in relative volume
          *
          * @param eIntV[in,out] : Internal energy per reference volume
          * length ne*nPassed
          * along ne :
          *    [eTotal, cold, eQ, etherms, ?, ?, deltrh, ?, deltz, eMelt]
          * on input, all but the eTotal (first entry) should be zero, and eTotal is beginning-of-step;
          * on output, eTotal is updated to end-of-step
          *
          * @param stressSvecPV[in,out] : Cauchy stress components
          * length nsvp*nPassed
          * first six components are the deviatoric part (zero trace)
          * along nsvp :
          *    [sx, sy, sz, tyz, txz, txy, p]
          * beginning-of-step on input, end-of-step on output
          *
          * @param histV[in,out] : History (eg, state) variables
          * length numHist*nPassed if hIndx is NULL
          * along numHist : order is as indicated by getHistInfo
          * beginning-of-step on input, end-of-step on output
          *
          * @param tkelvV[out] : end-of-step temperature
          * length nPassed
          *
          * @param sddV[out] : other output quantities
          * length nsdd*nPassed
          * along nsdd, index by i_sdd_* (eg, i_sdd_gmod)
          *
          * @param mtanSDV[out] : tangent stiffness
          * length nsvec*nsvec*nPassed
          *
          * @param hIndx[in] : index array for histV, or NULL for close-packed
          * length nPassed, or NULL
          *    pass a negative value for a host-code point that is to be 'skipped'
          */

         __ecmech_host__
         virtual void getResponse(const double & dt,
                                  const double * defRateV,
                                  const double * spinV,
                                  const double * volRatioV,
                                  double * eIntV,
                                  double * stressSvecPV,
                                  double * histV,
                                  double * tkelvV,
                                  double * sddV,
                                  double * mtanSDV,
                                  const int & nPassed) const = 0;

         /**
          * @brief
          * Get the history variable information.
          */
         virtual void getHistInfo(std::vector<std::string>& histNames,
                                  std::vector<double>& initVals,
                                  std::vector<bool>& plot,
                                  std::vector<bool>& state) const = 0;


         /**
          * @brief
          * Get number of history variables
          */
         __ecmech_hdev__
         virtual int getNumHist( ) const = 0;

         /**
          *  @brief
          *  Set the accelerator to be used for getResponse
          *  The available options are CPU, OPENMP, and CUDA.
          */
         __ecmech_host__
         virtual void setAccelerator(ecmech::Accelerator accel) = 0;

         /**
          * @brief Get the reference density
          */
         __ecmech_hdev__
         virtual double getRhoRef() const {
            if (_rho0 < 0.0) { // want to be able to call this before _complete
               ECMECH_FAIL(__func__, "rho0 does not appear to have been set");
            }
            return _rho0;
         };


         /**
          * @brief
          * May end up requiring this to be called before the model may be used; and probably want to redefine this
          */
         __ecmech_hdev__
         virtual void complete() { _complete = true; };

         /**
          * @brief
          * Return whether or not complete has been called
          */
         __ecmech_hdev__
         virtual bool isComplete() { return _complete; };
   }; // class matModelBase
} // ecmech namespace

#endif // ECMech_matModelBase_include
