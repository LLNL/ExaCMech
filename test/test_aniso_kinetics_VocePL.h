// -*-c++-*-

#ifndef TEST_ANISO_KINETICS_VOCEPL_H
#define TEST_ANISO_KINETICS_VOCEPL_H

#include <cassert>
#include <cmath>

namespace ecmech {
   /**
    * slip and hardening kinetics
    *
    * power-law slip kinetics with Voce hardening law -- meant to be about as simple as it gets
    */
   template<bool nonlinear,
            int nslip>
   class KineticsAnisoVocePL
   {
      public:
         static const int nH = nslip;
         static const int _nslip = nslip;
         static const int nParams = 3 + 5 + nH + (nonlinear ? 1 : 0);
         static const int nVals = nslip;
         static const int nEvolVals = 2;

         // constructor
         __ecmech_hdev__
         KineticsAnisoVocePL(int nslip_) {};
         // deconstructor
         __ecmech_hdev__
         ~KineticsAnisoVocePL() {}

         __ecmech_host__
         inline void setParams(const std::vector<double> & params // const double* const params
                               ) {
            std::vector<double>::const_iterator parsIt = params.begin();

            //////////////////////////////
            // power-law stuff

            _mu = *parsIt; ++parsIt;
            _xm = *parsIt; ++parsIt;
            _gam_w = *parsIt; ++parsIt;

            // CALL fill_power_law(pl)
            // xmm  = xm - one ;
            _xnn = one / _xm;
            _xn = _xnn - one;
            // xMp1 = xnn + one
            //
            // CALL set_t_min_max(pl)
            _t_min = pow(ecmech::gam_ratio_min, _xm);
            _t_max = pow(ecmech::gam_ratio_ovf, _xm);

            //////////////////////////////
            // Voce hardening stuff

            _h0 = *parsIt; ++parsIt;
            for(int i = 0; i < nslip; i++)
            {
               _tausi[i] = *parsIt; ++parsIt;
            }
            _taus0 = *parsIt; ++parsIt;

            if (nonlinear) {
               _xmprime = *parsIt; ++parsIt;
               _xmprime1 = _xmprime - one;
            }
            else {
               _xmprime = one;
               _xmprime1 = zero;
            }
            _xms = *parsIt; ++parsIt;
            _gamss0 = *parsIt; ++parsIt;

            //////////////////////////////
            // nH

            _hdn_init = *parsIt; ++parsIt;

            //////////////////////////////

            int iParam = parsIt - params.begin();
            assert(iParam == nParams);
         }

         __ecmech_host__
         inline void getParams(std::vector<double> & params
                               ) const {
            // do not clear params in case adding to an existing set
            int paramsStart = params.size();

            //////////////////////////////
            // power-law stuff

            params.push_back(_mu);
            params.push_back(_xm);
            params.push_back(_gam_w);

            //////////////////////////////
            // Voce hardening stuff

            params.push_back(_h0);
            for(int i = 0; i < _nslip; i++) {
               params.push_back(_tausi[i]);
            }
            params.push_back(_taus0);
            params.push_back(_xms);
            params.push_back(_gamss0);

            //////////////////////////////
            // nH

            params.push_back(_hdn_init);

            //////////////////////////////

            int iParam = params.size() - paramsStart;
            assert(iParam == nParams);
         }

         __ecmech_host__
         void getHistInfo(std::vector<std::string> & names,
                          std::vector<double>       & init,
                          std::vector<bool>        & plot,
                          std::vector<bool>        & state) const {
            names.push_back("h");
            init.push_back(_hdn_init);
            plot.push_back(true);
            state.push_back(true);
         }

      private:

         // static const _nXnDim = nH*nH ; // do not bother

         //////////////////////////////
         // power-law stuff

         // parameters
         double _mu; // may evetually set for current conditions
         double _xm;
         double _gam_w; // pl%adots, adots0

         // derived from parameters
         double _t_max, _t_min, _xn, _xnn;

         //////////////////////////////
         // Voce hardening stuff

         double _h0, _taus0, _xms, _gamss0;
         double _tausi[_nslip];
         double _xmprime, _xmprime1;

         //////////////////////////////

         double _hdn_init;

      public:

         __ecmech_hdev__
         inline double getFixedRefRate(const double* const // vals, not used
                                       ) const
         {
            return _gam_w;
         }

         __ecmech_hdev__
         inline
         double
         getVals(double* const vals,
                 double, // p, not currently used
                 double, // tK, not currently used
                 const double* const h_state
                 ) const
         {
            double mVals = ecmech::zero;
            for (int iSlip = 0; iSlip < _nslip; ++iSlip) {
               vals[iSlip] = h_state[iSlip]; // _gAll
               mVals += vals[iSlip];
               assert(vals[iSlip] > zero);
            }

            mVals /= _nslip;

            return mVals;
         }

         __ecmech_hdev__
         inline
         void
         evalGdots(double* const gdot,
                   double* const dgdot_dtau,
                   double* const dgdot_dg,
                   const double* const tau,
                   const double* const vals
                   ) const
         {
            for (int iSlip = 0; iSlip<this->_nslip; ++iSlip) {
               bool l_act;
               double gAll = vals[iSlip];
               this->evalGdot(gdot[iSlip], l_act, dgdot_dtau[iSlip], dgdot_dg[iSlip],
                              gAll,
                              tau[iSlip],
                              _mu // gss%ctrl%mu(islip)
                              );
            }
         }

         /**
          * see kinetics_pl_d
          */
         __ecmech_hdev__
         inline
         void
         evalGdot(
            double & gdot,
            bool  & l_act,
            double & dgdot_dtau, // wrt resolved shear stress
            double & dgdot_dg, // wrt slip system strength
#if MORE_DERIVS
            double & dgdot_dmu, // wrt shear modulus, not through g
            double & dgdot_dgamo, // wrt reference rate for thermal part
            double & dgdot_dgamr, // wrt reference rate for drag limited part
            double & dgdot_dtK, // wrt temperature, with other arguments fixed
#endif
            double   gIn,
            double   tau,
            double // mu, not currently used
#if MORE_DERIVS
            ,
            double   tK
#endif
            ) const
         {
            // zero things so that can more easily just return in inactive
            //// gdot_w = zero; gdot_r = zero; ! not used by l_linear or l_pl
            gdot = zero;
            //
            dgdot_dtau = zero;
            dgdot_dg = zero;
#if MORE_DERIVS
            dgdot_dmu = zero;
            dgdot_dgamo = zero;
            dgdot_dgamr = zero;
            dgdot_dtK = zero;
#endif
            l_act = false;

            double g_i = one / gIn; // assume have checked gIn>0 elsewhere
            double t_frac = tau * g_i; // has sign of tau
            double at = fabs(t_frac);

            if (at > _t_min) {
               //
               l_act = true;

               if (at > _t_max) {
                  // ierr = IERR_OVF_p
                  // set gdot big, evpp may need this for recovery
                  gdot = ecmech::gam_ratio_ovffx * _gam_w;
                  gdot = copysign(gdot, tau);
                  // do not set any of deriviatives (they are, in truth, zero)
               }
               else {
                  double abslog = log(at);
                  double blog = _xn * abslog;
                  double temp = _gam_w * exp(blog);

                  gdot = temp * t_frac;

                  dgdot_dtau = temp * _xnn * g_i; // note: always positive, = xnn * gdot/t
                  dgdot_dg = -dgdot_dtau * t_frac; // = - gdot * xnn * g_i
#if MORE_DERIVS
                  // dgdot_dmu   =  zero ; // already done
                  dgdot_dgamo = gdot / _gam_w;
                  // dgdot_dgamr = zero ; // already done
                  // dgdot_dtK   = zero ; // already done

                  // IF (pl%tdp%T_dep) THEN
                  // dgdot_dtK   =  gdot * pl%tdp%qoverr / (tK * tK)
                  // END IF
                  // IF (pl%tdp%T_dep_m) THEN !  .AND. pl%xm < XM_UB_p
                  // dxn_dtK = -(pl%xnn*pl%xnn) * pl%tdp%dxm_dtK_current
                  // dgdot_dtK = dgdot_dtK + &
                  // & gdot * abslog * dxn_dtK
                  // END IF
#endif
               }
            }
         } // evalGdot

         __ecmech_hdev__
         inline
         int
         updateH(double* const hs_u,
                 const double* const hs_o,
                 double dt,
                 const double* const gdot,
                 int outputLevel = 0) const
         {
            double hs_u_1[_nslip];
            int nFEvals = updateHN<KineticsAnisoVocePL>(this,
                                                        &hs_u_1[0], &hs_o[0], dt, gdot,
                                                        outputLevel);

            for(int i = 0; i < _nslip; i++) {
               hs_u[i] = hs_u_1[i];
            }

            return nFEvals;
         }

         __ecmech_hdev__
         inline
         void
         getEvolVals(double* const evolVals,
                     const double* const gdot
                     ) const
         {
            // recompute effective shear rate here versus using a stored value
            double shrate_eff = vecsssumabs_n(gdot, _nslip); // could switch to template if template class on _nslip

            double sv_sat = _taus0;
            if (shrate_eff > ecmech::idp_tiny_sqrt) {
               sv_sat = _taus0 * pow((shrate_eff / _gamss0), _xms);
            }
            evolVals[0] = shrate_eff;
            evolVals[1] = sv_sat;
         }

         __ecmech_hdev__
         inline
         void
         getSdotN(double *sdot,
                  double *dsdot_ds,
                  const double* const h,
                  const double* const evolVals) const
         {
            double shrate_eff = evolVals[0];
            double sv_sat = evolVals[1];

            for(int i = 0; i < _nslip * _nslip; i++) {
               dsdot_ds[i] = 0.0;
            }

            for(int iSlip = 0; iSlip < _nslip; iSlip++) {
               //When the below ternary op is true then sdot and dsdot_ds remain zero.
               double temp2 = (sv_sat <= _tausi[iSlip]) ? zero : one / (sv_sat - _tausi[iSlip]);

               // IF (PRESENT(dfdtK)) THEN
               // dfdtK(1) = zero
               // END IF
               // Just throwing this in here for the tests
               // in reality we would need to set dsdot_ds in another section
               // after checking if it's a nullptr or not
               assert(dsdot_ds != nullptr);

               if (nonlinear) {
                  double temp1 = pow((sv_sat - h[iSlip]) * temp2, _xmprime1);
                  sdot[iSlip] = _h0 * temp1 * (sv_sat - h[iSlip]) * temp2 * shrate_eff;
                  dsdot_ds[iSlip * _nslip + iSlip] = -_h0 * temp2 * shrate_eff * _xmprime * temp1;
               }
               else {
                  double temp1 = _h0 * ((sv_sat - h[iSlip]) * temp2);
                  sdot[iSlip] = temp1 * shrate_eff;
                  // double dfdshr = temp1 + _h0 * ( (h - _tausi) / (temp2*temp2)) * _xms * sv_sat ;
                  dsdot_ds[iSlip * _nslip + iSlip] = -_h0 * temp2 * shrate_eff;
               }
            }
         }
   }; // class KineticsVocePL
} // namespace ecmech

#endif // ECMECH_KINETICS_VOCEPL_H
