// -*-c++-*-

#ifndef ECMECH_KINETICS_BCCMD_H
#define ECMECH_KINETICS_BCCMD_H

#include <cassert>
#include <cmath>

#define ECMECH_NN_INDX(p, q, nDim) (p) * (nDim) + (q)

// I've left out some of the additional functions in this class which are related to
// extra derivatives that are needed for fully implicit solve of the updated state (elastic strain, lattice rotation, and hardening state)
// as calculating those extra derivatives usually aren't worth the hassle if we are just playing
// around with models.
// Next for the slip geometry class I will need to make some changes over there
// so the new slip geometries can just plug and play with the existing solvers.

namespace ecmech {
   /**
    * slip and hardening kinetics
    *
    * power-law slip kinetics with some hardening law -- used as a template
    * Template on the slip geometry class
    */
   template<class SlipGeom>
   class KineticsBCCMD
   {
      public:
         /// Number of hardening state variables
         /// The hardening state can be either the CRSS or it could be something like
         /// the DD content or someting else
         static const int nH = SlipGeom::nslip;
         /// Number of slip systems we're dealing with if it that is something useful
         static const int m_num_slip = SlipGeom::nslip;
         /// Number of parameters the model needs to be instantiated
         static const int nParams = 4+3+1;
         /// Number of slip kinetic related-variables outputted
         /// Think of this as things like the CRSS values, evolving reference
         /// slip rates for both thermal and phonon drag contributions, and potentially
         /// other evolving variables that we can calculate at the beginning of time
         /// step and not have to recalculate every iterations of our coupled solve
         /// of the elastic strain and lattice rotation
         static const int nVals = SlipGeom::nslip;
         static const int nValsDerivs = SlipGeom::nslip;
         /// These are variables that the hardening equation would need to solve for
         /// its update but the variables are not constant themselves.
         /// A common set would be for example in a voce model, the updated
         /// saturation strength (g^{sat}_0 (\frac{\sum_{i = 0}^{number of slip systems} |\dot{\gamma}_i| }{constant})^m')
         /// as the saturation strength evolves based on the sum of the absolute value of the gammadots.
         /// In the orowan model as another example, we need the signed mobile dislocation scalar velocity
         /// as an input.
         static const int nEvolVals = nH;

         // Generally  don't using anything other than the default here
         __ecmech_hdev__
         KineticsBCCMD(int nslip) : _nslip(nslip) {};
         // deconstructor
         __ecmech_hdev__
         ~KineticsBCCMD() {}

         /// In ExaCMech each class will be handed the parameters that they said they needed
         /// It is up to the modeller to iterate through this vector and  pull out the parameters
         /// and put them where they need to go.
         /// Additionally, modellers could also generate other parameters based on the inputted ones
         /// that the model will use later on.
         __ecmech_host__
         inline void setParams(const std::vector<double> & params // const double* const params
                               ) {
            std::vector<double>::const_iterator parsIt = params.begin();

            //////////////////////////////
            // power-law stuff
            // shear modulus in case the model uses it
            _mu = *parsIt; ++parsIt;
            // Burgers vector magnitude
            _bmag = *parsIt; ++parsIt;
            // This would be the power law exponent term
            _xm = *parsIt; ++parsIt;
            // This would be the references slip rate term
            _gam_w = *parsIt; ++parsIt;

            // CALL fill_power_law(pl)
            // xmm  = xm - one ;
            // These are terms that are constant during the simulation and we don't
            // really need to calculate them every time we call slip kinetics portion
            // of the class
            _xnn = one / _xm;
            _xn = _xnn - one;
            // xMp1 = xnn + one
            //
            // CALL set_t_min_max(pl)
            // For numerics, we define a minimum and maximum (rss / crss) value
            // that translates to either a slip rate that is essentially zero
            // or slip rate that is going off to infinity but we really want to
            // cap it to some large number
            _t_min = pow(ecmech::gam_ratio_min, _xm);
            _t_max = pow(ecmech::gam_ratio_ovf, _xm);

            //////////////////////////////
            // Hardening parameters
            _alpha = *parsIt; ++parsIt;
            _k1 = *parsIt; ++parsIt;
            _k2 = *parsIt; ++parsIt;
            
            //////////////////////////////
            // nH
            // All the terms related to our hardening state
            // You'll often see the initial state provided to us saved off as well.
            // This is just so when you call getParam the initial state can be provided back
            // Also, if you'd like to based on the initial state you could provide a lower
            // bound for example related to the dislocation content that you don't want the
            // model to go under when you start try to update your hardening state.
            _hdn_init = *parsIt; ++parsIt;
            _hdn_min = 1e-4 * _hdn_init;

            //////////////////////////////

            assert((parsIt - params.begin()) == nParams);
         }

         /// Here you'll just return all the parameters that were provided to you
         /// up above in the same order as up above as well.
         __ecmech_host__
         inline void getParams(std::vector<double> & params
                               ) const {
#ifdef ECMECH_DEBUG
            // do not clear params in case adding to an existing set
            int paramsStart = params.size();
#endif

            //////////////////////////////
            // power-law stuff

            params.push_back(_mu);
            params.push_back(_bmag);
            params.push_back(_xm);
            params.push_back(_gam_w);

            //////////////////////////////
            // hardening stuff
            params.push_back(_alpha);
            params.push_back(_k1);
            params.push_back(_k2);

            //////////////////////////////
            // nH

            params.push_back(_hdn_init);

            //////////////////////////////
#ifdef ECMECH_DEBUG
            assert((params.size() - paramsStart) == nParams);
#endif
         }

         /// Return the nH related variables along with any names you might want
         /// associated with them
         __ecmech_host__
         void getHistInfo(std::vector<std::string> & names,
                          std::vector<double>       & init,
                          std::vector<bool>        & plot,
                          std::vector<bool>        & state) const {
            
            for (int iSlip = 0; iSlip < SlipGeom::nslip; iSlip++) {
               names.push_back("rho_" + std::to_string(iSlip));
               //init.push_back(h_state[iSlip]);
               init.push_back(_hdn_init);
               plot.push_back(true);
               state.push_back(true);
            }
         }

      private:

         // static const _nXnDim = nH*nH ; // do not bother
         const int _nslip; // could template on this if there were call to do so

         //////////////////////////////
         // Power-law stuff

         double _mu, _bmag;
         double _xm;
         double _gam_w;

         // derived from parameters
         double _t_max, _t_min, _xn, _xnn;

         //////////////////////////////
         // Hardening
         double _alpha;
         double _k1, _k2;

         //////////////////////////////
         // nH
         double _hdn_init;
         double _hdn_min;

      public:

         /// This is used to help scale portions of the elastic strain and lattice
         /// rotation solve. Traditionally, it contains the reference slip rate
         /// For models where that value is slip system dependent,
         /// you can do something similar to what I did with the orowan model
         /*      // Thermal activation + phonon ref slip rate = (1/(f_D * \bar{L}/b * sqrt(qM_0)/sqrt(qM)) + 1/(gammadot_r0 * qM))^-1
               const double isqrth = 1.0 / sqrt(vals[1 + _nslip + iVal]);
               const double rate = 1.0 / ((1.0 / (_lbar_b * _fD * isqrth)) + (1.0 / (_gam_ro * vals[1 + _nslip + iVal])));
               if (rate > maxRefRate) {
                  maxRefRate = rate;
               }
         */
         __ecmech_hdev__
         inline double getFixedRefRate(const double* const // vals, not used
                                       ) const
         {
            return _gam_w;
         }

         /// This is where all of those kinetic values are evaluated
         /// vals are the kinetic values - which can contain things like the
         /// reference slip rates, CRSS values, or a constant term that is divided by
         /// temperature
         /// p down below is the pressure term and tK is the temperature
         /// h_state is the hardness state (CRSS for voce model and DD content for orowan model)
         /// Also, it returns the average flow strength (CRSS value) across all slip systems
         __ecmech_hdev__
         inline
         double
         getVals(double* const vals,
                 double, // p, not currently used
                 double, // tK, not currently used
                 const double* const h_state,
                 double* const val_derivs = nullptr
                 ) const
         {
            assert(val_derivs == nullptr);
            double crss = ecmech::zero;
            for (int iSlip = 0; iSlip < _nslip; ++iSlip) {
               crss += h_state[iSlip];
               assert(h_state[iSlip] > zero);
            }
            crss = _alpha*_mu*_bmag*sqrt(crss);
            
            double mVals = ecmech::zero;
            for (int iSlip = 0; iSlip < _nslip; ++iSlip) {
               vals[iSlip] = crss;
               mVals += vals[iSlip];
               assert(vals[iSlip] > zero);
            }
            mVals /= _nslip;

            return mVals;
         }

         /// Evaluates our slip rate and its derivatives when provided the RSS value across all slip systems
         /// and the kinetic values calculated in getVals
         /// The derivatives we need are the derivative of the slip rate wrt the RSS and
         /// the derivative of the slip rate wrt the CRSS
         __ecmech_hdev__
         inline
         void
         evalGdots(double* const gdot,
                   double* const dgdot_dtau,
                   double* const dgdot_dg,
                   const double* const tau,
                   const double* const vals,
                   const bool dgdot_dh_conv = false,
                   const double* const val_drivs = nullptr
                   ) const
         {
            assert(dgdot_dh_conv == false);
            assert(val_drivs == nullptr);

            for (int iSlip = 0; iSlip < _nslip; ++iSlip) {
               bool l_act;
               double taua = tau[iSlip];
               double chia = tau[SlipGeom::nslip+iSlip];
               
               //printf("sys[%d] tau = %e, chi = %e\n",iSlip,taua,chia*180.0/M_PI);
               
               double gAll = vals[iSlip];
               // traditionally we have a separate function that will calculate everything
               // for only one slip system
               this->evalGdot(gdot[iSlip], l_act, dgdot_dtau[iSlip], dgdot_dg[iSlip],
                              gAll, taua, _mu);
            }
         }

         /// Calculates the slip rate and derivatives for a given slip system
         /// The MORE_DERIVS portion of things isn't used at this point by ECMech
         /// so we can probably just set them to 0 within another ifdef down below
         /// or just ignore them completely.
         /// l_act can just be ignored we don't actually use it.
         __ecmech_hdev__
         inline
         void
         evalGdot(
            double & gdot,
            bool  & l_act,
            double & dgdot_dtau, // wrt resolved shear stress
            double & dgdot_dg, // wrt slip system strength
            double   gIn,
            double   tau,
            double // mu, not currently used
            ) const
         {
            // zero things so that can more easily just return in inactive
            //// gdot_w = zero; gdot_r = zero; ! not used by l_linear or l_pl
            gdot = zero;
            //
            dgdot_dtau = zero;
            dgdot_dg = zero;
            
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
               }
            }
         } // evalGdot

         /// This is called externally  by the portion of code that does the
         /// elastic strain and lattice rotation update. However, it's only
         /// called at the beginning of time step and is not called iteratively
         /// so the inputs are all begining of time step values
         /// hs_u is our updated hardening state variable
         /// hs_o is our beginning of time step hardening state variable
         /// dt is our delta time step
         /// gdot is our beginning of time step slip rate
         /// output level is optional but it can be used for logging purposes
         __ecmech_hdev__
         inline
         int
         updateH(double* const hs_u,
                 const double* const hs_o,
                 double dt,
                 const double* const gdot,
                 int outputLevel = 0) const
         {
            double log_hs_u[SlipGeom::nslip];
            double log_hs_o[SlipGeom::nslip];
            double gdotabs[SlipGeom::nslip];
            
            for(int islip = 0; islip < SlipGeom::nslip; islip++) {
               log_hs_o[islip] = log(fmax(hs_o[0], _hdn_min));
               gdotabs[islip] = abs(gdot[islip]);
            }

            // If the equation is incredibly  stiff it's possible this won't solve
            int nFEvals = updateHN<KineticsBCCMD>(this,
                                                  log_hs_u, log_hs_o, dt, gdotabs,
                                                  outputLevel);

            for(int islip = 0; islip < SlipGeom::nslip; islip++) {
               hs_u[islip] = exp(log_hs_u[islip]);;
            }

            return nFEvals;
         }
         
         __ecmech_hdev__
         inline
         void
         setH0Ext(double *const h0) const
         {
            for (int i = 0; i < nH; i++) {
               h0[i] = log(fmax(h0[i], _hdn_min));
            }
            return;
         }
         
         __ecmech_hdev__
         inline
         void
         getHUpdate(const double *const h0,
                    const double *const del_h,
                    const double *const del_h_scale,
                    double *const       h,
                    const bool /*updateFinal*/) const
         {
            // We always return the non-log form of h even though
            // we get the log form in as we need to make use of the
            // regular form within the kinetics update and gdot eval
            // calculations
            for (int i = 0; i < nH; i++) {
               const double factor = h0[i] + del_h[i] * del_h_scale[i];
               h[i] = exp(factor);
            }
         }
         
         __ecmech_hdev__
         inline
         void
         getExtDerivs(double* const hdot,
                      double* const dhdot_dh,
                      double* const dhdot_dgdot,
                      double* const /*dgdot_dh*/,
                      double* const hard,
                      const double* const gdot) const
         {
            double gdotabs[SlipGeom::nslip];
            double evolVals[nEvolVals];
            // Transform this back into the log form for the later residual calculation
            for(int islip = 0; islip < SlipGeom::nslip; islip++) {
               hard[islip] = log(hard[islip]);
               gdotabs[islip] = abs(gdot[islip]);
            }
            getEvolVals(evolVals, gdotabs);
            getSdotN(hdot, dhdot_dh, hard, evolVals, dhdot_dgdot);
         }

         /// This calculates the variables I'd mentioned up above and now again down below
         /// related to the hardening state
         /// These are variables that the hardening equation would need to solve for
         /// its update but the variables are not constant themselves.
         /// A common set would be for example in a voce model, the updated
         /// saturation strength (g^{sat}_0 (\frac{\sum_{i = 0}^{number of slip systems} |\dot{\gamma}_i| }{constant})^m')
         /// as the saturation strength evolves based on the sum of the absolute value of the gammadots.
         /// In the orowan model as another example, we need the signed mobile dislocation scalar velocity
         /// as an input. 
         __ecmech_hdev__
         inline
         void
         getEvolVals(double* const evolVals,
                     const double* const gdotabs
                     ) const
         {
            for (int i = 0; i < _nslip; i++) {
                evolVals[i] = gdotabs[i];
            }
         }

         /// This function does not have a great name.
         /// It calculates time rate of change of the hardening state (sdot) 
         /// and its derivatives dsdot_ds which is the derivative of the
         /// time rate of change of the hardening state wrt the hardening state
         /// Input values are h - hardening state
         /// evolVals - which are values calculated from getEvolVals
         /// This function is called from updateHN<KineticsBCCMD>
         /// and used as part of the nonlinear solve for the updated state
         __ecmech_hdev__
         inline
         void
         getSdotN(double *sdot,
                  double *dsdot_ds,
                  const double* const h,
                  const double* const evolVals,
                  double* const dsdot_dgdot = nullptr // optional parameter
                ) const
         {
            assert(dsdot_dgdot == nullptr);
            {
               // we normally just assume  this value always exists
               const int nslip2 = SlipGeom::nslip * SlipGeom::nslip;
               for (int i = 0; i < nslip2; i++) {
                  dsdot_ds[i] = ecmech::zero;
               }
            }
            // h = log(DD)
            // dDD / dt = DD * dh / dt
            // dh / dt = dDD / dt * 1 / DD
            // d DD_i / dt = (k1 * sqrt(A_{ij} DD_j) - k2 * DD_i) * gammadot_i
            // dh / dt = (k1 * sqrt(A_{ij} DD_j) / DD_i - k2) * gammadot_i
            // specialized here for the A_{ij} = I
            // dh / dt = (k1 / sqrt(DD_i) - k2) * gammadot_i
            // specialized case
            // \dot{h} / dh = -1/2 * k_1 * (DD)^{-1/2}
            // more general case I believe if I did the derivs correctly...
            // \dot{h^i} / dh_j = \dot{h^i} / d DD_j * d DD^j / d h_k
            // d DD^j / d h_k = DD_j when j == k and 0 for j neq k
            // for i neq j
            // 1/2 *  \frac{k_1 * A_{ij}}{DD_i * \sqrt(A_{ij}DD_j)} * gammadot_i * DD_j
            // for i == j
            // (1/2 *  \frac{k_1 * A_{ij}}{DD_i * \sqrt(A_{ij}DD_j)} - \frac{k1 * \sqrt(A_{ij}DD_j)}{DD_i^2} ) * gammadot_i * DD_j
            // = (\frac{k1 A_{ij} DD_i - 2 k1 * A_ij DD_j}{2 * DD^2_i * sqrt(A_{ij} * DD_j)}) gammadot_i * DD_j
            // when A_ij = I this reduces down to
            // -k1 / 2 * (DD_i)^{-1/2} * gammadot_i 
            // which is what we get out in the regular MTS KM model so that's a good sign
            // I did something right and the off diagonal terms are zero
            for (int islip = 0; islip < SlipGeom::nslip; islip++) {
               double temp_hs_a = exp(-onehalf * h[islip]);
               double temp1 = _k1 * temp_hs_a - _k2;
               sdot[islip] = temp1 * evolVals[islip];
               dsdot_ds[ECMECH_NN_INDX(islip, islip, SlipGeom::nslip)] = (-_k1 * onehalf * temp_hs_a) * evolVals[islip];
            }
         }
   }; // class KineticsBCCMD
} // namespace ecmech

#endif // ECMECH_KINETICS_BCCMD_H
