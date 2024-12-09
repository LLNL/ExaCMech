// -*-c++-*-

#ifndef ECMECH_KINETICS_BCCMD_H
#define ECMECH_KINETICS_BCCMD_H

#include <cassert>
#include <cmath>

#define ECMECH_NN_INDX(p, q, nDim) (p) * (nDim) + (q)

namespace ecmech {
   /**
    * slip and hardening kinetics
    *
    * Template on the slip geometry class
    */
   template<class SlipGeom>
   class KineticsBCCMD
   {
      public:
         /// Number of hardening state variables
         /// The hardening state can be either the CRSS or it could be something like
         /// the DD content or someting else
         static constexpr int nH = SlipGeom::nslip;
         /// Number of slip systems we're dealing with if it that is something useful
         static constexpr int m_num_slip = SlipGeom::nslip;
         /// Number of parameters the model needs to be instantiated
         static constexpr int nParams = 8+4+1+3;
         /// Number of slip kinetic related-variables outputted
         static constexpr int nVals = 2 * SlipGeom::nslip + 1;
         /// These are variables that the hardening equation would need to solve for
         /// its update but the variables are not constant themselves.
         static constexpr int nEvolVals = nH + 1;

         // Generally  don't using anything other than the default here
         __ecmech_hdev__
         KineticsBCCMD(int) {}
         // deconstructor
         ~KineticsBCCMD() = default;

         __ecmech_hdev__
         KineticsBCCMD(const double* const params, int)
         {
            setParams(params);
         }

         __ecmech_host__
         inline void setParams(const std::vector<double> & params)
         {
            setParams(params.data());
         }

         __ecmech_hdev__
         inline
         void setParams(const double* const params) {
            const double* parsIt = params;

            //////////////////////////////
            // power-law stuff
            // shear modulus in case the model uses it
            m_mu = *parsIt; ++parsIt;
            // Burgers vector magnitude
            m_bmag = *parsIt; ++parsIt;
            // This would be the power law exponent term
            m_xm = *parsIt; ++parsIt;
            // This would be the references slip rate term
            m_gam_w0 = *parsIt; ++parsIt;
            // Peierls stress
            m_tau_p = *parsIt; ++parsIt;
            // alpha Peierls
            m_alpha_p = *parsIt; ++parsIt;
            // Shear velocity
            m_vmax = *parsIt; ++parsIt;
            // Drag stress
            m_tau_drag = *parsIt; ++parsIt;

            // These are terms that are constant during the simulation and we don't
            // really need to calculate them every time we call slip kinetics portion
            // of the class
            m_xnn = one / m_xm;
            m_xn = m_xnn - one;
            //
            // CALL set_t_min_max(pl)
            // For numerics, we define a minimum and maximum (rss / crss) value
            // that translates to either a slip rate that is essentially zero
            // or slip rate that is going off to infinity but we really want to
            // cap it to some large number
            m_t_min = pow(ecmech::gam_ratio_min, m_xm);
            m_t_max = pow(ecmech::gam_ratio_ovf, m_xm);

            //////////////////////////////
            // Hardening parameters
            m_alpha = *parsIt; ++parsIt;
            m_k1 = *parsIt; ++parsIt;
            m_k2 = *parsIt; ++parsIt;
            m_krelax = *parsIt; ++parsIt;
            m_gdot_0 = *parsIt; ++parsIt;
            m_tkelv0 = *parsIt; ++parsIt;
            m_ak = *parsIt; ++parsIt;
            //////////////////////////////
            // nH
            // All the terms related to our hardening state
            // You'll often see the initial state provided to us saved off as well.
            // This is just so when you call getParam the initial state can be provided back
            // Also, if you'd like to based on the initial state you could provide a lower
            // bound for example related to the dislocation content that you don't want the
            // model to go under when you start try to update your hardening state.
            m_hdn_init = *parsIt; ++parsIt;
            m_hdn_min = 1e-4 * m_hdn_init;

            //////////////////////////////

#if defined(ECMECH_DEBUG)
            int iParam = parsIt - params;
            if (iParam != nParams) {
               ECMECH_FAIL(__func__, "iParam != nParams");
            }
#endif
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

            params.push_back(m_mu);
            params.push_back(m_bmag);
            params.push_back(m_xm);
            params.push_back(m_gam_w0);
            params.push_back(m_tau_p);
            params.push_back(m_alpha_p);
            params.push_back(m_vmax);
            params.push_back(m_tau_drag);

            //////////////////////////////
            // hardening stuff
            params.push_back(m_alpha);
            params.push_back(m_k1);
            params.push_back(m_k2);
            params.push_back(m_krelax);
            params.push_back(m_gdot_0);
            params.push_back(m_tkelv0);
            params.push_back(m_ak);
            //////////////////////////////
            // nH

            params.push_back(m_hdn_init);

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
               init.push_back(m_hdn_init);
               plot.push_back(true);
               state.push_back(true);
            }
         }

      private:

         //////////////////////////////
         // Power-law stuff

         double m_mu, m_bmag;
         double m_xm;
         double m_gam_w0;
         double m_tau_p, m_alpha_p;
         double m_vmax, m_tau_drag;

         // derived from parameters
         double m_t_max, m_t_min, m_xn, m_xnn;

         //////////////////////////////
         // Hardening
         double m_alpha;
         double m_k1, m_k2, m_krelax;
         double m_gdot_0, m_tkelv0, m_ak;
         //////////////////////////////
         // nH
         double m_hdn_init;
         double m_hdn_min;

      public:

         /// This is used to help scale portions of the elastic strain and lattice
         /// rotation solve. Traditionally, it contains the reference slip rate
         /// For models where that value is slip system dependent,
         /// you can do something similar to what I did with the orowan model
         /*      // Thermal activation + phonon ref slip rate = (1/(f_D * \bar{L}/b * sqrt(qM_0)/sqrt(qM)) + 1/(gammadot_r0 * qM))^-1
               const double isqrth = 1.0 / sqrt(vals[1 + m_num_slip + iVal]);
               const double rate = 1.0 / ((1.0 / (_lbar_b * _fD * isqrth)) + (1.0 / (_gam_ro * vals[1 + m_num_slip + iVal])));
               if (rate > maxRefRate) {
                  maxRefRate = rate;
               }
         */
         __ecmech_hdev__
         inline double getFixedRefRate(const double* const // vals, not used
                                       ) const
         {
            return m_gam_w0;
         }

         /// This is where all of those kinetic values are evaluated
         /// vals are the kinetic values - which can contain things like the
         /// reference slip rates, CRSS values, or a constant term that is divided by
         /// temperature
         /// p down below is the pressure term and tkelv is the temperature
         /// h_state is the hardness state (CRSS for voce model and DD content for orowan model)
         /// Also, it returns the average flow strength (CRSS value) across all slip systems
         __ecmech_hdev__
         inline
         double
         getVals(double* const vals,
                 double, // p, not currently used
                 double tkelv,
                 const double* const h_state
                 ) const
         {
            double crss = ecmech::zero;
            for (int iSlip = 0; iSlip < m_num_slip; ++iSlip) {
               crss += h_state[iSlip];
               assert(h_state[iSlip] > zero);
            }
            crss = m_alpha * m_mu * m_bmag * sqrt(crss);
            
            double mVals = ecmech::zero;
            for (int iSlip = 0; iSlip < m_num_slip; ++iSlip) {
               vals[iSlip] = crss;
               vals[m_num_slip + iSlip] = h_state[iSlip];
               mVals += vals[iSlip];
               assert(vals[iSlip] > zero);
            }
            constexpr double inv_nslip = 1.0 / m_num_slip;
            mVals *= inv_nslip;
            
            vals[2*m_num_slip] = tkelv;

            return mVals;
         }

         /// Evaluates our slip rate and its derivatives when provided the RSS value across all slip systems
         /// and the kinetic values calculated in getVals
         /// The derivative we need is the derivative of the slip rate wrt the RSS         
         __ecmech_hdev__
         inline
         void
         evalGdots(double* const gdot,
                   double* const dgdot_dtau,
                   const double* const tau,
                   const double* const vals
                   ) const
         {     
            double tkelv = 0.0;

            for (int iSlip = 0; iSlip < m_num_slip; ++iSlip) {
               bool l_act;
               double taua = tau[iSlip];
               double chia = tau[SlipGeom::nslip + iSlip];
               
               double crss = vals[iSlip];
               double rhoa = vals[SlipGeom::nslip + iSlip];
               // traditionally we have a separate function that will calculate everything
               // for only one slip system
               this->evalGdot(gdot[iSlip], l_act, dgdot_dtau[iSlip],
                              crss, rhoa, taua, chia, tkelv);
            }
         }

         /// Calculates the slip rate and derivatives for a given slip system
         __ecmech_hdev__
         inline
         void
         evalGdot(
            double & gdot,
            bool  & l_act,
            double & dgdot_dtau, // wrt resolved shear stress
            double   crss,
            double   rho,
            double   tau,
            double   chi,
            double   /*tkelv*/
            ) const
         {
            // zero things so that can more easily just return in inactive
            //// gdot_w = zero; gdot_r = zero; ! not used by l_linear or l_pl
            gdot = zero;
            //
            dgdot_dtau = zero;            
            l_act = false;

            double tau_p = m_tau_p / cos(chi-m_alpha_p);

            double t_eff = fmax(fabs(tau) - tau_p, 0.0);

            double xnn = m_xnn;
            double xn = m_xn;
            double gam_w0 = m_gam_w0;

            double v0_T = m_gam_w0;
            double v0_AT = v0_T * 0.1;
            gam_w0 = v0_T + (M_PI/6.0 + chi) * 3.0/M_PI * (v0_AT - v0_T);

            double g_i = one / crss; // assume have checked gIn>0 elsewhere
            double t_frac = t_eff * g_i;
            t_frac = copysign(t_frac, tau); // has sign of tau
            double at = fabs(t_frac);
            
            double gam_w;
            if (gam_w0 < 0.0) gam_w = fabs(gam_w0);
            else gam_w = rho * m_bmag * gam_w0;
            
            double gmax = rho * m_bmag * m_vmax * (1.0-exp(-t_eff/m_tau_drag));

            if (at > m_t_min) {
               l_act = true;
               if (at > m_t_max) {
                  // ierr = IERR_OVF_p
                  // set gdot big, evpp may need this for recovery
                  gdot = ecmech::gam_ratio_ovffx * gam_w;
                  gdot = copysign(gdot, tau);
                  // do not set any of deriviatives (they are, in truth, zero)
               }
               else {
                  double abslog = log(at);
                  double blog = xn * abslog;
                  double temp = gam_w * exp(blog);

                  gdot = temp * t_frac;

                  dgdot_dtau = temp * xnn * g_i; // note: always positive, = xnn * gdot/t
                  // Smooth capping to gmax with Lorentz-like factor
                   double ac = 10.0;
                   double gfrac_a = pow(fabs(gdot)/gmax, ac);
                   double dfact = pow(1.0 + gfrac_a, -(ac+1.0)/ac);
                   
                   // Correct derivative accounting for gmax(tau)
                   double dgmax = rho * m_bmag * m_vmax / m_tau_drag * exp(-t_eff/m_tau_drag);
                   double A = fabs(gdot) * dgmax * gfrac_a;
                   dgdot_dtau = dfact * (A + gmax * dgdot_dtau) / gmax;
                   gdot = gdot / pow(1.0 + gfrac_a, 1.0/ac);
               }
            }
         } // evalGdot

         /// This is called externally by the portion of code that does the
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
                 const double* const hvals,
                 double tkelv,
                 int outputLevel = 0) const
         {
            double log_hs_u[SlipGeom::nslip];
            double log_hs_o[SlipGeom::nslip];
            double evolVals[nEvolVals] = {};
            
            for(int islip = 0; islip < SlipGeom::nslip; islip++) {
               log_hs_o[islip] = log(fmax(hs_o[islip], m_hdn_min));
            }
            getEvolVals(evolVals, gdot);
            // If the equation is incredibly  stiff it's possible this won't solve
            int nFEvals = updateHN<KineticsBCCMD, true>(this,
                                                  log_hs_u, log_hs_o, dt, evolVals, hvals, tkelv,
                                                  outputLevel);

            for(int islip = 0; islip < SlipGeom::nslip; islip++) {
               hs_u[islip] = exp(log_hs_u[islip]);
            }
            return nFEvals;
         }

         /// This calculates the variables I'd mentioned up above and now again down below
         /// related to the hardening state
         /// These are variables that the hardening equation would need to solve for
         /// its update but the variables are not constant themselves.
         /// A common set would be for example in a voce model, the updated
         /// saturation strength (g^{sat}_0 (\frac{\sum_{i = 0}^{number of slip systems} |\dot{\gamma}_i| }{constant})^m')
         /// as the saturation strength evolves based on the sum of the absolute value of the gammadots.
         __ecmech_hdev__
         inline
         void
         getEvolVals(double* const evolVals,
                     const double* const gdot
                     ) const
         {
            for (int i = 0; i < m_num_slip; i++) {
                evolVals[i] = abs(gdot[i]);
                evolVals[nEvolVals - 1] += abs(gdot[i]);
            }
         }

         /// This function calculates time rate of change of the hardening state (sdot) 
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
                  const double* const h_i,
                  const double* const evolVals,
                  const double* const hvals,
                  double tkelv
                ) const
         {
            constexpr bool LOGFORM = true;
            constexpr size_t nslip = SlipGeom::nslip;
            constexpr size_t JDIM = 2;
            constexpr size_t nDimSys = SlipGeom::nslip;
            constexpr size_t h_content = (LOGFORM) ? nslip : 1;

            double hexp[h_content];
            if (LOGFORM) {
               for (size_t iDD = 0; iDD < nslip; iDD++) {
                  // Prevent divide by 0 errors...
                  hexp[iDD] = fmax(exp(h_i[iDD]), ecmech::idp_eps);
               }
            }
            const double* const h = (LOGFORM) ? &hexp[0] : h_i;

            const double gamma = evolVals[nEvolVals - 1];
            auto k1_func = [=](const double xi) -> double {
               /*
               Current version does things a bit different from the paper.
               double chia = hvals[islip];
               // If we are in the AT zone, then we need to increase k1
               // to account for the fact that dislocations do take
               // a longer path and thus are likely to multiply more
               double amin = 0.95;
         
               //double a = fmin(1.0 + (amin - 1.0) * xi * 6.0 / M_PI, 1.0);
               //a = 1.0 / a;
         
               double a = 1.0/(1.0/cos(M_PI/6.0 - m_alpha_p)-1.0) * (1.0/amin - 1.0);
               a = 1.0 + a * (1.0 / cos(xi - m_alpha_p) - 1.0);
         
               return m_k1 * a;
               */
               // Pure paper implementation of things
               return m_k1 * (1.0 +(m_ak / (cos(xi - m_alpha_p))));
            };

            auto k2_func = [=] () -> double {
               const double gamma_ratio = (gamma > ecmech::gam_ratio_min) ? (gamma / m_gdot_0) : ecmech::gam_ratio_min;
               return m_k2 * log(gamma_ratio) * log(tkelv / m_tkelv0);
            };

            auto f_func = [=] (const double abs_gamma_dot) -> double {
               constexpr double A = 100.0;
               constexpr double t = 0.01;
               const double gamma_ratio = (gamma > ecmech::gam_ratio_min) ? (abs_gamma_dot / gamma) : gam_ratio_max;
               const double exp_inner = -A * (gamma_ratio - t);
               return 1.0 - (1.0 / ( 1.0 + exp(exp_inner)));
            };

            auto k_relax_func = [=] (const double h) -> double {
               const double relax_term = 1.0 - exp(-(h - m_hdn_min) / m_hdn_min);
               return relax_term;
            };

            // From the paper if sqrt(a_ij * rho_j) does not correspond to
            // sqrt(I_ij * rho_j) where I is the identity matrix
            // then we'd something like the below
            double m_a_mat[nslip * nslip] = {};
            {
               RAJA::View<double, RAJA::Layout<JDIM> > amat(&m_a_mat[0], nslip, nslip);
               for (size_t islip = 0; islip < nslip; islip++) {
                  amat(islip, islip) = 1.0;
               }
            }
            double amat_rho[SlipGeom::nslip] = {};
            vecsVMa<SlipGeom::nslip>(&amat_rho[0], &m_a_mat[0], &h[0]);

            for (int islip = 0; islip < SlipGeom::nslip; islip++) {
               const double k1   = k1_func(hvals[islip]) * evolVals[islip];
               const double k2   = k2_func() * evolVals[islip];
               const double fval = f_func(evolVals[islip]) * m_krelax * k_relax_func(h[islip]);
               amat_rho[islip] = sqrt(amat_rho[islip]);
               sdot[islip] = (k1 * amat_rho[islip] - k2 * h[islip]) - fval * h[islip];
            }

            if (LOGFORM) {
               for (int iDD = 0; iDD < SlipGeom::nslip; iDD++) {
                  sdot[iDD] *= (1.0 / h[iDD]);
               }
            }
            if (dsdot_ds)
            {
               // zero out dsdot_ds matrix
               for (size_t i = 0; i < nDimSys * nDimSys; i++) {
                  dsdot_ds[i] = ecmech::zero;
               }
               RAJA::View<const double, RAJA::Layout<JDIM> > amat(&m_a_mat[0], nslip, nslip);
               RAJA::View<double, RAJA::Layout<JDIM> > dsdot_ds_view(dsdot_ds, nDimSys, nDimSys);
               for (size_t islip = 0; islip < nslip; islip++) {
                  const double k1   = k1_func(hvals[islip]) * evolVals[islip];
                  const double k2   = k2_func() * evolVals[islip];
                  const double fval = f_func(evolVals[islip]) * m_krelax;
                  const double ratio = h[islip] / m_hdn_min;
                  // Don't want this getting too big or else we get an infinity error later on...
                  const double ratio_max = fmin(ratio, 80.0);
                  // Copied this from Wolfram alpha for the f_func() * k_relax * k_relax_func() * h
                  const double fval_der = fval - (fval * exp(1 - ratio_max) * (m_hdn_min - h[islip]))/m_hdn_min;
                  dsdot_ds_view(islip, islip) += -(k2 + fval_der);
                  for (size_t jslip = 0; jslip < nslip; jslip++) {
                     dsdot_ds_view(islip, jslip) += k1 * amat(islip, jslip) * 0.5 / amat_rho[jslip];
                  }
               }
               // Generic solution to transform over into log space
               if (LOGFORM) {
                  for (size_t iDD = 0; iDD < nslip; iDD++) {
                     dsdot_ds_view(iDD, iDD) -= sdot[iDD];
                  }
                  for (size_t iDD = 0; iDD < nslip; iDD++) {
                     for (size_t jDD = 0; jDD < nslip; jDD++) {
                        dsdot_ds_view(iDD, jDD) *= h[jDD] / h[iDD];
                     }
                  }
               }
            }
         }
   }; // class KineticsBCCMD
} // namespace ecmech

#endif // ECMECH_KINETICS_BCCMD_H
