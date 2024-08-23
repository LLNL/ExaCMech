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
         static const int nParams = 8+4+1;
         /// Number of slip kinetic related-variables outputted
         /// Think of this as things like the CRSS values, evolving reference
         /// slip rates for both thermal and phonon drag contributions, and potentially
         /// other evolving variables that we can calculate at the beginning of time
         /// step and not have to recalculate every iterations of our coupled solve
         /// of the elastic strain and lattice rotation
         static const int nVals = 2 * SlipGeom::nslip + 1;
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
         KineticsBCCMD(int) {};
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

            // CALL fill_power_law(pl)
            // xmm  = xm - one ;
            // These are terms that are constant during the simulation and we don't
            // really need to calculate them every time we call slip kinetics portion
            // of the class
            m_xnn = one / m_xm;
            m_xn = m_xnn - one;
            // xMp1 = xnn + one
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
               //init.push_back(h_state[iSlip]);
               init.push_back(m_hdn_init);
               plot.push_back(true);
               state.push_back(true);
            }
         }

      private:

         // static const _nXnDim = nH*nH ; // do not bother
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
         /// p down below is the pressure term and temp_k is the temperature
         /// h_state is the hardness state (CRSS for voce model and DD content for orowan model)
         /// Also, it returns the average flow strength (CRSS value) across all slip systems
         __ecmech_hdev__
         inline
         double
         getVals(double* const vals,
                 double, // p, not currently used
                 double temp_k,
                 const double* const h_state,
                 double* const val_derivs = nullptr
                 ) const
         {
            assert(val_derivs == nullptr);
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
            mVals /= m_num_slip;
            
            vals[2*m_num_slip] = temp_k;

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
                   const double* const tau,
                   const double* const vals
                   ) const
         {     
            double temp_k = vals[2 * SlipGeom::nslip];

            for (int iSlip = 0; iSlip < m_num_slip; ++iSlip) {
               bool l_act;
               double taua = tau[iSlip];
               double chia = tau[SlipGeom::nslip + iSlip];
               
               double crss = vals[iSlip];
               double rhoa = vals[SlipGeom::nslip + iSlip];
			   
               // traditionally we have a separate function that will calculate everything
               // for only one slip system
               this->evalGdot(gdot[iSlip], l_act, dgdot_dtau[iSlip],
                              crss, rhoa, taua, chia, temp_k);
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
            double   crss,
            double   rho,
            double   tau,
            double   chi,
            double   /*temp_k*/
            ) const
         {
            // zero things so that can more easily just return in inactive
            //// gdot_w = zero; gdot_r = zero; ! not used by l_linear or l_pl
            gdot = zero;
            //
            dgdot_dtau = zero;            
            l_act = false;
        
            //double u = cos(chi + M_PI/6.0); // 1 for T, 0.5 for AT
            //double tau_p = m_tau_p + 1.0*(1.0-u) * m_tau_p;
			
            double tau_p = m_tau_p / cos(chi-m_alpha_p);
			
            double t_eff = fmax(fabs(tau) - tau_p, 0.0);
            
            //double u = cos(chi + M_PI/6.0); // 1 for T, 0.5 for AT
            //double t_eff = fmax(fabs(u * tau) - tau_p, 0.0);
			
			
			double xnn = m_xnn;
			double xn = m_xn;
			double gam_w0 = m_gam_w0;
			/*
			// NEW //
			
			double k = 10.0;
			
			double n_exp_T = 5.0;
        	double n_exp_AT = 20.0;
        	//xnn = n_exp_T + (M_PI/6.0+chi)*3.0/M_PI * (n_exp_AT-n_exp_T);
			
			xnn = n_exp_T + 1.0/(1.0+exp(-k*chi)) * (n_exp_AT-n_exp_T);
			xn = xnn - one;
			
			double v0_T = 26.17;
        	double v0_AT = 0.094;
        	//gam_w0 = v0_T + (M_PI/6.0+chi)*3.0/M_PI * (v0_AT-v0_T);
			gam_w0 = v0_T + 1.0/(1.0+exp(-k*chi)) * (v0_AT-v0_T);
			// NEW //
			
			//printf("  xnn = %e, v0 = %e\n", xnn, gam_w0);
			*/
			
#if 1        
            double g_i = one / crss; // assume have checked gIn>0 elsewhere
            double t_frac = t_eff * g_i;
            t_frac = copysign(t_frac, tau); // has sign of tau
            double at = fabs(t_frac);
            
            double gam_w;
            if (gam_w0 < 0.0) gam_w = fabs(gam_w0);
            else gam_w = rho * m_bmag * gam_w0;
            
            double gmax = rho * m_bmag * m_vmax * (1.0-exp(-t_eff/m_tau_drag));
			
			if (at > m_t_min) {
               //
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
                  
                  if (fabs(gdot) > gmax) {
                      gdot = copysign(gmax, tau);
                      dgdot_dtau = zero;
                  }
               }
            }
#else
			{
				double Bdrag = 2e-8; // MBar.us
				double temp = rho * m_bmag * m_bmag / Bdrag;
				
				t_eff = fmax(fabs(tau) - tau_p - crss, 0.0);
				
				gdot = t_eff * temp;
				gdot = copysign(gdot, tau);
				
				dgdot_dtau = temp;
			}
#endif
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
                 const double* const hvals,
                 double temp_k,
                 int outputLevel = 0) const
         {
            double log_hs_u[SlipGeom::nslip];
            double log_hs_o[SlipGeom::nslip];
            double gdotabs[SlipGeom::nslip];
            
            for(int islip = 0; islip < SlipGeom::nslip; islip++) {
               log_hs_o[islip] = log(fmax(hs_o[islip], m_hdn_min));
               gdotabs[islip] = abs(gdot[islip]);
            }
            // If the equation is incredibly  stiff it's possible this won't solve
            int nFEvals = updateHN<KineticsBCCMD>(this,
                                                  log_hs_u, log_hs_o, dt, gdotabs, hvals, temp_k,
                                              outputLevel);

            // We need to check that none of our solutions became negative
            // If we did obtain something negative then we should abort
            // It means our time step was too large for this step.
            // If this is not desirable / possible then we should probably
            // do a terrible hack and cut the dt by some factor resolve things by
            // assuming a constant slip rate during the time step, and then
            // evolve the dd content. We would get a solution, but it wouldn't necessarily
            // be correct.
        #if 1
            bool flag = false;
            for (int islip = 0; islip < SlipGeom::nslip; islip++) {
               if(log_hs_u[islip] < one) {
                  flag = true;
                  break;
               }
            }
            if (flag || nFEvals < 0)
            {
               ECMECH_WARN(__func__, "Solver failed to converge, trying again by substepping through the solution");
               // This is pretty ad-hoc but it seems to work fairly well for a number of simple test cases.
               // It's definitely not the best way to probably do things though...
               int nsub = 10;
               while (nsub < 10000) {
                   const double dtnew = dt / nsub;
                   double log_hs_temp[SlipGeom::nslip];

                   for (int islip = 0; islip < SlipGeom::nslip; islip++) {
                      log_hs_u[islip] = log(fmax(hs_o[islip], m_hdn_min));
                   }

                   for (int i = 0; i < nsub; i++)
                   {
                      for (int islip = 0; islip < SlipGeom::nslip; islip++) {
                         log_hs_temp[islip] = fmax(log_hs_u[islip], log(m_hdn_min));
                      }
                      nFEvals += updateHN<KineticsBCCMD>(this,
                                                         log_hs_u, log_hs_temp, dtnew, gdotabs, hvals, temp_k,
                                                         outputLevel);
                      flag = false;
                      for (int islip = 0; islip < SlipGeom::nslip; islip++) {
                         if(log_hs_u[islip] < one) {
                            flag = true;
                            break;
                         }
                      }
                      
                      if (nFEvals < 0) {
                          flag = true;
                          break;
                      }
                   }
                   
                   if (!flag) break;
                   nsub *= 2;
               }

               if (flag)
               {
                  for (int islip = 0; islip < SlipGeom::nslip; islip++) {
                     printf("dd[%d]: %lf ", islip, exp(log_hs_u[islip]));
                  }
                  printf("\n");
                  ECMECH_FAIL(__func__, "Solver failed to converge!");
               }
            }
        #else
            if (nFEvals < 0) {
                ECMECH_FAIL(__func__, "Solver failed to converge!");
            }
        #endif
            
            
            for(int islip = 0; islip < SlipGeom::nslip; islip++) {
               hs_u[islip] = fmax(exp(log_hs_u[islip]), m_hdn_min);
            }
            //printf("dens = %e %e %e %e\n",hs_u[0]*1e4,hs_u[1]*1e4,hs_u[2]*1e4,hs_u[3]*1e4);

            return nFEvals;
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
            for (int i = 0; i < m_num_slip; i++) {
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
                  const double* const hvals,
                  double temp_k
                ) const
         {
            {
               // we normally just assume  this value always exists
               const int nslip2 = SlipGeom::nslip * SlipGeom::nslip;
               for (int i = 0; i < nslip2; i++) {
                  dsdot_ds[i] = ecmech::zero;
               }
            }
            
            double gtot = 0.0;
            double gdotmax = 0.0;
            double kfact[SlipGeom::nslip];
            for (int islip = 0; islip < SlipGeom::nslip; islip++) {
                gtot += evolVals[islip];
                gdotmax = fmax(evolVals[islip], gdotmax);
                
                // smoothing factor to prevent numerical instabilities in the solve
               //  double h0 = 1.5*m_hdn_init;
               //  double k = 15.0/m_hdn_init;
                kfact[islip] = 1.0;//1.0/(1.0+exp(-k*(exp(h[islip])-h0)));
            }
            
            double frel[SlipGeom::nslip] = { 0.0 };
            if (gdotmax > 1e-10) {
                for (int islip = 0; islip < SlipGeom::nslip; islip++) {
                    double ratio = evolVals[islip] / gdotmax;
                    double ratiothres = 0.01;
                    frel[islip] = 1.0-1.0/(1.0+exp(-100.0*(ratio-ratiothres)));
                    
                    //if (h[islip] < log(1e2*m_hdn_min)) frel[islip] = 0.0;
                    //if (h[islip] < log(2*m_hdn_init)) frel[islip] = 0.0;
                    
                    frel[islip] = frel[islip] * kfact[islip];
                }
            }
            
            // Define k1 as a function of the orientation
            double k1[SlipGeom::nslip];
            //printf("getSdotN:\n");
            for (int islip = 0; islip < SlipGeom::nslip; islip++) {
                k1[islip] = m_k1;
                if (SlipGeom::dynamic) {
                    double chia = hvals[islip];
                    // If we are in the AT zone, then we need to increase k1
                    // to account for the fact that dislocations do take
                    // a longer path and thus are likely to multiply more
                    
                    double amin = 0.95;
					
                    //double a = fmin(1.0 + (amin - 1.0) * chia * 6.0 / M_PI, 1.0);
                    //a = 1.0 / a;
					
                    double a = 1.0/(1.0/cos(M_PI/6.0 - m_alpha_p)-1.0) * (1.0/amin - 1.0);
                    a = 1.0 + a * (1.0 / cos(chia - m_alpha_p) - 1.0);
					
                    k1[islip] = m_k1 * a;
                    //printf("sys[%d] chi = %e, rho = %e, gdot = %e\n",islip,chia*180.0/M_PI,exp(h[islip]),evolVals[islip]);
                }
            }
            
            // Define k2 as a function of gdot and temp_k
            double k2_ref = m_k2; // reference k2 value for 2e8/s at 300K
            double k2_temp = 0.05756349443979855 * log(temp_k / 7.309541735840538e-06);
            //printf("temp = %e, k2_temp = %e\n",temp_k,k2_temp);
            
            double rate_cut = 1e4; //1e-3;
            double lograte = log(0.5 * gtot * 1e6 + rate_cut);
            double k2_rate = -0.3433061910379516 * lograte + 7.586381434140954;
            //double k2_rate = 6.75092510e-03 * lograte * lograte - 5.72927375e-01 * lograte + 9.50718982e+00;
            
            
            double k2[SlipGeom::nslip];
            for (int islip = 0; islip < SlipGeom::nslip; islip++) {
                k2[islip] = k2_ref * k2_rate * k2_temp * kfact[islip];
                if (gtot < 1e-10) k2[islip] = 0.0;
                //if (h[islip] < log(1e2*m_hdn_min)) k2[islip] = 0.0;
                //printf("sys[%d] k1 = %e, k2 = %e, frel = %e\n",islip,k1[islip],k2[islip],frel[islip]);
            }
            //printf("gtot = %e, k2 = %e\n",gtot,k2);
            
            
            // TEST: adjust values while keeping the same saturation ratio k1/k2
            
            // for (int islip = 0; islip < SlipGeom::nslip; islip++) {
            //     //double r = 0.2;
            //     //k1[islip] *= r;
            //     //k2[islip] *= r;
				// //printf("sys %d: chi = %e, k1 = %e, k2 = %e, rho = %e\n",islip,hvals[islip]*180.0/M_PI,k1[islip],k2[islip],exp(h[islip])*1e4);
            // }
			//printf("---\n");
            
            
            // h = log(DD)
            // dDD / dt = DD * dh / dt
            // dh / dt = dDD / dt * 1 / DD
            // d DD_i / dt = (k1 * sqrt(A_{ij} DD_j) - k2 * DD_i) * gammadot_i
            // dh / dt = (k1 * sqrt(A_{ij} DD_j) / DD_i - k2) * gammadot_i
            // specialized here for the A_{ij} = I
            // dh / dt = (k1 / sqrt(DD_i) - k2) * gammadot_i
            // specialized case
            // d\dot{h} / dh = -1/2 * k_1 * (DD)^{-1/2}
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
               double temp1 = k1[islip] * temp_hs_a - k2[islip];
               sdot[islip] = temp1 * evolVals[islip] - frel[islip] * m_krelax;
               dsdot_ds[ECMECH_NN_INDX(islip, islip, SlipGeom::nslip)] = (-k1[islip] * onehalf * temp_hs_a) * evolVals[islip];
            }
         }
   }; // class KineticsBCCMD
} // namespace ecmech

#endif // ECMECH_KINETICS_BCCMD_H
